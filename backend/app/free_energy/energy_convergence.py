"""Energy-side convergence diagnostics for MM-GBSA ΔG_bind series.

This is the *energy* analogue of the structural C1 convergence check in
``scripts/corrected_fe_gate.py`` (``evaluate_c1``: pass iff the eq-window
begin-vs-end mean pocket-RMSD differ by < 0.4 Å). Where C1 asks "once the
ligand entered its settled window, did the *structure* stop drifting?", this
module asks the same of the *binding energy*.

It operates directly on the per-frame ``delta_g_total`` series that the
MM-GBSA estimator already restricts to the last-30 % production window
(``backend/app/free_energy/mmgbsa_runner.py`` — ``eq_window_last_fraction =
0.30``). That series is the equilibrated/production window; **no further
windowing is applied here**, exactly as C1 evaluates over the already-sliced
eq window. At the default 5 ps snapshot cadence a 1 ns run yields ~60 windowed
frames, so the autocorrelation-based statistical inefficiency below is
necessarily a coarse estimate (see ``analyze_run`` notes).

All public functions are pure (numpy + stdlib only) and side-effect free so
they can be unit-tested hermetically. File/manifest I/O and the
cross-replicate calibration live beside the structural gate in
``scripts/energy_convergence_gate.py``.
"""

from __future__ import annotations

import dataclasses
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# === Energy-convergence thresholds (parallel to corrected_fe_gate.py) ===
# The eq-window begin/end slice fraction is shared with the structural C1
# gate (``corrected_fe_gate.EQ_BEGIN_END_FRAC = 0.10``) so the energy twin
# reads the same first/last 10 % of the (already-windowed) series.
EQ_BEGIN_END_FRAC = 0.10

# Linear-drift flatness: a run is "flat" iff |slope| <= SLOPE_FLAT_SE_MULT *
# SE(slope) — i.e. the fitted drift is statistically indistinguishable from
# zero at ~2σ.
SLOPE_FLAT_SE_MULT = 2.0

# Snapshot cadence default — mirrors md.job.MdRunSettings.snapshot_every_ps.
DEFAULT_SNAPSHOT_PS = 5.0

# n_eff floor below which the corrected SEM is flagged unreliable.
LOW_NEFF_FRAMES = 10.0
# Series shorter than this make the autocorrelation/g estimate coarse.
COARSE_G_FRAMES = 80

# Calibrated begin-to-end energy-flatness bound, kcal/mol. The AUTHORITATIVE
# value is the calibration artifact written by
# ``scripts/energy_convergence_gate.py`` and re-loaded via
# ``load_energy_flatness_threshold()``; this constant is only the documented
# fallback used when that artifact is absent (e.g. a fresh clone). It is a
# hand-kept mirror of the persisted ``energy_flatness_threshold_kcal_per_mol``
# (the same way ``corrected_fe_gate``'s constants are echoed into its report),
# calibrated as mean(drift)+3·SD over the 23-replicate β-tubulin/taxol set.
# NOTE: the structural gate stores its bound as a hardcoded literal and only
# *echoes* it into output; this threshold goes further and is genuinely
# load-from-data, because it is calibration-derived rather than asserted.
# Mirror of docs/energy_convergence_calibration.json →
# calibration.energy_flatness_threshold_kcal_per_mol (mean+3SD over all 23
# runs). Hand-synced; the persisted artifact remains authoritative.
ENERGY_FLATNESS_FALLBACK_KCAL: Optional[float] = 12.376

# Default location of the persisted calibration artifact (repo-root/docs/).
CALIBRATION_ARTIFACT = (
    Path(__file__).resolve().parents[3] / "docs" / "energy_convergence_calibration.json"
)


# ---------------------------------------------------------------------------
# Primitive estimators (pure)
# ---------------------------------------------------------------------------
def _slice_count(n: int, frac: float = EQ_BEGIN_END_FRAC) -> int:
    """First/last slice size — mirrors ``corrected_fe_gate.py`` line 187
    (``max(1, int(round(eq_count * EQ_BEGIN_END_FRAC)))``)."""
    return max(1, int(round(n * frac)))


def begin_end_drift(series: Sequence[float], frac: float = EQ_BEGIN_END_FRAC) -> float:
    """``|mean(first frac) − mean(last frac)|`` of the windowed ΔG series.

    Energy twin of ``corrected_fe_gate.evaluate_c1``'s eq_begin/eq_end on
    pocket-RMSD. Returns kcal/mol (always ≥ 0); ``nan`` for an empty series.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n == 0:
        return float("nan")
    k = _slice_count(n, frac)
    return abs(float(x[:k].mean()) - float(x[-k:].mean()))


def linear_drift(
    series: Sequence[float], dt_ps: float = DEFAULT_SNAPSHOT_PS
) -> Tuple[float, float, Optional[bool]]:
    """OLS of ΔG vs time. Returns ``(slope, slope_se, flat)`` where ``slope``
    is kcal/mol per ps, ``slope_se`` its standard error, and ``flat`` is True
    iff ``|slope| <= SLOPE_FLAT_SE_MULT * slope_se``. ``flat`` is ``None`` when
    N<3 (the residual variance / SE is undefined).

    The absolute time offset of the window does not affect the slope or its
    SE, so feeding the window-local frame index (×dt) is equivalent to the
    true trajectory time.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 3:
        return (float("nan"), float("nan"), None)
    t = np.arange(n, dtype=float) * float(dt_ps)
    tbar = t.mean()
    sxx = float(((t - tbar) ** 2).sum())
    if sxx == 0.0:
        return (float("nan"), float("nan"), None)
    xbar = x.mean()
    slope = float(((t - tbar) * (x - xbar)).sum() / sxx)
    intercept = xbar - slope * tbar
    resid = x - (intercept + slope * t)
    sse = float((resid ** 2).sum())
    sigma2 = sse / (n - 2)
    slope_se = math.sqrt(sigma2 / sxx) if sigma2 > 0 else 0.0
    flat = abs(slope) <= SLOPE_FLAT_SE_MULT * slope_se
    return (slope, slope_se, bool(flat))


def _autocorrelation_g(series: Sequence[float]) -> float:
    """Statistical inefficiency ``g = 1 + 2 Σ_t (1 − t/N) C(t)`` truncated at
    the first non-positive autocorrelation after lag 1 — the pure-numpy
    fallback used when pymbar is unavailable.

    ``C(t)`` is the normalized autocorrelation at lag t using the biased
    (population) variance, matching pymbar's convention. Returns ``g ≥ 1.0``.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 2:
        return 1.0
    dx = x - x.mean()
    var = float((dx * dx).mean())  # biased variance, pymbar convention
    if var <= 0.0:
        return 1.0  # constant series → no autocorrelation penalty
    g = 1.0
    for t in range(1, n):
        c_t = float((dx[: n - t] * dx[t:]).sum() / (n - t)) / var
        # "first non-positive autocorrelation after lag 1": lag 1 is always
        # included; from lag 2 on, stop at the first C(t) <= 0 (before adding).
        if c_t <= 0.0 and t > 1:
            break
        g += 2.0 * (1.0 - t / n) * c_t
    return max(1.0, g)


def statistical_inefficiency(series: Sequence[float]) -> Tuple[float, str]:
    """Return ``(g, source)``. Prefers ``pymbar.timeseries`` if importable,
    else the in-house ``_autocorrelation_g``. ``source`` is one of
    ``"pymbar"`` / ``"fallback"`` / ``"trivial"`` for provenance."""
    x = np.asarray(series, dtype=float)
    if x.size < 2:
        return 1.0, "trivial"
    try:
        from pymbar import timeseries as _ts  # type: ignore

        fn = getattr(_ts, "statistical_inefficiency", None) or getattr(
            _ts, "statisticalInefficiency"
        )
        g = float(fn(x))
        return max(1.0, g), "pymbar"
    except Exception:
        return _autocorrelation_g(x), "fallback"


def corrected_sem(series: Sequence[float]) -> Tuple[float, float, float, str]:
    """Return ``(g, n_eff, sem, g_source)`` with ``n_eff = N/g`` and
    ``sem = std(ddof=1) / sqrt(n_eff)`` (autocorrelation-corrected SEM)."""
    x = np.asarray(series, dtype=float)
    n = x.size
    g, source = statistical_inefficiency(x)
    g = max(1.0, float(g))
    n_eff = (n / g) if g > 0 else float(n)
    if n < 2:
        return g, n_eff, float("nan"), source
    sd = float(np.std(x, ddof=1))
    sem = sd / math.sqrt(n_eff) if n_eff > 0 else float("nan")
    return g, n_eff, sem, source


# ---------------------------------------------------------------------------
# Equilibration-onset detection & truncation helpers (pure)
# ---------------------------------------------------------------------------
def cumulative_mean(series: Sequence[float]) -> List[float]:
    """Running (cumulative) mean from t=0: element i is mean(series[:i+1])."""
    x = np.asarray(series, dtype=float)
    if x.size == 0:
        return []
    return (np.cumsum(x) / np.arange(1, x.size + 1)).tolist()


def last_fraction_window(
    series: Sequence[float], last_fraction: float = 0.30
) -> Tuple[int, List[float]]:
    """Return ``(start_index, windowed_slice)`` for the last ``last_fraction``
    of the series — the SAME convention mmgbsa_runner uses for the production
    window (``eq_window_start = int(N * 0.7)``)."""
    x = list(series)
    n = len(x)
    start = int(n * (1.0 - last_fraction))
    return start, x[start:]


def detect_equilibration(
    series: Sequence[float], nskip: int = 1
) -> Tuple[int, float, float]:
    """Locate the equilibration onset t0 (a frame index) that MAXIMIZES the
    number of post-t0 independent samples ``N_eff(t0) = (N − t0)/g(series[t0:])``
    — i.e. the running-mean plateau onset, mirroring
    ``pymbar.timeseries.detectEquilibration``. Prefers pymbar if importable,
    else an in-house O(N²·t_cut) scan. Returns ``(t0, g_at_t0, n_eff_max)``.

    Applied to either signal (pocket-aligned RMSD trend or per-frame ΔG), this
    detects WHERE the running mean stabilizes; for RMSD the caller uses only
    the onset (trend), never the absolute distance from the docked pose, so
    relaxation is excluded rather than mistaken for non-equilibrium.
    """
    x = np.asarray(series, dtype=float)
    n = int(x.size)
    if n < 3:
        return 0, 1.0, float(n)
    try:
        from pymbar import timeseries as _ts  # type: ignore

        fn = getattr(_ts, "detect_equilibration", None) or getattr(
            _ts, "detectEquilibration"
        )
        t0, g, neff = fn(x)
        return int(t0), float(g), float(neff)
    except Exception:
        pass
    best_t0, best_neff, best_g = 0, -1.0, 1.0
    step = max(1, int(nskip))
    # Leave at least 3 samples in the tail so g is estimable.
    for t0 in range(0, n - 2, step):
        g = _autocorrelation_g(x[t0:])
        neff = (n - t0) / g
        if neff > best_neff:
            best_neff, best_t0, best_g = neff, t0, g
    return best_t0, best_g, best_neff


# ---------------------------------------------------------------------------
# Per-run analysis
# ---------------------------------------------------------------------------
@dataclass
class RunConvergence:
    """Per-run energy-convergence diagnostics. The fields the user asked be
    added to the result-artifact schema are: ``begin_end_drift``, ``slope``,
    ``slope_se``, ``slope_flat``, ``g``, ``n_eff``, ``corrected_sem``."""

    n_frames: int
    dt_ps: float
    window_mean: float
    begin_end_drift: float
    slope: float
    slope_se: float
    slope_flat: Optional[bool]
    g: float
    n_eff: float
    corrected_sem: float
    g_source: str
    low_neff: bool
    stored_mean: Optional[float] = None
    mean_mismatch: Optional[float] = None
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


def analyze_run(
    series: Sequence[float],
    dt_ps: float = DEFAULT_SNAPSHOT_PS,
    stored_mean: Optional[float] = None,
) -> RunConvergence:
    """Full per-run convergence diagnostics over the already-windowed ΔG
    series. ``stored_mean`` (e.g. ``result.delta_g_mean_kcal_per_mol``) is
    cross-checked against the recomputed window mean when supplied."""
    x = np.asarray(series, dtype=float)
    n = int(x.size)
    notes: List[str] = []

    if n == 0:
        notes.append("empty ΔG series.")
        return RunConvergence(
            n_frames=0, dt_ps=float(dt_ps), window_mean=float("nan"),
            begin_end_drift=float("nan"), slope=float("nan"),
            slope_se=float("nan"), slope_flat=None, g=float("nan"),
            n_eff=float("nan"), corrected_sem=float("nan"),
            g_source="trivial", low_neff=True, stored_mean=stored_mean,
            mean_mismatch=None, notes=notes,
        )

    wmean = float(x.mean())
    bed = begin_end_drift(x)
    slope, slope_se, flat = linear_drift(x, dt_ps)
    g, n_eff, sem, g_source = corrected_sem(x)
    low_neff = bool(n_eff < LOW_NEFF_FRAMES)

    if n <= COARSE_G_FRAMES:
        notes.append(
            f"N={n} frames: autocorrelation/statistical-inefficiency estimate "
            f"is coarse (short series)."
        )
    if low_neff:
        notes.append(
            f"n_eff={n_eff:.1f} < {LOW_NEFF_FRAMES:.0f}: corrected SEM is "
            f"unreliable for this run."
        )

    mean_mismatch: Optional[float] = None
    if stored_mean is not None:
        mean_mismatch = abs(wmean - float(stored_mean))
        tol = 1e-6 * max(1.0, abs(float(stored_mean))) + 1e-6
        if mean_mismatch > tol:
            notes.append(
                f"window mean {wmean:.4f} != stored result mean "
                f"{float(stored_mean):.4f} (Δ={mean_mismatch:.2e}); the "
                f"per-frame series may not match the stored aggregate."
            )

    return RunConvergence(
        n_frames=n, dt_ps=float(dt_ps), window_mean=wmean,
        begin_end_drift=bed, slope=slope, slope_se=slope_se, slope_flat=flat,
        g=g, n_eff=n_eff, corrected_sem=sem, g_source=g_source,
        low_neff=low_neff, stored_mean=stored_mean,
        mean_mismatch=mean_mismatch, notes=notes,
    )


# ---------------------------------------------------------------------------
# Cross-replicate aggregation / calibration (pure)
# ---------------------------------------------------------------------------
def calibrate_flatness_threshold(
    drift_by_key: Dict[str, float], k_sd: float = 3.0
) -> Dict[str, Any]:
    """Calibrate the begin-to-end energy-flatness bound as
    ``mean(drift) + k_sd·SD(drift)`` over the replicate set, reported both
    with and without the single largest-drift outlier removed (mirroring how
    the structural drift bound was calibrated). ``drift_by_key`` maps a
    replicate key → its ``begin_end_drift``."""
    keys = list(drift_by_key.keys())
    drifts = np.array([float(drift_by_key[k]) for k in keys], dtype=float)
    n = int(drifts.size)
    if n == 0:
        raise ValueError("no drift values to calibrate from")

    mean = float(drifts.mean())
    sd = float(drifts.std(ddof=1)) if n >= 2 else float("nan")
    threshold = mean + k_sd * sd if math.isfinite(sd) else float("nan")

    imax = int(np.argmax(drifts))
    outlier_key = keys[imax]
    outlier_drift = float(drifts[imax])

    mask = np.ones(n, dtype=bool)
    mask[imax] = False
    d2 = drifts[mask]
    if d2.size >= 2:
        mean2 = float(d2.mean())
        sd2 = float(d2.std(ddof=1))
        threshold2 = mean2 + k_sd * sd2
    else:
        mean2 = float(d2.mean()) if d2.size else float("nan")
        sd2 = float("nan")
        threshold2 = float("nan")

    order = sorted(range(n), key=lambda i: drifts[i])
    distribution = [
        {"key": keys[i], "begin_end_drift": float(drifts[i])} for i in order
    ]

    return {
        "k_sd": float(k_sd),
        "n_runs": n,
        "mean_drift": mean,
        "sd_drift": sd,
        "energy_flatness_threshold_kcal_per_mol": threshold,
        "gross_outlier_key": outlier_key,
        "gross_outlier_drift": outlier_drift,
        "mean_drift_no_outlier": mean2,
        "sd_drift_no_outlier": sd2,
        "energy_flatness_threshold_no_outlier_kcal_per_mol": threshold2,
        "drift_distribution": distribution,
    }


def between_replicate_stats(
    run_means_by_pose: Dict[int, Sequence[float]],
) -> Tuple[Dict[int, Dict[str, float]], float]:
    """Per-pose between-replicate ``{n, mean, sd, sem}`` of run-mean ΔG
    (sd ddof=1, sem = sd/√n) plus the pooled SD over all runs."""
    out: Dict[int, Dict[str, float]] = {}
    allvals: List[float] = []
    for pose in sorted(run_means_by_pose.keys()):
        v = np.asarray(list(run_means_by_pose[pose]), dtype=float)
        nn = int(v.size)
        mean = float(v.mean()) if nn else float("nan")
        sd = float(v.std(ddof=1)) if nn >= 2 else float("nan")
        sem = sd / math.sqrt(nn) if (nn >= 2) else float("nan")
        out[pose] = {"n": nn, "mean": mean, "sd": sd, "sem": sem}
        allvals.extend(v.tolist())
    a = np.asarray(allvals, dtype=float)
    pooled_sd = float(a.std(ddof=1)) if a.size >= 2 else float("nan")
    return out, pooled_sd


def within_vs_between_verdict(
    median_within_sem: float, between_sd: float
) -> Tuple[str, float]:
    """One-line verdict on the limiting uncertainty. Returns
    ``(verdict, ratio)`` with ``ratio = between_sd / median_within_sem``."""
    if not (math.isfinite(median_within_sem) and math.isfinite(between_sd)) \
            or median_within_sem <= 0:
        return ("indeterminate (non-finite inputs)", float("nan"))
    ratio = between_sd / median_within_sem
    if ratio >= 2.0:
        verdict = (
            "limiting uncertainty is SLOW between-seed sampling "
            "(between-replicate SD ≫ within-run SEM): add replicates, not length"
        )
    elif ratio <= 0.5:
        verdict = (
            "limiting uncertainty is FAST within-run sampling "
            "(within-run SEM ≫ between-replicate SD): extend trajectories"
        )
    else:
        verdict = (
            "within-run and between-seed contributions are comparable "
            "(neither dominates)"
        )
    return (verdict, ratio)


# ---------------------------------------------------------------------------
# Threshold loader + gate (the "as data" path)
# ---------------------------------------------------------------------------
def load_energy_flatness_threshold(
    path: Path = CALIBRATION_ARTIFACT, *, use_outlier_excluded: bool = False
) -> Optional[float]:
    """Read the calibrated begin-to-end energy-flatness bound from the
    persisted calibration artifact (the authoritative, load-from-data
    source). Falls back to ``ENERGY_FLATNESS_FALLBACK_KCAL`` when the artifact
    is missing/unreadable."""
    key = (
        "energy_flatness_threshold_no_outlier_kcal_per_mol"
        if use_outlier_excluded
        else "energy_flatness_threshold_kcal_per_mol"
    )
    try:
        d = json.loads(Path(path).read_text())
        val = d.get("calibration", d).get(key)
        if val is not None and math.isfinite(float(val)):
            return float(val)
    except Exception:
        pass
    return ENERGY_FLATNESS_FALLBACK_KCAL


# NOTE (Stage 2): the energy-flatness GATE (`evaluate_energy_convergence`, the
# begin-end-drift < threshold PASS/FAIL) has been REMOVED. Energy no longer
# qualifies or declines a run — the gate is the STRUCTURAL verdict
# (md.analyze.classify_stability). `begin_end_drift` above survives as a reported
# diagnostic only. `load_energy_flatness_threshold` / `calibrate_flatness_threshold`
# remain as calibration utilities (scripts/energy_convergence_gate.py) but feed
# no gate.


# ---------------------------------------------------------------------------
# Running-mean landing / drift-noise test (Stage 1b primitives)
# ---------------------------------------------------------------------------
# These generic time-series primitives live here (next to corrected_sem /
# statistical_inefficiency, which they use) and are imported by BOTH the
# structural gate (md.analyze: pose-RMSD landing, k_sigma + Å effect-size floor)
# AND the energy path (t0_energy on the ΔG series, kcal/mol effect-size floor).
# The effect-size floor's absolute minimum is a caller-supplied parameter so the
# same test serves either unit; the defaults below reproduce the structural Å
# values so md.analyze's behavior is unchanged by the relocation.
DRIFT_NOISE_K_SIGMA = 2.0          # |Δ|/SE(Δ) ≤ this ⇒ statistically flat
LANDING_MIN_WINDOW_FRAMES = 20     # shortest window the drift test will judge
LANDING_END_GRID = 12              # sub-window end-points sampled per candidate
EFFECT_SIZE_FRAC = 0.25            # fraction of within-window stddev
EFFECT_SIZE_ABS_MIN = 0.20         # default absolute floor (structural Å value)

# Energy-side absolute floor for the ΔG running-mean landing test (kcal/mol):
# a begin-to-end change in the windowed ΔG mean below this is "nobody cares"
# (negligible vs MM-GBSA SEMs of ~0.3-0.9 kcal/mol). Self-scaling means the
# 0.25·stddev term (≈1 kcal/mol on these ~4 kcal/mol-stddev series) usually
# dominates this floor; it is the small backstop for low-noise series.
ENERGY_EFFECT_SIZE_ABS_MIN_KCAL = 0.5


@dataclass
class DriftNoiseResult:
    converged: Optional[bool]   # None when the window is too short to judge
    delta: float                # mean(last frac) − mean(first frac)
    se: float                   # autocorrelation-corrected SE(Δ)
    z: float                    # |delta| / se
    g: float                    # statistical inefficiency over the window
    n: int
    effect_floor: float = float("nan")  # |Δ| at/below which the change is trivial


def drift_noise_test(
    window: Sequence[float],
    k_sigma: float = DRIFT_NOISE_K_SIGMA,
    frac: float = 0.10,
    *,
    effect_size_frac: float = EFFECT_SIZE_FRAC,
    effect_size_abs_min: float = EFFECT_SIZE_ABS_MIN,
) -> DriftNoiseResult:
    """Drift/noise test for a signal window:

        Δ      = mean(last `frac`) − mean(first `frac`)
        SE(Δ)  = sqrt(SE(mean_first)² + SE(mean_last)²)   [autocorr-corrected]
        converged  iff  (|Δ|/SE ≤ k_sigma)  OR  (|Δ| ≤ effect_floor)

    The end-slice (not whole-window) variance makes it a true DRIFT detector;
    the effect-size floor `max(effect_size_abs_min, effect_size_frac·stddev)`
    keeps a chemically/energetically trivial change from being flagged as drift
    even on a low-noise signal (where SE is tiny). Trend-safe: for a linear
    climb Δ≈3.46·stddev ≫ 0.25·stddev, so the floor cannot mask a real climb.
    NaNs dropped first; converged=None when < 4 finite samples remain."""
    x = np.asarray([v for v in window
                    if v is not None and math.isfinite(v)], dtype=float)
    n = int(x.size)
    if n < 4:
        return DriftNoiseResult(None, float("nan"), float("nan"),
                                float("nan"), float("nan"), n, float("nan"))
    k = min(max(2, int(round(n * frac))), n // 2)
    first, last = x[:k], x[-k:]
    delta = float(last.mean() - first.mean())
    g_f, _, sem_f, _ = corrected_sem(first)
    g_l, _, sem_l, _ = corrected_sem(last)
    se = math.sqrt(
        (0.0 if math.isnan(sem_f) else sem_f) ** 2
        + (0.0 if math.isnan(sem_l) else sem_l) ** 2
    )
    if se == 0.0:
        z = 0.0 if delta == 0.0 else float("inf")
    else:
        z = abs(delta) / se
    scale = float(np.std(x, ddof=1)) if n > 1 else 0.0
    effect_floor = max(float(effect_size_abs_min), float(effect_size_frac) * scale)
    converged = bool((z <= k_sigma) or (abs(delta) <= effect_floor))
    return DriftNoiseResult(converged, delta, se, z,
                            max(float(g_f), float(g_l)), n, effect_floor)


def _landing_end_grid(s: int, m: int, min_window: int, grid: int) -> List[int]:
    """End-points e for the sub-windows [s:e] tested at candidate start s —
    a `grid`-point sample of [s+min_window, m], always including m."""
    lo = s + min_window
    if lo >= m:
        return [m]
    pts = {lo + int(round(i * (m - lo) / (grid - 1))) for i in range(grid)}
    pts.add(m)
    return sorted(e for e in pts if lo <= e <= m)


def detect_landing(
    series: Sequence[float],
    k_sigma: float = DRIFT_NOISE_K_SIGMA,
    *,
    min_window: int = LANDING_MIN_WINDOW_FRAMES,
    end_grid: int = LANDING_END_GRID,
    effect_size_frac: float = EFFECT_SIZE_FRAC,
    effect_size_abs_min: float = EFFECT_SIZE_ABS_MIN,
    tail_tolerance_frac: float = 0.0,
) -> Optional[int]:
    """Landing frame t0 = the earliest frame after which the running mean
    passes the drift/noise test for EVERY sub-window [t0:e] from there to the
    end (reaches AND holds — a transient flattening that later resumes drifting
    fails). Returns the frame index into the ORIGINAL series (NaN frames are
    dropped for the test but the returned index maps back). None when the run
    never lands. `effect_size_abs_min` is the floor in the SIGNAL's units
    (Å for pose-RMSD, kcal/mol for ΔG).

    `tail_tolerance_frac` (Stage 1c — ROBUST HOLDS) judges convergence over the
    BULK of the post-landing window only, dropping the terminal
    `tail_tolerance_frac · (m − s)` frames from every candidate start `s`. A
    brief terminal excursion (≤ that fraction) is therefore NOT counted as loss
    of plateau, while a SUSTAINED late climb necessarily extends earlier than
    the tail, lands inside the gated bulk, and still fails. `0.0` (default)
    recovers the exact reaches-and-holds-to-the-very-end behavior (used by the
    energy operative-window path, which must not tolerate any tail drift)."""
    finite = [(i, v) for i, v in enumerate(series)
              if v is not None and math.isfinite(v)]
    m = len(finite)
    if m < min_window:
        return None
    fidx = [i for i, _ in finite]
    fval = [v for _, v in finite]
    for s in range(0, m - min_window + 1):
        # Robust holds: cap the judged window at the bulk end-point, tolerating
        # the terminal excursion. tail_tolerance_frac=0 ⇒ e_top = m (exact).
        drop = (int(math.floor(tail_tolerance_frac * (m - s)))
                if tail_tolerance_frac > 0 else 0)
        e_top = m - drop
        if e_top - s < min_window:
            continue
        if all(drift_noise_test(fval[s:e], k_sigma,
                                effect_size_frac=effect_size_frac,
                                effect_size_abs_min=effect_size_abs_min).converged
               for e in _landing_end_grid(s, e_top, min_window, end_grid)):
            return fidx[s]
    return None


# ---------------------------------------------------------------------------
# Window-choice comparison + OPERATIVE window (Stage 2)
# ---------------------------------------------------------------------------
# Operates on the stored from-frame-0 per-frame ``delta_g_total`` series — no MD.
# The OPERATIVE window for the ΔG average is the energy running-mean landing
# window [t0_energy:end], where t0_energy is detected by the Stage-1b
# ``detect_landing`` (drift/noise test WITH the effect-size floor) applied to the
# ΔG series with a kcal/mol floor. The legacy ``last_30pct`` and the
# ``fixed_discard`` window are still reported for transparency, but only
# ``detected_t0_energy`` is operative. The old max(t0_RMSD, t0_energy) rule and
# the last-30% default are gone. Begin-end drift is reported as a DIAGNOSTIC,
# never a gate (energy no longer gates — structure does).
DEFAULT_FIXED_DISCARD_PS = 700.0


@dataclass
class WindowStat:
    """Window mean + autocorrelation-corrected uncertainty for one window
    choice over the from-frame-0 ΔG series."""

    label: str
    t0_frame: int
    t0_ps: float
    n_frames: int
    window_mean: float
    g: float
    n_eff: float
    corrected_sem: float
    g_source: str
    low_neff: bool

    def to_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


def _window_stat(label: str, series: Sequence[float], t0_frame: int,
                 dt_ps: float) -> WindowStat:
    x = np.asarray(series, dtype=float)
    t0_frame = max(0, int(t0_frame))
    sl = x[t0_frame:]
    t0_ps = t0_frame * float(dt_ps)
    if sl.size == 0:
        return WindowStat(label, t0_frame, t0_ps, 0, float("nan"), float("nan"),
                          float("nan"), float("nan"), "trivial", True)
    g, n_eff, sem, src = corrected_sem(sl)
    return WindowStat(label, t0_frame, t0_ps, int(sl.size), float(sl.mean()),
                      g, n_eff, sem, src, bool(n_eff < LOW_NEFF_FRAMES))


def compare_equilibration_windows(
    series: Sequence[float],
    dt_ps: float = DEFAULT_SNAPSHOT_PS,
    *,
    fixed_discard_ps: float = DEFAULT_FIXED_DISCARD_PS,
    effect_size_abs_min_kcal: float = ENERGY_EFFECT_SIZE_ABS_MIN_KCAL,
    n_eff_floor: float = LOW_NEFF_FRAMES,
) -> Dict[str, Any]:
    """Window mean + N_eff + corrected SEM under three window choices over the
    FULL from-frame-0 ΔG ``series``, plus the OPERATIVE window and its
    sampling-sufficiency stats.

    Windows reported (transparency): ``last_30pct``, ``fixed_discard``, and
    ``detected_t0_energy`` (= the energy running-mean landing window
    [t0_energy:end]). The OPERATIVE window is ``detected_t0_energy``.

    Sampling-sufficiency (the only surviving energy-side check, an ADEQUACY
    guard — not a convergence gate): the operative window should span
    N_eff ≥ ``n_eff_floor`` (≈ m·τ with m=10). ``low_neff`` flags shortfall.
    Begin-end drift is reported as a diagnostic only."""
    x = np.asarray(series, dtype=float)
    n = int(x.size)

    last30 = _window_stat("last_30pct", x, int(n * 0.7), dt_ps)
    disc_frames = max(0, int(round(float(fixed_discard_ps) / float(dt_ps))))
    fixed = _window_stat(
        f"fixed_discard_{int(round(fixed_discard_ps))}ps", x, disc_frames, dt_ps,
    )

    # OPERATIVE: energy running-mean landing (reaches & holds), kcal/mol floor.
    t0e = detect_landing(x, effect_size_abs_min=effect_size_abs_min_kcal)
    energy_landed = t0e is not None
    op_start = t0e if energy_landed else 0
    detected = _window_stat("detected_t0_energy", x, op_start, dt_ps)

    # Sampling-sufficiency over the operative window. τ_int = (g−1)/2 frames
    # (g = 1 + 2τ_int); N_eff = N/g. low_neff iff N_eff < floor (≈ m·τ, m=10).
    g = detected.g
    tau_int_frames = (g - 1.0) / 2.0 if math.isfinite(g) else float("nan")
    operative = {
        "window": "detected_t0_energy",
        "energy_landed": energy_landed,
        "t0_energy_frame": (int(t0e) if energy_landed else None),
        "t0_energy_ps": (int(t0e) * float(dt_ps) if energy_landed else None),
        "effect_size_abs_min_kcal": float(effect_size_abs_min_kcal),
        "n_frames": detected.n_frames,
        "delta_g_mean": detected.window_mean,
        "g": g,
        "g_source": detected.g_source,
        "tau_int_frames": tau_int_frames,
        "tau_int_ps": (tau_int_frames * float(dt_ps)
                       if math.isfinite(tau_int_frames) else float("nan")),
        "n_eff": detected.n_eff,
        "corrected_sem": detected.corrected_sem,
        "n_eff_floor": float(n_eff_floor),
        "low_neff": bool(detected.low_neff),
        "note": ("ADEQUACY guard (N_eff ≥ floor), NOT a convergence gate; "
                 "the gate is the structural verdict. ΔG averaged over "
                 "[t0_energy:end]; full series if energy never landed."),
    }

    return {
        "n_frames_total": n,
        "dt_ps": float(dt_ps),
        "fixed_discard_ps": float(fixed_discard_ps),
        "operative_window": "detected_t0_energy",
        "windows": {
            "last_30pct": last30.to_dict(),
            "fixed_discard": fixed.to_dict(),
            "detected_t0_energy": detected.to_dict(),
        },
        "operative": operative,
        # DIAGNOSTIC ONLY (not a gate): begin-end drift over the operative window.
        "begin_end_drift_diagnostic_kcal": begin_end_drift(x[op_start:]),
    }
