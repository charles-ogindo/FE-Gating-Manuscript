"""Recompute an MD run's MM-GBSA ΔG_bind under the 3 ns windowing convention.

The legacy ``compute_md_fe`` averages ΔG over the last 30 % of frames. The 3 ns
extension runs instead discard a fixed equilibration period (default 0.7 ns) and
average all later frames (~2.3 ns on a full 3 ns run). To make the window choice
auditable we evaluate the FULL from-frame-0 per-frame ΔG series ONCE (MM-GBSA
single-points over every frame — this re-runs the *estimator*, NOT the MD) and
then report the window mean + N_eff + autocorrelation-corrected SEM under three
choices:

  (i)   last_30pct      — legacy.
  (ii)  fixed_discard   — discard first ``--discard-ps`` (default 700), average rest.
  (iii) detected_t0_run — t0_run = max(t0_RMSD, t0_energy); RMSD onset trend-only.

The full from-0 series + the three windows are persisted to
``jobs/<md_id>/free_energy/windows.json`` so a later gate/aggregation reads them
back. Reusable for every extension run going forward:

    python scripts/fe_windows.py <md_id> [<md_id> ...] [--discard-ps 700]

Defaults to pose0_rep2 (8b70b23b) when no md_id is given.
"""

from __future__ import annotations

# CRITICAL preload order (mirror run_length_convergence.py): libexpat
# (RTLD_GLOBAL) then openmm, BOTH before any backend import.
import ctypes as _ctypes
import os as _os

_libexpat = _os.path.join(_os.environ.get("CONDA_PREFIX", ""), "lib", "libexpat.so.1")
if _os.path.exists(_libexpat):
    _ctypes.CDLL(_libexpat, mode=_ctypes.RTLD_GLOBAL)
import openmm as _preload_openmm  # noqa: F401,E402

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Dict, List, Optional  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR  # noqa: E402
from backend.app.free_energy.energy_convergence import (  # noqa: E402
    compare_equilibration_windows,
    DEFAULT_FIXED_DISCARD_PS,
)
from backend.app.free_energy.gating import validate  # noqa: E402
from backend.app.free_energy.mmgbsa import estimate_mmgbsa  # noqa: E402
from backend.app.free_energy.mmgbsa_runner import (  # noqa: E402
    _FF_IMPLICIT,
    _FF_PROTEIN_AND_IONS,
    _FF_SMALL_MOL,
    _collect_frame_positions,
    prepare_mmgbsa_systems,
)

DEFAULT_IDS = ["8b70b23b-b0ed-44c6-ab2a-9e1c7891c0f8"]  # pose0_rep2


def _stored_series(md_id: str) -> Optional[Dict[str, object]]:
    """Reuse a previously-written from-frame-0 ΔG series (no MM-GBSA re-eval).
    The task operates on the STORED series; recomputation is only needed if it
    is absent."""
    p = JOBS_DIR / md_id / "free_energy" / "windows.json"
    if not p.is_file():
        return None
    try:
        d = json.loads(p.read_text())
        series = d.get("per_frame_delta_g_total_from_frame0")
        if not series:
            return None
        md_summary = json.loads((JOBS_DIR / md_id / "md" / "summary.json").read_text())
        return {
            "series": [float(v) for v in series],
            "dt_ps": float(d.get("dt_ps") or 5.0),
            "verdict": md_summary.get("verdict"),
            "n_frames_skipped": 0,
            "source": "stored windows.json (reused; no MM-GBSA re-eval)",
        }
    except Exception:
        return None


def full_delta_g_series(md_id: str) -> Dict[str, object]:
    """Evaluate MM-GBSA over EVERY frame → the from-frame-0 ΔG_total series."""
    md_dir = JOBS_DIR / md_id / "md"
    summary = json.loads((md_dir / "summary.json").read_text())
    md_settings = summary.get("settings") or {}
    solvent = md_settings.get("solvent", "implicit")

    (cpx_sys, rec_sys, lig_sys, cpx_top, rec_top, lig_top,
     ridx, lidx) = prepare_mmgbsa_systems(md_id, use_cache=True)

    positions = _collect_frame_positions(
        md_dir, expected_topology=cpx_top, solvent_mode=solvent,
    )  # ALL frames, from frame 0 — no last-30% slice, no head discard.

    gate = validate(md_id).to_dict()
    method_meta = {
        "name": "single-trajectory MM-GBSA (full from-frame-0 series)",
        "implicit_solvent": "OBC2",
        "force_fields": _FF_PROTEIN_AND_IONS + [_FF_IMPLICIT],
        "small_molecule_forcefield": _FF_SMALL_MOL,
        "nonbonded_method": "CutoffNonPeriodic",
        "nonbonded_cutoff_nm": 1.0,
        "md_trajectory_solvent": solvent,
    }
    result = estimate_mmgbsa(
        complex_system=cpx_sys, complex_topology=cpx_top,
        receptor_system=rec_sys, receptor_topology=rec_top,
        ligand_system=lig_sys, ligand_topology=lig_top,
        frame_positions=positions,
        receptor_idx=ridx, ligand_idx=lidx,
        gate_can_run=gate["can_run"],
        gate_reason=(gate["reasons"][0]["message"] if gate["reasons"] else None),
        method_meta=method_meta,
    )
    series = [fe.delta_g_total for fe in result.per_frame]
    return {
        "series": series,
        "dt_ps": float(md_settings.get("snapshot_every_ps") or 5.0),
        "verdict": summary.get("verdict"),
        "n_frames_skipped": result.n_frames_skipped,
        "source": "recomputed (MM-GBSA single-points over all frames)",
    }


def _row(label: str, ws: Dict[str, object], operative: bool = False) -> str:
    star = " ◀ OPERATIVE" if operative else ""
    return (f"  {label:<24}{ws['t0_ps']/1000:>7.2f}{ws['n_frames']:>7}"
            f"{ws['window_mean']:>11.3f}{ws['n_eff']:>9.1f}"
            f"{ws['corrected_sem']:>9.3f}{('  LOW Neff' if ws['low_neff'] else '')}{star}")


def run_one(md_id: str, discard_ps: float, *, reuse: bool = True) -> None:
    print(f"\n{'='*74}\nFE windows — {md_id[:8]}\n{'='*74}")
    t0 = time.perf_counter()
    data = _stored_series(md_id) if reuse else None
    if data is None:
        data = full_delta_g_series(md_id)
    wall = time.perf_counter() - t0
    series = data["series"]
    dt = data["dt_ps"]

    win = compare_equilibration_windows(series, dt_ps=dt, fixed_discard_ps=discard_ps)
    op = win["operative"]
    total_ps = (len(series) - 1) * dt
    drifting = (data["verdict"] == "drifting")
    short = total_ps < 2900  # < ~2.9 ns → did not reach the full 3 ns

    print(f"  series: {len(series)} frames from frame 0; trajectory "
          f"{total_ps/1000:.3f} ns  [{data['source']}]  ({wall:.0f} s)")
    print(f"  GATE = structural verdict (md.classify_stability): "
          f"{data['verdict']}  (energy no longer gates)")
    print(f"  t0_energy landing: "
          + (f"frame {op['t0_energy_frame']} = {op['t0_energy_ps']/1000:.3f} ns "
             f"(abs_min={op['effect_size_abs_min_kcal']} kcal/mol)"
             if op["energy_landed"]
             else "NONE — energy running mean never lands; operative window = full series"))
    print(f"\n  {'window':<24}{'t0(ns)':>7}{'N':>7}{'mean':>11}{'Neff':>9}{'cSEM':>9}")
    print(f"  {'-'*67}")
    print(_row("(i) last-30%", win["windows"]["last_30pct"]))
    print(_row(f"(ii) fixed {discard_ps:.0f}ps disc", win["windows"]["fixed_discard"]))
    print(_row("(iii) detected t0_energy", win["windows"]["detected_t0_energy"],
               operative=True))
    print(f"\n  OPERATIVE (t0_energy→end): ΔG = {op['delta_g_mean']:.3f} ± "
          f"{op['corrected_sem']:.3f} kcal/mol  over N={op['n_frames']} frames")
    print(f"    g={op['g']:.2f}  τ_int={op['tau_int_ps']:.0f} ps  "
          f"N_eff={op['n_eff']:.1f}  (floor {op['n_eff_floor']:.0f}; "
          f"{'LOW_NEFF — under-sampled' if op['low_neff'] else 'adequate'})")
    print(f"    begin-end drift (DIAGNOSTIC only, not a gate): "
          f"{win['begin_end_drift_diagnostic_kcal']:.3f} kcal/mol")
    if drifting or short:
        flags = ([f"structural verdict=drifting"] if drifting else []) + \
                ([f"only {total_ps/1000:.3f} ns (< 3 ns)"] if short else [])
        print(f"  ⚠ PRELIMINARY: {', '.join(flags)}.")

    out_dir = JOBS_DIR / md_id / "free_energy"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "md_job_id": md_id,
        "schema_version": 2,
        "computed_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "signal": "delta_g_total",
        "dt_ps": dt,
        "n_frames_total": len(series),
        "trajectory_ns": total_ps / 1000.0,
        "series_source": data["source"],
        # The GATE is the structural verdict; energy contributes NO pass/fail.
        "gate": "structural_verdict",
        "structural_verdict": data["verdict"],
        "operative_window": win["operative_window"],
        "operative": op,
        "delta_g_operative_kcal": op["delta_g_mean"],
        "corrected_sem_kcal": op["corrected_sem"],
        "n_eff": op["n_eff"],
        "low_neff": op["low_neff"],
        "begin_end_drift_diagnostic_kcal": win["begin_end_drift_diagnostic_kcal"],
        "preliminary": bool(drifting or short or op["low_neff"]),
        "preliminary_reasons": (["drifting"] if drifting else [])
                               + (["short_trajectory"] if short else [])
                               + (["low_neff"] if op["low_neff"] else []),
        "windows": win["windows"],
        "per_frame_delta_g_total_from_frame0": series,
    }
    (out_dir / "windows.json").write_text(json.dumps(payload, indent=2))
    print(f"  → wrote {out_dir / 'windows.json'}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("md_ids", nargs="*", default=None,
                    help="MD job id(s); defaults to pose0_rep2")
    ap.add_argument("--discard-ps", type=float, default=DEFAULT_FIXED_DISCARD_PS,
                    help="fixed equilibration discard for window (ii), ps")
    ap.add_argument("--recompute", action="store_true",
                    help="re-run MM-GBSA over all frames instead of reusing the "
                         "stored from-frame-0 series")
    args = ap.parse_args(argv)
    ids = args.md_ids or DEFAULT_IDS
    for md_id in ids:
        run_one(md_id, args.discard_ps, reuse=not args.recompute)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
