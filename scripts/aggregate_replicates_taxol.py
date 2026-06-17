"""Per-pose + pooled aggregation of the taxol replicate sweep.

Reads `docs/replicates_taxol_manifest.json`, runs the corrected stability
gate on every replicate, runs single-trajectory MM-GBSA on every gate-
qualifying replicate, then aggregates:

  - PER-POSE: mean ΔG_bind across the pose's 3 replicates + the between-
    replicate σ (the deliverable that replaces the misleading intra-run
    SEM; n=3 is statistically tiny but that IS the experimental design).
  - POOLED across the 9 taxol replicates: a non-parametric bootstrap
    95% CI on the pooled mean ΔG_bind (10 000 resamples; bias-corrected
    percentile CI).

Resumable: skips MD jobs whose FE summary already exists, so partial
sweeps can be re-aggregated.

Outputs:
  docs/replicates_taxol_aggregation.md   — gate + FE + per-pose stats +
                                            pooled bootstrap CI
  docs/replicates_taxol_aggregation.csv  — one row per replicate
                                            (pose, rep, seed, md_id,
                                             gate verdict, ΔG, σ_intra)
"""

from __future__ import annotations

# Pre-load openmm to win the libstdc++ ABI race against the backend's
# psycopg/sqlalchemy chain (see scripts/run_replicates_taxol.py for the
# full rationale). Without this, mmgbsa_runner.compute_md_fe fails with
# obscure ImportErrors inside SystemGenerator.
import openmm as _preload_openmm  # noqa: F401  pylint: disable=unused-import

import csv
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR

DOCS = Path(__file__).resolve().parents[1] / "docs"
DOCS.mkdir(parents=True, exist_ok=True)
MANIFEST_PATH = DOCS / "replicates_taxol_manifest.json"
AGG_MD = DOCS / "replicates_taxol_aggregation.md"
AGG_CSV = DOCS / "replicates_taxol_aggregation.csv"

BOOTSTRAP_N = 10_000
BOOTSTRAP_SEED = 20260616  # deterministic so the CI re-prints identically


def load_manifest() -> Dict[str, Dict[str, Any]]:
    if not MANIFEST_PATH.is_file():
        raise FileNotFoundError(
            f"{MANIFEST_PATH} not found — run scripts/run_replicates_taxol.py first"
        )
    return json.loads(MANIFEST_PATH.read_text())


def fe_summary_path(md_id: str) -> Path:
    return JOBS_DIR / md_id / "free_energy" / "summary.json"


def run_fe_for(md_id: str) -> Optional[Dict[str, Any]]:
    """Run MM-GBSA via compute_md_fe if the FE summary doesn't yet exist.

    Returns the FE result block (the 'result' subdict) or None on failure.
    """
    sp = fe_summary_path(md_id)
    if sp.is_file():
        try:
            f = json.loads(sp.read_text())
            return f
        except Exception:
            pass
    from backend.app.free_energy.mmgbsa_runner import compute_md_fe
    print(f"  computing FE for {md_id[:8]} ...", end="", flush=True)
    t0 = time.perf_counter()
    try:
        block = compute_md_fe(md_id)
    except Exception as e:
        print(f" ✗ {type(e).__name__}: {e}")
        return None
    wall = time.perf_counter() - t0
    sp.parent.mkdir(parents=True, exist_ok=True)
    sp.write_text(json.dumps(block, indent=2))
    print(f" ✓ ΔG={block['result']['delta_g_mean_kcal_per_mol']:+.2f} "
          f"(wall {wall:.0f}s)")
    return block


def run_corrected_gate() -> Dict[str, Dict[str, Any]]:
    """Run scripts/corrected_fe_gate.py's analyse_one on every replicate.

    Returns {md_id: {qualifies, eq_begin_end_delta_a, c1..c5, ...}}.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from corrected_fe_gate import analyze_one, derive_pocket_residues  # type: ignore

    pocket = derive_pocket_residues()
    out: Dict[str, Dict[str, Any]] = {}
    manifest = load_manifest()
    for key, e in manifest.items():
        md_id = e.get("md_job_id")
        if not md_id:
            continue
        r = analyze_one(md_id, pocket)
        out[md_id] = r
    return out


def bootstrap_ci(values: List[float], n_resamples: int = BOOTSTRAP_N,
                 alpha: float = 0.05, seed: int = BOOTSTRAP_SEED) -> Dict[str, float]:
    """Non-parametric percentile bootstrap CI on the mean."""
    rng = random.Random(seed)
    means = []
    n = len(values)
    for _ in range(n_resamples):
        sample = [values[rng.randrange(n)] for _ in range(n)]
        means.append(sum(sample) / n)
    means.sort()
    lo_idx = int(n_resamples * (alpha / 2))
    hi_idx = int(n_resamples * (1 - alpha / 2)) - 1
    return {
        "mean": sum(values) / n,
        "ci_lo": means[lo_idx],
        "ci_hi": means[hi_idx],
        "n_replicates": n,
        "n_resamples": n_resamples,
    }


def main():
    manifest = load_manifest()
    print(f"=== Replicate sweep aggregation ===")
    print(f"  manifest entries : {len(manifest)}")

    # 1) Corrected stability gate per replicate
    gate_results = run_corrected_gate()

    # 2) MM-GBSA per gate-qualifying replicate
    fe_results: Dict[str, Dict[str, Any]] = {}
    for key, e in manifest.items():
        md_id = e.get("md_job_id")
        if not md_id:
            continue
        gate = gate_results.get(md_id, {})
        if "skipped" in gate:
            print(f"  {key}: skipped ({gate['skipped']})")
            continue
        if not gate.get("qualifies"):
            print(f"  {key} ({md_id[:8]}): DOES NOT QUALIFY — running FE anyway "
                  f"(replicate-stat exercise still needs the ΔG number)")
        block = run_fe_for(md_id)
        if block:
            fe_results[md_id] = block

    # 3) Per-replicate rows (CSV)
    rows: List[Dict[str, Any]] = []
    for key, e in sorted(manifest.items()):
        md_id = e.get("md_job_id")
        if not md_id:
            continue
        gate = gate_results.get(md_id, {})
        fe = fe_results.get(md_id, {})
        fe_r = (fe.get("result") or {}) if fe else {}
        rows.append({
            "key": key,
            "pose": e.get("pose"),
            "rep": e.get("rep_idx"),
            "seed": e.get("seed"),
            "md_job_id": md_id,
            "n_frames": e.get("n_frames"),
            "verdict_md": e.get("verdict"),
            "qualifies_corrected_gate": gate.get("qualifies"),
            "eq_begin_end_delta_a": gate.get("eq_begin_end_delta_a"),
            "top5_persistence_eq": gate.get("top5_persistence_eq"),
            "in_pocket_frac_eq": gate.get("in_pocket_frac_eq"),
            "com_final_a": gate.get("com_final_a"),
            "n_frames_fe": fe_r.get("n_frames_evaluated"),
            "delta_g_mean_kcal_per_mol": fe_r.get("delta_g_mean_kcal_per_mol"),
            "delta_g_sem_kcal_per_mol": fe_r.get("delta_g_sem_kcal_per_mol"),
            "delta_g_stddev_kcal_per_mol": fe_r.get("delta_g_stddev_kcal_per_mol"),
            "nonbonded_mean": (fe_r.get("components_mean_kcal_per_mol") or {}).get("nonbonded"),
            "solvation_mean": (fe_r.get("components_mean_kcal_per_mol") or {}).get("solvation"),
        })
    with AGG_CSV.open("w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in rows:
                w.writerow(r)

    # 4) Per-pose aggregation
    per_pose: Dict[int, Dict[str, Any]] = {}
    pooled: List[float] = []
    for pose in (0, 1, 2):
        dgs = [r["delta_g_mean_kcal_per_mol"] for r in rows
               if r["pose"] == pose and r["delta_g_mean_kcal_per_mol"] is not None]
        if not dgs:
            per_pose[pose] = {"n": 0}
            continue
        per_pose[pose] = {
            "n": len(dgs),
            "values": dgs,
            "mean": statistics.fmean(dgs),
            "stdev_between_replicates": (statistics.stdev(dgs) if len(dgs) >= 2 else float("nan")),
        }
        pooled.extend(dgs)
    pooled_stats = bootstrap_ci(pooled) if pooled else None

    # 5) Markdown writeup
    lines = [
        "# Replicate sweep aggregation — 3 reps × 3 taxol poses",
        "",
        "Per-replicate MM-GBSA ΔG_bind (kcal/mol) aggregated into between-",
        "replicate σ + pooled bootstrap 95% CI. Replaces the misleading intra-",
        "run SEM from the single-replicate sweep (`docs/free_energy_qualifying.md`).",
        "",
        "## Per-replicate table",
        "",
        "| Pose | Rep | seed | md_id | gate | Δ(begin−end) | top-5 % | ΔG_bind | σ_intra | n_FE |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        dg = r["delta_g_mean_kcal_per_mol"]
        dg_s = f"{dg:+.2f}" if dg is not None else "—"
        sig = r["delta_g_stddev_kcal_per_mol"]
        sig_s = f"{sig:.2f}" if sig is not None else "—"
        delta = r["eq_begin_end_delta_a"]
        delta_s = f"{delta:.2f}" if delta is not None else "—"
        top5 = r["top5_persistence_eq"]
        top5_s = f"{top5*100:.1f}%" if top5 is not None else "—"
        qual = ("✅" if r["qualifies_corrected_gate"]
                else ("❌" if r["qualifies_corrected_gate"] is False else "—"))
        lines.append(
            f"| {r['pose']} | {r['rep']} | {r['seed']} | `{r['md_job_id'][:8]}` "
            f"| {qual} | {delta_s} Å | {top5_s} | **{dg_s}** | {sig_s} "
            f"| {r['n_frames_fe'] or '—'} |"
        )

    lines += ["", "## Per-pose aggregation", "",
              "| Pose | n_replicates | mean ΔG_bind | σ_between_replicates |",
              "|---|---|---|---|"]
    for pose in (0, 1, 2):
        p = per_pose[pose]
        if p.get("n", 0) == 0:
            lines.append(f"| {pose} | 0 | — | — |")
            continue
        sigma = p["stdev_between_replicates"]
        sigma_s = (f"{sigma:.2f}" if isinstance(sigma, float) and not math.isnan(sigma)
                   else "—")
        lines.append(
            f"| {pose} | {p['n']} | **{p['mean']:+.2f} kcal/mol** | {sigma_s} kcal/mol |"
        )

    lines += ["", "## Pooled bootstrap — 9 taxol replicates", ""]
    if pooled_stats:
        lines.append(
            f"Pooled mean ΔG_bind (n={pooled_stats['n_replicates']} replicates): "
            f"**{pooled_stats['mean']:+.2f} kcal/mol**"
        )
        lines.append(
            f"95% bootstrap CI (n_resamples={pooled_stats['n_resamples']}, "
            f"seed={BOOTSTRAP_SEED}): "
            f"[**{pooled_stats['ci_lo']:+.2f}**, **{pooled_stats['ci_hi']:+.2f}**] kcal/mol"
        )
    else:
        lines.append("_No qualifying ΔG values to bootstrap — sweep incomplete?_")

    lines += [
        "",
        "## Interpretation",
        "",
        "- **σ_between_replicates** is the headline replication uncertainty —"
        " replaces the misleading intra-run SEM (which only captures the per-"
        "frame variance within a single trajectory and dramatically understates"
        " the true uncertainty of single-trajectory MM-GBSA).",
        "- **Pooled bootstrap 95% CI** is the supportable confidence interval"
        " for paclitaxel's ΔG_bind on this pocket under this force-field stack,"
        " over the population of 9 independent 1-ns explicit-TIP3P replicates"
        " across 3 docked poses.",
        "- The single-replicate FE numbers in `docs/free_energy_qualifying.md`"
        " (5bc61f59 / 34840aa1 / 80e53d8a) are NOT in this aggregation —"
        " they're separate reference runs that established the qualifying"
        " baseline. The 9 replicates here are deliberately independent draws"
        " from the same protocol.",
    ]
    AGG_MD.write_text("\n".join(lines))
    print(f"\nWrote {AGG_MD}")
    print(f"Wrote {AGG_CSV}")


if __name__ == "__main__":
    main()
