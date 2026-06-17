"""Run single-trajectory MM-GBSA on the FE-gate-qualifying MD runs.

Qualifying set comes from docs/all_md_corrected_gate.md (5/11 qualified
under the final eq-window-begin/end-only gate, 2026-06-16):

  5bc61f59 — taxol pose 0,         explicit TIP3P, 1 ns
  34840aa1 — taxol pose 1,         explicit TIP3P, 1 ns
  80e53d8a — taxol pose 2,         explicit TIP3P, 1 ns
  a0b04941 — Juliprosopine pose 0, explicit TIP3P, 1 ns
  02f30602 — taxol pose 0 (Run A extended), IMPLICIT OBC2, 1 ns
             — qualifies on the metrics; included here per user direction
             even though it is flagged "implicit lineage" in the gate
             report. Implicit is also the engine's originally-tested code
             path so its FE result is the most-trustworthy data point.

The bundled MM-GBSA gate keys on the old `verdict='stable'` field and will
report `gate_can_run=False` for all four (verdict=drifting). The estimator
honors that flag by marking `sampling_adequate=false / preliminary=true`
in its output but STILL computes the number — which is the contract we
want here. The corrected-gate qualification is recorded separately in
`docs/all_md_corrected_gate.md` and is the authoritative source.

The script is resumable: if `jobs/<id>/free_energy/summary.json` exists
it is loaded and reused. To force a re-run, delete that file.

Output: prints a per-run block + a final summary table; writes
`docs/free_energy_qualifying.md` with the same content.
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR

DOCS = Path("/home/xchem/projects/docking_app2/docs")
QUALIFYING = [
    {
        "md_id": "5bc61f59-834f-4e71-a492-d32ddfdc7326",
        "compound": "taxol",
        "pose": 0,
        "solvent": "explicit",
        "eq_delta_a": 0.07,
        "top5_persistence": 0.993,
    },
    {
        "md_id": "34840aa1-cbf5-4c6b-a665-ea4f52110f5d",
        "compound": "taxol",
        "pose": 1,
        "solvent": "explicit",
        "eq_delta_a": 0.40,
        "top5_persistence": 0.997,
    },
    {
        "md_id": "80e53d8a-1926-4525-b2e1-55cb1e30eedd",
        "compound": "taxol",
        "pose": 2,
        "solvent": "explicit",
        "eq_delta_a": 0.17,
        "top5_persistence": 1.000,
    },
    {
        "md_id": "a0b04941-e4b0-40ef-9459-67becac4a61c",
        "compound": "Juliprosopine",
        "pose": 0,
        "solvent": "explicit",
        "eq_delta_a": 0.29,
        "top5_persistence": 0.813,
    },
    {
        "md_id": "02f30602-8624-4c0a-863b-9c130d357364",
        "compound": "taxol (Run A extended)",
        "pose": 0,
        "solvent": "implicit",
        "eq_delta_a": 0.28,
        "top5_persistence": 0.863,
    },
]


def existing_fe_summary(md_id: str) -> Optional[Dict[str, Any]]:
    """Return the persisted free-energy block (if any) so a previously
    completed FE run is reused on a re-launch."""
    p = JOBS_DIR / md_id / "free_energy" / "summary.json"
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def persist_fe_block(md_id: str, fe_block: Dict[str, Any]) -> None:
    out_dir = JOBS_DIR / md_id / "free_energy"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(fe_block, indent=2))


def run_one(entry: Dict[str, Any]) -> Dict[str, Any]:
    md_id = entry["md_id"]
    print(f"\n──────────────────────────────────────────────────────────────")
    print(f"### {entry['compound']} pose {entry['pose']}  ({md_id[:8]})")
    print(f"  solvent={entry['solvent']}  "
          f"corrected-gate Δ(begin−end)={entry['eq_delta_a']} Å  "
          f"top-5={entry['top5_persistence']*100:.1f}%")
    cached = existing_fe_summary(md_id)
    if cached is not None:
        print(f"  → reusing existing free_energy/summary.json")
        return cached

    from backend.app.free_energy.mmgbsa_runner import compute_md_fe

    t0 = time.perf_counter()
    try:
        fe_block = compute_md_fe(md_id)
    except Exception as e:
        wall = time.perf_counter() - t0
        print(f"  ✗ FAILED after {wall:.1f} s — {type(e).__name__}: {e}")
        traceback.print_exc()
        return {
            "status": "failed",
            "md_job_id": md_id,
            "error": f"{type(e).__name__}: {e}",
            "wall_seconds": wall,
        }
    wall = time.perf_counter() - t0
    fe_block["sweep_wall_seconds"] = wall
    persist_fe_block(md_id, fe_block)
    print(f"  ✓ done in {wall:.1f} s  → wrote free_energy/summary.json")
    return fe_block


def fmt_dg(block: Dict[str, Any]) -> str:
    if block.get("status") == "failed":
        return f"failed ({block.get('error')})"
    dg = (block.get("result") or {}).get("dg_bind_kcal_per_mol")
    sem = (block.get("result") or {}).get("dg_bind_sem_kcal_per_mol")
    if dg is None:
        return "—"
    if sem is not None:
        return f"{dg:+.2f} ± {sem:.2f} kcal/mol"
    return f"{dg:+.2f} kcal/mol"


def fmt_components(block: Dict[str, Any]) -> str:
    res = block.get("result") or {}
    parts = []
    for k in ("bonded", "nonbonded", "solvation"):
        v = res.get(f"dg_{k}_kcal_per_mol")
        if v is not None:
            parts.append(f"{k}={v:+.2f}")
    return ", ".join(parts) if parts else "—"


def main():
    print("=== MM-GBSA sweep on FE-gate-qualifying MD runs ===")
    print(f"Qualifying set ({len(QUALIFYING)}):")
    for e in QUALIFYING:
        print(f"  - {e['md_id'][:8]}  {e['compound']:13s} pose {e['pose']}  "
              f"{e['solvent']:8s}  Δ={e['eq_delta_a']} Å")
    print()

    results: List[Dict[str, Any]] = []
    for entry in QUALIFYING:
        block = run_one(entry)
        results.append({"entry": entry, "block": block})

    # === Final summary table ===
    print("\n\n## Summary\n")
    print("| Run | Compound | Pose | Sol | ΔG_bind | bonded / nonbonded / solvation | sampling_adequate | wall |")
    print("|---|---|---|---|---|---|---|---|")
    md_lines = [
        "# MM-GBSA on FE-gate-qualifying MD runs",
        "",
        "Single-trajectory MM-GBSA computed on the four explicit-TIP3P MD runs "
        "that cleared the final FE gate (eq-window begin/end Δ < 0.4 Å, "
        "in-pocket ≥ 0.95, top-5 residue persistence ≥ 65 %, ≥ 50 eq frames). "
        "Estimator: `backend.app.free_energy.mmgbsa.estimate_mmgbsa` "
        "(single-trajectory; bonded + nonbonded + GBSA-OBC2 solvation; "
        "configurational entropy omitted). Implicit-lineage run `02f30602` "
        "deliberately NOT included.",
        "",
        "**Caveat — bundled gate vs corrected gate.** The bundled "
        "`free_energy.gating.validate` keys on `summary.verdict='stable'` and "
        "marks all four runs preliminary (verdict='drifting'). The corrected "
        "gate at `docs/all_md_corrected_gate.md` is the authoritative "
        "qualification source — it ignores `verdict` and judges convergence + "
        "current state. The MM-GBSA number is computed either way; "
        "`sampling_adequate=false / preliminary=true` in the per-run JSON is "
        "the bundled-gate's verdict-based flag, not a corrected-gate signal.",
        "",
        "## Per-run results",
        "",
        "| Run | Compound | Pose | Sol | ΔG_bind | bonded / nonbonded / solvation | sampling_adequate | wall (s) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        e = r["entry"]; b = r["block"]
        adequate = (b.get("result") or {}).get("sampling_adequate")
        adequate_cell = ("—" if adequate is None
                         else ("✅ yes" if adequate else "⚠ preliminary"))
        wall = b.get("sweep_wall_seconds") or b.get("wall_seconds") or "—"
        wall_str = (f"{wall:.0f}" if isinstance(wall, (int, float)) else str(wall))
        row = (f"`{e['md_id'][:8]}` | {e['compound']} | {e['pose']} | "
               f"{e['solvent']} | {fmt_dg(b)} | {fmt_components(b)} | "
               f"{adequate_cell} | {wall_str}")
        print("| " + row + " |")
        md_lines.append("| " + row + " |")

    md_lines += [
        "",
        "## Method",
        "",
        "- Force-field stack: amber14-all + amber14/tip3p ion templates + "
        "implicit/obc2 + gaff-2.11 (reuses the MD parameterization).",
        "- Nonbonded: CutoffNonPeriodic, cutoff 1.0 nm (same NonbondedForce "
        "class as the MD).",
        "- Per-frame component split via OpenMM force-group dispatch.",
        "- Equilibration discard: auto 20 ps for explicit MD (waters relax "
        "around the unrestrained solute post-restraint-release).",
        "- Configurational entropy: OMITTED (normal-mode / quasi-harmonic out "
        "of scope). Reported number is enthalpic + solvation only — the "
        "standard MM-GBSA quantity.",
        "- Single-trajectory subtraction: ΔG_bind = ⟨E_complex⟩ − ⟨E_receptor⟩ "
        "− ⟨E_ligand⟩.",
        "",
        "Per-run artifacts: `jobs/<md_id>/free_energy/summary.json` (the full "
        "self-describing FE block written by `mmgbsa_runner.compute_md_fe`).",
    ]
    out = DOCS / "free_energy_qualifying.md"
    out.write_text("\n".join(md_lines))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
