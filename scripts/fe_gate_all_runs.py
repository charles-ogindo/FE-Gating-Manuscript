"""Apply the free-energy (FE) gate to every openmm_full MD run in the
catalogue. Writes two outputs under docs/:
  - all_md_fe_gate.csv (machine-readable)
  - all_md_fe_gate.md  (paste-friendly)

Honest reporting: every openmm_full run is included, even smoke-scale
ones and the wrong-pocket d96d719a-lineage rows that have no frames on
disk. The gate's own messages explain WHY each row passes or declines.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR
from backend.app.db.session import SessionLocal
from backend.app.db.models import MDRun
from backend.app.free_energy.gating import validate
from backend.app.free_energy import (
    DEFAULT_MIN_TOTAL_FRAMES, DEFAULT_MIN_WINDOW_FRAMES,
    LIGAND_MAX_RMSD_WARN_A, LIGAND_MAX_RMSD_BLOCK_A,
    BACKBONE_FINAL_RMSD_WARN_A, WINDOW_DURATION_WARN_PS,
)

DOCS = Path("/home/xchem/projects/docking_app2/docs")

# Same compound-naming map as all_md_diagnostics.py
DOCKING_TO_LIGAND_TO_COMMON = {
    "4a37bb0c-7655-411e-a238-1ff5b7bab910": {
        "lig_000000": "taxol",
        "lig_000001": "Juliprosopine",
        "lig_000002": "Primaquine",
    },
    "d96d719a": {"lig_000000": "taxol (wrong-pocket)"},
}


def compound_for(docking_job_id, ligand_name) -> str:
    if docking_job_id is None or ligand_name is None:
        return ligand_name or "?"
    key = str(docking_job_id)
    for dock_id, lig_map in DOCKING_TO_LIGAND_TO_COMMON.items():
        if key.startswith(dock_id) or dock_id.startswith(key):
            return lig_map.get(ligand_name, ligand_name)
    return ligand_name


def main():
    # Sweep all md_runs marked openmm_full (engine column in md_runs).
    # We also include rows where engine is None (could be a failed openmm run
    # that didn't reach the summary writeback) so the table is fully honest.
    with SessionLocal() as db:
        all_runs = db.query(MDRun).order_by(MDRun.id).all()
        targets = []
        for m in all_runs:
            mid = str(m.id)
            # Filter: include if engine=='openmm_full' OR engine is None and the
            # filesystem summary.json says openmm_full. Skip pure surrogate
            # rows (they don't reach the FE pipeline by gate design — but
            # we'll let the gate's own MD_NOT_OPENMM_FULL reason fire if a
            # caller accidentally passes them).
            summary_engine = None
            sp = JOBS_DIR / mid / "md" / "summary.json"
            if sp.is_file():
                try:
                    summary_engine = (json.loads(sp.read_text())
                                      .get("engine") or {}).get("kind")
                except Exception:
                    pass
            if (m.engine or summary_engine) == "openmm_full":
                targets.append(m)
    print(f"openmm_full MD runs to gate: {len(targets)}")

    print(f"\nGate thresholds (from backend.app.free_energy.__init__):")
    print(f"  DEFAULT_MIN_TOTAL_FRAMES        = {DEFAULT_MIN_TOTAL_FRAMES}")
    print(f"  DEFAULT_MIN_WINDOW_FRAMES       = {DEFAULT_MIN_WINDOW_FRAMES}")
    print(f"  LIGAND_MAX_RMSD_WARN_A          = {LIGAND_MAX_RMSD_WARN_A}")
    print(f"  LIGAND_MAX_RMSD_BLOCK_A         = {LIGAND_MAX_RMSD_BLOCK_A}")
    print(f"  BACKBONE_FINAL_RMSD_WARN_A      = {BACKBONE_FINAL_RMSD_WARN_A}")
    print(f"  WINDOW_DURATION_WARN_PS         = {WINDOW_DURATION_WARN_PS}")

    rows: List[Dict[str, Any]] = []
    for m in targets:
        mid = str(m.id)
        compound = compound_for(m.docking_job_id, m.ligand_name)
        try:
            g = validate(mid).to_dict()
        except Exception as e:
            g = {"can_run": False, "reasons": [
                {"key": "GATE_CRASH", "message": str(e)[:200]}],
                "warnings": [], "md_status": None, "md_verdict": None,
                "resolved_sampling_plan": None}
        # Pose-RMSD the gate looked at (from md/summary.json):
        sp = JOBS_DIR / mid / "md" / "summary.json"
        rmsd_pose_max = None
        rmsd_pose_final = None
        rmsd_bb_final = None
        n_frames = None
        production_ps = None
        if sp.is_file():
            try:
                s = json.loads(sp.read_text())
                metrics = s.get("metrics") or {}
                rmsd_pose_max = metrics.get("rmsd_ligand_pose_max_a")
                rmsd_pose_final = metrics.get("rmsd_ligand_pose_final_a")
                rmsd_bb_final = metrics.get("rmsd_backbone_final_a")
                n_frames = s.get("n_frames")
                production_ps = (s.get("settings") or {}).get("production_ps")
            except Exception:
                pass
        plan = g.get("resolved_sampling_plan") or {}
        rows.append({
            "md_job_id": mid,
            "compound": compound,
            "ligand_name": m.ligand_name,
            "pose_rank": m.pose_rank,
            "docking_job_id": str(m.docking_job_id) if m.docking_job_id else None,
            "md_engine": m.engine or "—",
            "md_status": g.get("md_status") or "—",
            "md_verdict": g.get("md_verdict") or "—",
            "production_ps": production_ps,
            "n_frames": n_frames,
            "rmsd_ligand_pose_max_a": rmsd_pose_max,
            "rmsd_ligand_pose_final_a": rmsd_pose_final,
            "rmsd_backbone_final_a": rmsd_bb_final,
            "can_run": g.get("can_run", False),
            "n_hard_gate_reasons": len(g.get("reasons") or []),
            "n_soft_gate_warnings": len(g.get("warnings") or []),
            "hard_gate_keys": [r.get("key") for r in (g.get("reasons") or [])],
            "soft_gate_keys": [w.get("key") for w in (g.get("warnings") or [])],
            "hard_gate_reasons_text": " | ".join(
                f"{r.get('key')}: {r.get('message')[:120]}"
                for r in (g.get("reasons") or [])
            ),
            "soft_gate_warnings_text": " | ".join(
                f"{w.get('key')}: {w.get('message')[:120]}"
                for w in (g.get("warnings") or [])
            ),
            "n_frames_sampled": (plan or {}).get("n_frames_sampled"),
            "stride": (plan or {}).get("stride"),
            "window_last_fraction": (plan or {}).get("window_last_fraction"),
        })

    # Sort: passing runs first (so the few that can_run get top billing), then
    # by compound to group taxol / Juli / Prim / wrong-pocket.
    def sort_key(r):
        passing = 0 if r["can_run"] else 1
        comp = r["compound"] or ""
        if "wrong-pocket" in comp:
            comp_order = 9
        elif comp == "taxol":
            comp_order = 0
        elif comp == "Juliprosopine":
            comp_order = 1
        elif comp == "Primaquine":
            comp_order = 2
        else:
            comp_order = 5
        pose = r.get("pose_rank") if r.get("pose_rank") is not None else 99
        return (passing, comp_order, pose, r["md_job_id"])
    rows.sort(key=sort_key)

    csv_path = DOCS / "all_md_fe_gate.csv"
    cols = [
        "md_job_id", "compound", "ligand_name", "pose_rank", "docking_job_id",
        "md_engine", "md_status", "md_verdict",
        "production_ps", "n_frames", "n_frames_sampled", "stride",
        "window_last_fraction",
        "rmsd_ligand_pose_max_a", "rmsd_ligand_pose_final_a",
        "rmsd_backbone_final_a",
        "can_run", "n_hard_gate_reasons", "n_soft_gate_warnings",
        "hard_gate_keys", "soft_gate_keys",
        "hard_gate_reasons_text", "soft_gate_warnings_text",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            row = []
            for c in cols:
                v = r.get(c)
                if isinstance(v, list):
                    v = ";".join(v) if v else ""
                row.append(v if v is not None else "")
            w.writerow(row)
    print(f"\nWrote {csv_path}")

    # Markdown
    def fmt(v):
        if v is None or v == "":
            return "—"
        if isinstance(v, bool):
            return "**pass**" if v else "decline"
        if isinstance(v, float):
            if abs(v) < 1e-4 and v != 0:
                return f"{v:.1e}"
            return f"{v:.2f}"
        if isinstance(v, list):
            return ", ".join(v) if v else "—"
        return str(v)

    md_path = DOCS / "all_md_fe_gate.md"
    lines = [
        "# FE gate applied to every openmm_full MD run",
        "",
        "Loops `backend.app.free_energy.gating.validate(md_job_id)` over every "
        f"MD run flagged engine=`openmm_full` in the catalogue ({len(rows)} "
        "runs). Sourced from db.MDRun ∪ filesystem summary.json (engine).",
        "",
        "## Gate thresholds (verbatim from `backend.app.free_energy.__init__`)",
        "",
        f"- `DEFAULT_MIN_TOTAL_FRAMES`   = **{DEFAULT_MIN_TOTAL_FRAMES}**  "
        "  (hard gate — ≥ 200 frames before any windowing)",
        f"- `DEFAULT_MIN_WINDOW_FRAMES`  = **{DEFAULT_MIN_WINDOW_FRAMES}**   "
        "  (hard gate — ≥ 50 frames in the sampled window)",
        f"- `LIGAND_MAX_RMSD_WARN_A`     = **{LIGAND_MAX_RMSD_WARN_A}** Å  "
        " (soft gate — warns above this)",
        f"- `LIGAND_MAX_RMSD_BLOCK_A`    = **{LIGAND_MAX_RMSD_BLOCK_A}** Å  "
        " (soft gate — *blocks* if ligand max RMSD ≥ this — \"unbound\")",
        f"- `BACKBONE_FINAL_RMSD_WARN_A` = **{BACKBONE_FINAL_RMSD_WARN_A}** Å"
        "  (soft gate — backbone-not-equilibrated warning)",
        f"- `WINDOW_DURATION_WARN_PS`    = **{WINDOW_DURATION_WARN_PS}** ps "
        "  (soft gate — short post-eq window warning)",
        "",
        "The gate's hard gates (must pass all):",
        "  A `MD_SUMMARY_MISSING` — md/summary.json present + parseable",
        "  B `MD_NOT_COMPLETE`    — summary.status == 'completed'",
        "  C `MD_NOT_STABLE`      — summary.verdict == 'stable'",
        "  D `MD_NOT_OPENMM_FULL` — summary.engine.kind == 'openmm_full'",
        "  E `MD_FRAMES_MISSING` / `MD_FRAMES_EMPTY` / `MD_FRAME_UNREADABLE`",
        "  F `INSUFFICIENT_FRAMES` — n_total ≥ 200 AND n_sampled ≥ 50",
        "",
        "Soft gates emit warnings but do NOT block: `LIGAND_DRIFT` "
        "(2-5 Å warn), `BACKBONE_NOT_EQUILIBRATED` (BB Cα RMSD final ≥ "
        "2.5 Å), `OUTLIER_FRAMES`, `SHORT_WINDOW_TIME` (< 1 ns sampled), "
        "`STRIDE_MANUAL`, `EXPERIMENTAL_ENGINE` (standing flag — every "
        "passing run gets it; this is the \"experimental / approximate\" "
        "chip the dashboard shows on the FE estimate).",
        "",
        "## Result",
        "",
        f"Out of {len(rows)} openmm_full runs, **"
        f"{sum(1 for r in rows if r['can_run'])} pass** the FE gate.",
        "",
        "| job_id (8) | compound | pose | dur(ps) | n_frm | ligand pose RMSD (fin/max) | BB final | gate | hard-gate reasons | soft-gate warnings |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        rmsd_pair = f"{fmt(r['rmsd_ligand_pose_final_a'])} / {fmt(r['rmsd_ligand_pose_max_a'])}"
        lines.append("| " + " | ".join([
            f"`{r['md_job_id'][:8]}`",
            str(r["compound"] or "—"),
            str(r.get("pose_rank") if r.get("pose_rank") is not None else "—"),
            fmt(r["production_ps"]),
            fmt(r["n_frames"]),
            rmsd_pair,
            fmt(r["rmsd_backbone_final_a"]),
            ("✅ **pass**" if r["can_run"] else "❌ decline"),
            fmt(r.get("hard_gate_keys") or []),
            fmt(r.get("soft_gate_keys") or []),
        ]) + " |")

    lines += [
        "",
        "## Per-run gate-message detail (the why)",
        "",
    ]
    for r in rows:
        h = "; ".join(
            f"`{k}`" for k in (r.get("hard_gate_keys") or [])) or "—"
        s = "; ".join(
            f"`{k}`" for k in (r.get("soft_gate_keys") or [])) or "—"
        lines.append(
            f"- **`{r['md_job_id'][:8]}`** — {r['compound']} pose "
            f"{r.get('pose_rank') if r.get('pose_rank') is not None else '—'} "
            f"({fmt(r['production_ps'])} ps, {fmt(r['n_frames'])} frames, "
            f"verdict=`{r['md_verdict']}`): "
            + ("✅ **PASS**" if r["can_run"] else "❌ **DECLINE**")
            + f"  hard=[{h}], soft=[{s}]"
        )
        if r.get("hard_gate_reasons_text"):
            lines.append(f"  - hard-gate msgs: {r['hard_gate_reasons_text']}")
        if r.get("soft_gate_warnings_text"):
            lines.append(f"  - soft-gate msgs: {r['soft_gate_warnings_text']}")
    md_path.write_text("\n".join(lines))
    print(f"Wrote {md_path}")

    # Brief summary on stdout
    print(f"\nQuick stats:")
    by_gate_key = {}
    for r in rows:
        for k in (r.get("hard_gate_keys") or []):
            by_gate_key[k] = by_gate_key.get(k, 0) + 1
    for k, c in sorted(by_gate_key.items(), key=lambda kv: -kv[1]):
        print(f"  hard-gate {k}: {c} runs blocked by this key")
    if any(r["can_run"] for r in rows):
        print("\nPassing runs:")
        for r in rows:
            if r["can_run"]:
                print(f"  ✅ {r['md_job_id'][:8]}  {r['compound']} pose "
                      f"{r.get('pose_rank')}  "
                      f"({fmt(r['production_ps'])} ps, {fmt(r['n_frames'])} fr)")
    else:
        print("\nNo run passes the gate currently.")


if __name__ == "__main__":
    main()
