"""All-MD-runs diagnostics consolidator.

Sweeps the entire MD catalogue (DB md_runs ∪ filesystem jobs/<id>/md/) and
extracts per-run metrics into:
  - docs/all_md_diagnostics.csv  (machine-readable, committed alongside the
                                  B9 stability diagnostics)
  - docs/all_md_diagnostics.md   (human-readable, grouped)

Honest reporting policy:
  - Include EVERY md_run ever queued (DB ∪ FS), including failed, empty,
    surrogate, smoke-test, and wrong-pocket runs. Each gets a Notes
    column flagging its caveat.
  - For runs with ≥ 11 frames AND engine=openmm_full, run the
    per-run pocket-aligned analysis (same pocket-Cα Kabsch convention as
    scripts/stability_diagnostics.py — pocket = receptor residues with
    any heavy atom within 5 Å of any ligand heavy atom in the run's own
    frame 0, intersected with the pocket-residue set derived from the
    canonical reference run, 5bc61f59 explicit taxol pose 0).
  - For shorter / surrogate / DB-only / failed runs, fill in what the
    summary.json carries and mark the deeper columns "—".
  - Convergence test: 1st-half vs 2nd-half block-averaged pose-pkt RMSD.
    "yes" if |Δ| ≤ 0.3 Å AND last-half slope |≤ 1e-3 Å/ps|.
    "≈" if borderline. "no" otherwise. Failed/empty runs: "—".
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR
from backend.app.db.session import SessionLocal
from backend.app.db.models import MDRun

# Reuse helpers from stability_diagnostics
from stability_diagnostics import (  # type: ignore
    parse_pdb_frame,
    kabsch_rotation,
    align_via,
    rmsd_after_alignment,
    POCKET_RADIUS_A,
    CONTACT_CUTOFF_A,
    BURIED_CUTOFF_A,
    IN_POCKET_COM_A,
    EQ_WINDOW_LAST_FRAC,
)

DOCS_ROOT = Path(__file__).resolve().parents[1] / "docs"
DOCS_ROOT.mkdir(parents=True, exist_ok=True)

# Canonical reference for the pocket-residue set definition.
REFERENCE_RUN_ID = "5bc61f59-834f-4e71-a492-d32ddfdc7326"  # explicit taxol pose 0

# Compound identification — md_runs.ligand_name + docking_job_id lineage.
DOCKING_TO_LIGAND_TO_COMMON = {
    "4a37bb0c-7655-411e-a238-1ff5b7bab910": {
        "lig_000000": "taxol",
        "lig_000001": "Juliprosopine",
        "lig_000002": "Primaquine",
    },
    # d96d719a is the wrong-pocket lineage — kept ligand_name=lig_000000
    # but the receptor is a different construct; flagged in notes.
    "d96d719a": {"lig_000000": "taxol (wrong-pocket)"},
}


def common_name_for(docking_job_id: Optional[str], ligand_name: Optional[str]
                    ) -> str:
    if docking_job_id is None or ligand_name is None:
        return ligand_name or "?"
    key = str(docking_job_id)
    # Match either full UUID or its prefix
    for dock_id, lig_map in DOCKING_TO_LIGAND_TO_COMMON.items():
        if key.startswith(dock_id) or dock_id.startswith(key):
            if ligand_name in lig_map:
                return lig_map[ligand_name]
    return ligand_name


def derive_pocket_residues_from_reference() -> set:
    """Pocket residue set = receptor residues with any heavy atom within
    POCKET_RADIUS_A of any ligand heavy atom in the canonical reference
    run's frame 0. Used as the consistent definition across runs."""
    ref_dir = JOBS_DIR / REFERENCE_RUN_ID / "md"
    f0 = parse_pdb_frame(ref_dir / "frames" / "frame_000.pdb")
    lig = np.asarray([a["xyz"] for a in f0.ligand_heavy])
    rec_heavy_keys = []
    rec_heavy_pos = []
    for a in f0.receptor_heavy:
        rec_heavy_keys.append((a["chain"], a["resseq"]))
        rec_heavy_pos.append(a["xyz"])
    rec_heavy_pos = np.asarray(rec_heavy_pos)
    dmat = np.sqrt(((lig[:, None, :] - rec_heavy_pos[None, :, :]) ** 2)
                   .sum(axis=2))
    close = dmat.min(axis=0) < POCKET_RADIUS_A
    keys = set()
    for i, is_close in enumerate(close):
        if is_close:
            keys.add(rec_heavy_keys[i])
    return keys


def find_equilibration_frame(pose_pkt: List[float], min_frames: int = 20) -> Optional[int]:
    """Crude apparent-equilibration heuristic: walk a sliding 10-frame window
    over pose-pkt RMSD; return the first frame index where the window's mean
    is within 0.3 Å of the trajectory's final-10-frame mean, AND the window's
    std is ≤ 0.5 Å. Returns None if no such frame found (still drifting)."""
    if len(pose_pkt) < min_frames:
        return None
    arr = np.asarray(pose_pkt)
    final_mean = float(np.nanmean(arr[-10:]))
    for i in range(len(arr) - 10):
        window = arr[i:i + 10]
        if np.isnan(window).any():
            continue
        if (abs(np.mean(window) - final_mean) <= 0.3
                and np.std(window) <= 0.5):
            return i
    return None


def analyze_one_md(md_id: str, shared_pocket: set) -> Dict[str, Any]:
    """Compute the full diagnostic row. Falls back gracefully when frames
    are missing or short."""
    md_dir = JOBS_DIR / md_id / "md"
    row: Dict[str, Any] = {
        "md_job_id": md_id,
        "summary_present": (md_dir / "summary.json").is_file(),
        "frames_dir_present": (md_dir / "frames").is_dir(),
    }

    # Pull whatever summary.json carries first — works for surrogate, openmm,
    # failed runs, etc.
    summary: Dict[str, Any] = {}
    if row["summary_present"]:
        try:
            summary = json.loads((md_dir / "summary.json").read_text())
        except Exception:
            summary = {}
    settings = summary.get("settings") or {}
    metrics = summary.get("metrics") or {}
    engine = (summary.get("engine") or {}).get("kind") or "—"
    verdict = summary.get("verdict") or "—"
    duration_ps = settings.get("production_ps")
    snap_ps = settings.get("snapshot_every_ps")
    solvent = settings.get("solvent", "implicit")
    n_frames_summary = summary.get("n_frames")

    row.update({
        "engine": engine,
        "verdict": verdict,
        "duration_ps": duration_ps,
        "snapshot_every_ps": snap_ps,
        "solvent": solvent,
        "n_frames_summary": n_frames_summary,
        # Summary-side metrics for runs we can't deep-analyze.
        "summary_rmsd_pose_final_a": metrics.get("rmsd_ligand_pose_final_a"),
        "summary_rmsd_pose_max_a":   metrics.get("rmsd_ligand_pose_max_a"),
        "summary_rmsd_internal_final_a": metrics.get("rmsd_ligand_internal_final_a"),
        "summary_rmsd_backbone_final_a": metrics.get("rmsd_backbone_final_a"),
        "summary_hbond_persistence_frac": metrics.get("hbond_persistence_frac"),
        "wall_seconds": summary.get("wall_seconds"),
    })

    frame_paths = sorted((md_dir / "frames").glob("frame_*.pdb")) \
        if row["frames_dir_present"] else []
    row["n_frames_on_disk"] = len(frame_paths)

    # Eligibility for deep analysis: real frames + reasonable count.
    # The B9-class analysis assumes pocket Cα alignment is meaningful —
    # surrogate runs only perturb the ligand (rigid receptor), so the
    # pocket-aligned columns degrade to "ligand wiggle in fixed frame"
    # rather than binding-mode drift. Mark separately.
    can_deep = (engine == "openmm_full" and len(frame_paths) >= 11)

    if not can_deep:
        row.update({
            "pose_pkt_final_a": None, "pose_pkt_max_a": None,
            "com_displ_final_a": None, "com_displ_max_a": None,
            "converged": "—", "first_half_mean_a": None,
            "second_half_mean_a": None, "slope_a_per_ps": None,
            "in_pocket_frac_final": None, "buried_frac_final": None,
            "contacts_n_initial": None, "contacts_retained_frac_mean_eq": None,
            "equilibration_frame_idx": None,
        })
        return row

    # ---- deep analysis ----
    snapshots = [parse_pdb_frame(p) for p in frame_paths]
    n = len(snapshots)
    f0 = snapshots[0]
    rec0_ca = f0.receptor_ca
    lig0 = f0.ligand_heavy
    rec0_ca_pos = np.asarray([a["xyz"] for a in rec0_ca])
    lig0_pos = np.asarray([a["xyz"] for a in lig0])
    rec0_heavy = np.asarray([a["xyz"] for a in f0.receptor_heavy])
    rec0_heavy_keys = [(a["chain"], a["resseq"], a["resname"], a["name"])
                       for a in f0.receptor_heavy]

    rec_keys_in_run = {(a["chain"], a["resseq"]) for a in rec0_ca}
    pocket_keys = shared_pocket & rec_keys_in_run
    pocket_order = [(a["chain"], a["resseq"]) for a in rec0_ca
                    if (a["chain"], a["resseq"]) in pocket_keys]
    pocket0_pos = np.asarray([
        next(a["xyz"] for a in rec0_ca if (a["chain"], a["resseq"]) == k)
        for k in pocket_order
    ])

    # Initial contacts (ligand × protein ≤ 4 Å)
    dmat0 = np.sqrt(((lig0_pos[:, None, :] - rec0_heavy[None, :, :]) ** 2)
                    .sum(axis=2))
    contact_pairs = np.argwhere(dmat0 < CONTACT_CUTOFF_A)
    init_contacts = [(int(li), rec0_heavy_keys[int(ri)])
                     for li, ri in contact_pairs]
    contact_persist = [0] * len(init_contacts)

    pose_pkt: List[float] = []
    com_displ: List[float] = []
    in_pocket: List[bool] = []
    buried: List[float] = []
    contacts_retained: List[float] = []

    times_ps = [i * (snap_ps if snap_ps else 1.0) for i in range(n)]

    for idx, snap in enumerate(snapshots):
        rec_ca_pos = np.asarray([a["xyz"] for a in snap.receptor_ca])
        lig = np.asarray([a["xyz"] for a in snap.ligand_heavy])
        if (rec_ca_pos.shape != rec0_ca_pos.shape
                or lig.shape != lig0_pos.shape
                or len(pocket_order) == 0):
            pose_pkt.append(np.nan); com_displ.append(np.nan)
            in_pocket.append(False); buried.append(np.nan)
            contacts_retained.append(np.nan)
            continue
        # Pocket-Cα Kabsch fit; apply to ligand for pose RMSD + COM
        pocket_pos_i = np.asarray([
            next(a["xyz"] for a in snap.receptor_ca
                 if (a["chain"], a["resseq"]) == k)
            for k in pocket_order
        ])
        p_cent = pocket_pos_i.mean(axis=0); q_cent = pocket0_pos.mean(axis=0)
        P = pocket_pos_i - p_cent; Q = pocket0_pos - q_cent
        R = kabsch_rotation(P, Q)
        def xfm(arr):
            return (arr - p_cent) @ R.T + q_cent
        lig_xfm = xfm(lig)
        pose_pkt.append(float(np.sqrt(((lig_xfm - lig0_pos) ** 2)
                                      .sum(axis=1).mean())))
        com_i = lig_xfm.mean(axis=0); com_0 = lig0_pos.mean(axis=0)
        d = float(np.linalg.norm(com_i - com_0))
        com_displ.append(d)
        in_pocket.append(d < IN_POCKET_COM_A)
        rec_heavy = np.asarray([a["xyz"] for a in snap.receptor_heavy])
        # Buried = fraction of ligand heavy atoms with any protein heavy
        # atom within 4.5 Å.
        dmat_i = np.sqrt(((lig[:, None, :] - rec_heavy[None, :, :]) ** 2)
                         .sum(axis=2))
        buried.append(float((dmat_i.min(axis=1) < BURIED_CUTOFF_A).mean()))
        # Contacts retained from frame 0 set
        rec_lookup = {(a["chain"], a["resseq"], a["resname"], a["name"]): a["xyz"]
                      for a in snap.receptor_heavy}
        retained = 0
        for ci, (li_idx, key) in enumerate(init_contacts):
            ri_xyz = rec_lookup.get(key)
            if ri_xyz is None: continue
            l_xyz = lig[li_idx]
            dd = ((ri_xyz[0] - l_xyz[0]) ** 2 + (ri_xyz[1] - l_xyz[1]) ** 2
                  + (ri_xyz[2] - l_xyz[2]) ** 2) ** 0.5
            if dd < CONTACT_CUTOFF_A:
                retained += 1
                contact_persist[ci] += 1
        contacts_retained.append(
            retained / len(init_contacts) if init_contacts else np.nan)

    arr = np.asarray(pose_pkt)
    half = n // 2
    fh = float(np.nanmean(arr[:half]))
    sh = float(np.nanmean(arr[half:]))
    if half < n - 1:
        x = np.asarray(times_ps[half:])
        y = arr[half:]
        m = ~np.isnan(y)
        slope, _ = (np.polyfit(x[m], y[m], 1) if m.sum() >= 2
                    else (np.nan, np.nan))
    else:
        slope = np.nan

    # Convergence verdict
    if not np.isfinite(fh + sh + slope):
        conv = "—"
    elif abs(sh - fh) <= 0.3 and abs(slope) <= 1e-3:
        conv = "yes"
    elif abs(sh - fh) <= 0.7 and abs(slope) <= 3e-3:
        conv = "≈"
    else:
        conv = "no"

    eq_idx = find_equilibration_frame(pose_pkt)
    eq_window = int(round(n * (1.0 - EQ_WINDOW_LAST_FRAC)))
    in_pocket_eq_final = bool(in_pocket[-1])
    buried_final = float(buried[-1]) if buried else np.nan
    contacts_eq = float(np.nanmean(contacts_retained[eq_window:])) \
        if eq_window < n else float("nan")

    row.update({
        "pose_pkt_final_a": float(arr[-1]) if np.isfinite(arr[-1]) else None,
        "pose_pkt_max_a": float(np.nanmax(arr)) if np.isfinite(np.nanmax(arr)) else None,
        "com_displ_final_a": float(com_displ[-1]) if com_displ else None,
        "com_displ_max_a": float(np.nanmax(com_displ)) if com_displ else None,
        "converged": conv,
        "first_half_mean_a": fh, "second_half_mean_a": sh,
        "slope_a_per_ps": float(slope) if np.isfinite(slope) else None,
        "in_pocket_final": in_pocket_eq_final,
        "buried_frac_final": buried_final,
        "contacts_n_initial": len(init_contacts),
        "contacts_retained_frac_mean_eq": contacts_eq,
        "equilibration_frame_idx": eq_idx,
        "pocket_residues_intersected": len(pocket_order),
    })
    return row


def main():
    print(f"Reference run for pocket: {REFERENCE_RUN_ID}")
    shared_pocket = derive_pocket_residues_from_reference()
    print(f"Pocket residue set (from reference frame 0, 5 Å of any ligand "
          f"heavy atom): {len(shared_pocket)}")

    # Inventory: alive FS + DB md_runs (capture DB-only rows too).
    alive_ids = {p.name for p in JOBS_DIR.iterdir()
                 if p.is_dir() and p.name != "archive"
                 and (p / "md" / "summary.json").is_file()}
    with SessionLocal() as db:
        db_runs = db.query(MDRun).all()
        db_info = {str(m.id): {
            "docking_job_id": str(m.docking_job_id) if m.docking_job_id else None,
            "ligand_name": m.ligand_name,
            "pose_rank": m.pose_rank,
            "engine_db": m.engine, "verdict_db": m.verdict,
            "wall_seconds_db": m.wall_seconds,
            "production_ps_db": m.production_ps,
            "snapshot_every_ps_db": m.snapshot_every_ps,
        } for m in db_runs}

    all_ids = sorted(set(alive_ids) | set(db_info.keys()))
    print(f"Total MD ids to consider: {len(all_ids)} "
          f"(alive on FS: {len(alive_ids)}, in DB: {len(db_info)})")

    rows: List[Dict[str, Any]] = []
    for i, mid in enumerate(all_ids, 1):
        print(f"  [{i:>2}/{len(all_ids)}] {mid[:8]} ...", end="", flush=True)
        try:
            row = analyze_one_md(mid, shared_pocket)
        except Exception as e:
            row = {"md_job_id": mid, "analyze_error": str(e),
                   "engine": "—", "verdict": "—"}
            print(f"  ERROR: {e}")
            continue
        # Splice DB info (compound name, pose_rank, lineage)
        db = db_info.get(mid) or {}
        row["docking_job_id"] = db.get("docking_job_id")
        row["ligand_name"] = db.get("ligand_name")
        row["pose_rank"] = db.get("pose_rank")
        row["compound"] = common_name_for(row["docking_job_id"],
                                          row["ligand_name"])
        # If summary engine was "—" but DB has one (failed runs may have
        # written the row but no summary), use DB info as the visible engine.
        if row["engine"] in ("—", None) and db.get("engine_db"):
            row["engine"] = db["engine_db"]
        if row["verdict"] in ("—", None) and db.get("verdict_db"):
            row["verdict"] = db["verdict_db"]
        if not row.get("duration_ps") and db.get("production_ps_db"):
            row["duration_ps"] = db["production_ps_db"]
        if not row.get("snapshot_every_ps") and db.get("snapshot_every_ps_db"):
            row["snapshot_every_ps"] = db["snapshot_every_ps_db"]
        rows.append(row)
        print(f"  engine={row.get('engine'):<13} "
              f"verdict={str(row.get('verdict')):<10} "
              f"frames={row.get('n_frames_on_disk'):>3}  "
              f"pose_pkt_final={row.get('pose_pkt_final_a')}")

    # Compose Notes per row
    def make_note(r: Dict[str, Any]) -> str:
        notes = []
        if r.get("docking_job_id") and str(r["docking_job_id"]).startswith("d96d719a"):
            notes.append("WRONG POCKET (d96d719a lineage; do not use for binding analysis)")
        if not r.get("summary_present"):
            notes.append("no summary.json (DB-only / never started)")
        if r.get("engine") == "surrogate":
            notes.append("RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction)")
        if r.get("n_frames_on_disk", 0) > 0 and r.get("n_frames_on_disk", 0) < 11:
            notes.append("very short (< 11 frames)")
        if r.get("engine") == "openmm_full" and r.get("n_frames_on_disk", 0) >= 11 and r.get("n_frames_on_disk", 0) < 50:
            notes.append("smoke-test scale")
        if r.get("duration_ps") and r.get("duration_ps") < 50:
            notes.append(f"production_ps={r['duration_ps']}")
        if r.get("converged") == "yes":
            notes.append("plateaued")
        elif r.get("converged") == "no":
            notes.append("still drifting at end")
        return "; ".join(notes) if notes else ""

    for r in rows:
        r["notes"] = make_note(r)

    # Group: taxol on correct pocket → controls → wrong-pocket → DB-only / no-summary
    def group_key(r: Dict[str, Any]) -> Tuple[int, str]:
        comp = r.get("compound") or ""
        dock = str(r.get("docking_job_id") or "")
        if dock.startswith("d96d719a") or "wrong-pocket" in comp:
            return (3, comp)
        if not r.get("summary_present"):
            return (4, comp)
        if comp == "taxol":
            # Order taxol by pose_rank (None last)
            return (0, str(r.get("pose_rank") if r.get("pose_rank") is not None else 99))
        if comp in ("Juliprosopine", "Primaquine"):
            return (1, comp)
        return (2, comp)
    rows.sort(key=group_key)

    # Write CSV
    csv_path = DOCS_ROOT / "all_md_diagnostics.csv"
    cols = [
        "md_job_id", "compound", "ligand_name", "pose_rank",
        "docking_job_id", "engine", "verdict", "solvent",
        "duration_ps", "snapshot_every_ps", "n_frames_on_disk",
        "pose_pkt_final_a", "pose_pkt_max_a",
        "com_displ_final_a", "com_displ_max_a",
        "first_half_mean_a", "second_half_mean_a", "slope_a_per_ps",
        "converged", "equilibration_frame_idx",
        "in_pocket_final", "buried_frac_final",
        "contacts_n_initial", "contacts_retained_frac_mean_eq",
        "summary_hbond_persistence_frac",
        "wall_seconds", "notes",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])
    print(f"\nWrote {csv_path}")

    # Write Markdown
    md_path = DOCS_ROOT / "all_md_diagnostics.md"

    def fmt(v):
        if v is None or v == "":
            return "—"
        if isinstance(v, bool):
            return "Y" if v else "N"
        if isinstance(v, float):
            if abs(v) < 1e-4 and v != 0:
                return f"{v:.1e}"
            return f"{v:.2f}"
        return str(v)

    def fmt_pair(r, a, b):
        return f"{fmt(r.get(a))} / {fmt(r.get(b))}"

    lines = [
        "# All MD diagnostics — every run ever queued (READ-ONLY)",
        "",
        f"Sweep of the full MD catalogue (DB md_runs ∪ jobs/<id>/md/) at this "
        f"branch state. {len(rows)} MD jobs total. Pocket residue set "
        f"derived ONCE from the canonical reference run "
        f"({REFERENCE_RUN_ID}, explicit taxol pose 0) frame 0 — "
        f"{len(shared_pocket)} residues within 5 Å of any ligand heavy atom — "
        f"and intersected with each run's own receptor for the deep "
        f"per-run analysis.",
        "",
        "Honest reporting: every run is listed, including failed / empty / "
        "surrogate / smoke-test / wrong-pocket. Caveat per row in the Notes "
        "column. `Converged?` = yes when |1st-half − 2nd-half mean| ≤ 0.3 Å "
        "AND |last-half slope| ≤ 1e-3 Å/ps; ≈ if borderline; no otherwise.",
        "",
        "**FE gate verdicts are deliberately NOT applied here.** Just the raw "
        "metrics, grouped for readability.",
        "",
        "## Table",
        "",
        "| job_id (8) | compound | engine | sol | dur(ps) | frames "
        "| pose-pkt fin/max | COM fin/max | conv | hb/contact persist "
        "| in-pkt fin | buried fin | notes |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]

    for r in rows:
        # hb persistence — prefer the deep "contacts_retained_frac_mean_eq"
        # for the long openmm_full runs; fall back to the engine's own
        # summary hbond_persistence_frac for shorter / surrogate runs so
        # the column is never blank when something WAS recorded.
        if r.get("contacts_retained_frac_mean_eq") is not None:
            persist = f"{r['contacts_retained_frac_mean_eq']*100:.0f}%  "\
                      f"(of {r.get('contacts_n_initial') or '?'} init contacts)"
        elif r.get("summary_hbond_persistence_frac") is not None:
            persist = f"hbond {r['summary_hbond_persistence_frac']*100:.0f}%"
        else:
            persist = "—"
        lines.append(
            "| " + " | ".join([
                f"`{r['md_job_id'][:8]}`",
                str(r.get("compound") or "—"),
                str(r.get("engine") or "—"),
                str(r.get("solvent") or "—"),
                fmt(r.get("duration_ps")),
                fmt(r.get("n_frames_on_disk")),
                fmt_pair(r, "pose_pkt_final_a", "pose_pkt_max_a"),
                fmt_pair(r, "com_displ_final_a", "com_displ_max_a"),
                str(r.get("converged") or "—"),
                persist,
                fmt(r.get("in_pocket_final")),
                fmt(r.get("buried_frac_final")),
                r.get("notes") or "—",
            ]) + " |"
        )

    lines += [
        "",
        "## Columns",
        "",
        "- **engine**: `openmm_full` (real MD), `surrogate` (RDKit MMFF "
        "wiggle in rigid receptor — pose-RMSD ≡ internal-RMSD by "
        "construction), or `—` (DB row only, no run output).",
        "- **sol**: implicit (OBC2) or explicit (TIP3P NPT octahedron per "
        "the [B7] / [B9] protocol).",
        "- **pose-pkt fin/max**: ligand pose RMSD in pocket-Cα-aligned "
        "frame (5 Å pocket Cα fit), final + trajectory max (Å). The "
        "primary binding-mode drift metric.",
        "- **COM fin/max**: ligand center-of-mass displacement after pocket "
        "alignment, final + max (Å).",
        "- **conv**: yes / ≈ / no per the convergence rule above.",
        "- **hb/contact persist**: for openmm_full runs with frames, the "
        "fraction of the frame-0 (ligand × protein heavy atom ≤ 4 Å) "
        "contacts retained per frame, mean over the equilibrated window "
        "(last 30 %). For surrogate / short runs, the engine's own "
        "`hbond_persistence_frac` is reported instead (different metric).",
        "- **in-pkt fin**: Y if the ligand COM in the last frame is < 5 Å "
        "from frame-0 COM after pocket alignment.",
        "- **buried fin**: fraction of ligand heavy atoms with any "
        "protein heavy atom within 4.5 Å, in the last frame.",
        "- **notes**: short caveat list; full per-row JSON in "
        "`docs/all_md_diagnostics.csv`.",
        "",
        "## Grouping",
        "",
        "Order: taxol on the verified taxane pocket (4a37bb0c), grouped by "
        "pose_rank; controls (Juliprosopine, Primaquine); other / "
        "unidentified runs; wrong-pocket lineage (d96d719a — receptor is a "
        "different construct, do NOT use for any binding analysis); "
        "DB-only rows (md_run entry exists but no summary.json on disk; "
        "typically a queued job that crashed before any output).",
    ]

    md_path.write_text("\n".join(lines))
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
