"""Apply the FINAL FE-gate criteria to every openmm_full MD run.

Simplified 2026-06-16 — C1 collapses to a single eq-window begin-vs-end
check. Slope (a) and quartile trend (c) are gone; (b) was already
measuring the same drift more transparently. The full gate is now four
criteria — converged, in-pocket, interactions, eq-frames — and a run
QUALIFIES iff all four PASS:

  1. CONVERGED   — eq-window beginning vs end mean pocket-RMSD differ by
                    < 0.4 Å. eq_begin = mean of the first 10 % of the
                    eq-window frames (last 30 % of the trajectory);
                    eq_end = mean of the last 10 %. Asks: once the ligand
                    entered its settled window, did RMSD stop changing?
  [2. SELF-STABLE — REMOVED.]
  3. IN-POCKET   — in-pocket fraction in eq window ≥ 0.95 AND COM
                    displacement at the final frame ≤ 3.0 Å.
  4. INTERACTIONS — top-5 residue contact persistence (eq window).
                    Top-5 = the 5 residues with the highest per-frame
                    contact fraction across the whole trajectory; per
                    frame, count how many of those 5 are in contact
                    (any-heavy-atom ≤ 4 Å); metric = mean over the eq
                    window of (count ÷ 5) ≥ 0.65 (65 %).
  5. EQ-WINDOW FRAMES — ≥ 50 frames in the last 30 % of the trajectory.

summary.verdict is still ignored. Implicit-solvent runs stay in the
table flagged as "implicit lineage — declined strategy".
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

from stability_diagnostics import (  # type: ignore
    parse_pdb_frame, kabsch_rotation, POCKET_RADIUS_A,
)

DOCS = Path(__file__).resolve().parents[1] / "docs"
DOCS.mkdir(parents=True, exist_ok=True)
REFERENCE_RUN = "5bc61f59-834f-4e71-a492-d32ddfdc7326"

# Same compound mapping used elsewhere in this repo
DOCKING_TO_LIGAND_TO_COMMON = {
    "4a37bb0c-7655-411e-a238-1ff5b7bab910": {
        "lig_000000": "taxol",
        "lig_000001": "Juliprosopine",
        "lig_000002": "Primaquine",
    },
    "d96d719a": {"lig_000000": "taxol (wrong-pocket)"},
}

# === Thresholds ===
EQ_BEGIN_END_DELTA_A       = 0.4      # C1: |eq_end - eq_begin| < this
EQ_BEGIN_END_FRAC          = 0.10     # eq-window slice (first/last 10 %)
IN_POCKET_FRAC_MIN         = 0.95
COM_FINAL_MAX_A            = 3.0
INTERACTIONS_PERSIST_MIN   = 0.65
CONTACT_CUTOFF_A           = 4.0       # heavy-atom contact distance
EQ_WINDOW_LAST_FRAC        = 0.30
EQ_WINDOW_MIN_FRAMES       = 50


def compound_for(docking_job_id, ligand_name) -> str:
    if docking_job_id is None or ligand_name is None:
        return ligand_name or "?"
    key = str(docking_job_id)
    for dock_id, lig_map in DOCKING_TO_LIGAND_TO_COMMON.items():
        if key.startswith(dock_id) or dock_id.startswith(key):
            return lig_map.get(ligand_name, ligand_name)
    return ligand_name


def derive_pocket_residues():
    """Single canonical reference — 5 Å pocket of frame 0 of REFERENCE_RUN."""
    ref_dir = JOBS_DIR / REFERENCE_RUN / "md"
    f0 = parse_pdb_frame(ref_dir / "frames" / "frame_000.pdb")
    lig = np.asarray([a["xyz"] for a in f0.ligand_heavy])
    rec_keys = []
    rec_pos = []
    for a in f0.receptor_heavy:
        rec_keys.append((a["chain"], a["resseq"]))
        rec_pos.append(a["xyz"])
    rec_pos = np.asarray(rec_pos)
    dmat = np.sqrt(((lig[:, None, :] - rec_pos[None, :, :]) ** 2).sum(axis=2))
    close = dmat.min(axis=0) < POCKET_RADIUS_A
    return {rec_keys[i] for i, v in enumerate(close) if v}


def evaluate_c1(eq_begin: float, eq_end: float) -> Tuple[bool, float, str]:
    """C1 — eq-window beginning vs end mean pocket-RMSD.

    Returns (passes, |eq_end - eq_begin|, "PASS" / "FAIL").
    """
    if not np.isfinite(eq_begin) or not np.isfinite(eq_end):
        return (False, float("nan"), "FAIL")
    diff = abs(eq_end - eq_begin)
    ok = diff < EQ_BEGIN_END_DELTA_A
    return (ok, diff, "PASS" if ok else "FAIL")


# ────────────────────────────────────────────────────────────────────────
def analyze_one(mid: str, shared_pocket: set) -> Optional[Dict[str, Any]]:
    md_dir = JOBS_DIR / mid / "md"
    if not (md_dir / "summary.json").is_file():
        return {"md_job_id": mid, "skipped": "no summary.json"}
    summary = json.loads((md_dir / "summary.json").read_text())
    engine = (summary.get("engine") or {}).get("kind")
    if engine != "openmm_full":
        return {"md_job_id": mid, "skipped": f"engine={engine}"}
    frame_paths = sorted((md_dir / "frames").glob("frame_*.pdb"))
    n = len(frame_paths)
    if n < 11:
        return {"md_job_id": mid, "skipped": f"only {n} frames"}

    settings = summary.get("settings") or {}
    snap_ps = float(settings.get("snapshot_every_ps") or 5.0)
    duration_ps = settings.get("production_ps")
    solvent = (settings.get("solvent") or "").lower() or None

    # ---- per-frame metrics (pocket-Cα Kabsch fit) ----
    snapshots = [parse_pdb_frame(p) for p in frame_paths]
    f0 = snapshots[0]
    rec0_ca = f0.receptor_ca
    rec0_keys = {(a["chain"], a["resseq"]) for a in rec0_ca}
    pocket_keys = shared_pocket & rec0_keys
    pocket_order = [(a["chain"], a["resseq"]) for a in rec0_ca
                    if (a["chain"], a["resseq"]) in pocket_keys]
    if not pocket_order:
        return {"md_job_id": mid, "skipped": "no pocket residues match canonical set"}
    pocket0_pos = np.asarray([
        next(a["xyz"] for a in rec0_ca if (a["chain"], a["resseq"]) == k)
        for k in pocket_order
    ])
    lig0 = np.asarray([a["xyz"] for a in f0.ligand_heavy])
    rec0_heavy_pos = np.asarray([a["xyz"] for a in f0.receptor_heavy])
    com0 = lig0.mean(axis=0)

    pose_pkt: List[float] = []
    lig_xfm_series: List[np.ndarray] = []
    rec_xfm_series: List[np.ndarray] = []
    com_series: List[float] = []
    in_pocket: List[bool] = []
    for snap in snapshots:
        rec_ca_pos = np.asarray([a["xyz"] for a in snap.receptor_ca])
        lig = np.asarray([a["xyz"] for a in snap.ligand_heavy])
        rec_heavy = np.asarray([a["xyz"] for a in snap.receptor_heavy])
        if (rec_ca_pos.shape != np.asarray([a["xyz"] for a in rec0_ca]).shape
                or lig.shape != lig0.shape
                or rec_heavy.shape != rec0_heavy_pos.shape):
            pose_pkt.append(np.nan)
            lig_xfm_series.append(np.full_like(lig0, np.nan))
            rec_xfm_series.append(np.full_like(rec0_heavy_pos, np.nan))
            com_series.append(np.nan); in_pocket.append(False); continue
        pocket_pos_i = np.asarray([
            next(a["xyz"] for a in snap.receptor_ca
                 if (a["chain"], a["resseq"]) == k)
            for k in pocket_order
        ])
        p_cent = pocket_pos_i.mean(axis=0); q_cent = pocket0_pos.mean(axis=0)
        P = pocket_pos_i - p_cent; Q = pocket0_pos - q_cent
        R = kabsch_rotation(P, Q)
        lig_xfm = (lig - p_cent) @ R.T + q_cent
        rec_xfm = (rec_heavy - p_cent) @ R.T + q_cent
        lig_xfm_series.append(lig_xfm)
        rec_xfm_series.append(rec_xfm)
        pose_pkt.append(float(np.sqrt(((lig_xfm - lig0) ** 2).sum(axis=1).mean())))
        com = lig_xfm.mean(axis=0)
        d = float(np.linalg.norm(com - com0))
        com_series.append(d); in_pocket.append(d < 5.0)

    arr = np.asarray(pose_pkt)

    # ---- equilibrated window ----
    eq_start = int(round(n * (1.0 - EQ_WINDOW_LAST_FRAC)))
    eq_count = n - eq_start

    # ---- eq-window BEGINNING vs END means (C1) ----
    slice_count = max(1, int(round(eq_count * EQ_BEGIN_END_FRAC)))
    eq_begin_idx = list(range(eq_start, min(eq_start + slice_count, n)))
    eq_end_idx = list(range(max(eq_start, n - slice_count), n))
    eq_begin_mean = (float(np.nanmean(arr[eq_begin_idx]))
                     if eq_begin_idx else float("nan"))
    eq_end_mean = (float(np.nanmean(arr[eq_end_idx]))
                   if eq_end_idx else float("nan"))

    in_pocket_frac_eq = float(np.mean(in_pocket[eq_start:])) if eq_count > 0 else float("nan")
    com_final = float(com_series[-1]) if com_series else float("nan")

    # ---- INTERACTIONS (top-5 residue contact persistence over eq window) ----
    # Step 1 — map each receptor heavy atom to its residue:
    #   residue_atoms[(chain, resseq)] = (list_of_atom_indices, resname)
    residue_atoms: Dict[Tuple[str, int], Tuple[List[int], str]] = {}
    for atom_idx, a in enumerate(f0.receptor_heavy):
        k = (a["chain"], a["resseq"])
        rn = a.get("resname") or "?"
        if k not in residue_atoms:
            residue_atoms[k] = ([], rn)
        residue_atoms[k][0].append(atom_idx)
    res_keys = list(residue_atoms.keys())
    n_res = len(res_keys)

    # Step 2 — per-frame per-residue "in contact" boolean
    #   (any ligand heavy atom ≤ 4 Å of any of the residue's heavy atoms)
    res_contact_mat = np.zeros((n, n_res), dtype=bool)
    for f_idx in range(n):
        L = lig_xfm_series[f_idx]; Rcv = rec_xfm_series[f_idx]
        if np.isnan(L).any() or np.isnan(Rcv).any():
            continue
        # Per-receptor-atom: closest distance to any ligand heavy atom
        diffs = L[:, None, :] - Rcv[None, :, :]
        close_atom = (np.sqrt((diffs * diffs).sum(axis=2)) <= CONTACT_CUTOFF_A)  # (lig, rec)
        any_lig_close = close_atom.any(axis=0)  # (rec,)
        for r_idx, k in enumerate(res_keys):
            atom_idx_list, _ = residue_atoms[k]
            res_contact_mat[f_idx, r_idx] = bool(any_lig_close[atom_idx_list].any())

    # Step 3 — trajectory-wide per-residue contact frequency, pick top-5
    res_freq = res_contact_mat.mean(axis=0)  # over ALL frames
    order = np.argsort(-res_freq)
    top5_idx = order[:5]
    top5 = [
        {
            "chain": res_keys[i][0],
            "resseq": int(res_keys[i][1]),
            "resname": residue_atoms[res_keys[i]][1],
            "trajectory_frac": float(res_freq[i]),
        }
        for i in top5_idx
    ]

    # Step 4 — eq-window mean of (#-of-top-5-in-contact / 5)
    top5_per_frame = res_contact_mat[:, top5_idx].sum(axis=1) / 5.0
    top5_persistence_eq = (float(top5_per_frame[eq_start:].mean())
                           if eq_count > 0 else float("nan"))
    n_init = int(res_freq.sum() > 0)  # diagnostic only (kept for column compat)

    # ---- 4-criterion gate (C1 = eq-begin vs eq-end; C2 removed) ----
    c1_converged, eq_delta, c1_label = evaluate_c1(eq_begin_mean, eq_end_mean)

    c3_in_pocket = (
        np.isfinite(in_pocket_frac_eq) and in_pocket_frac_eq >= IN_POCKET_FRAC_MIN
        and np.isfinite(com_final) and com_final <= COM_FINAL_MAX_A
    )
    c4_interactions = (
        np.isfinite(top5_persistence_eq)
        and top5_persistence_eq >= INTERACTIONS_PERSIST_MIN
    )
    c5_eq_frames = eq_count >= EQ_WINDOW_MIN_FRAMES
    qualifies = bool(c1_converged and c3_in_pocket
                     and c4_interactions and c5_eq_frames)

    return {
        "md_job_id": mid,
        "engine": engine,
        "solvent": solvent,
        "verdict_ignored": summary.get("verdict"),
        "n_frames": n,
        "duration_ps": duration_ps,
        "snapshot_every_ps": snap_ps,
        "eq_window_count": eq_count,
        "eq_slice_count": slice_count,
        "eq_begin_mean_a": eq_begin_mean,
        "eq_end_mean_a": eq_end_mean,
        "eq_begin_end_delta_a": eq_delta if np.isfinite(eq_delta) else None,
        "c1_label": c1_label,
        "in_pocket_frac_eq": in_pocket_frac_eq,
        "com_final_a": com_final,
        "top5_persistence_eq": top5_persistence_eq,
        "top5_residues": top5,
        "c1_converged": c1_converged,
        "c3_in_pocket": c3_in_pocket,
        "c4_interactions": c4_interactions,
        "c5_eq_frames": c5_eq_frames,
        "qualifies": qualifies,
    }


def main():
    print("=== Corrected FE-gate evaluator (revised 2026-06-16) ===\n")
    shared_pocket = derive_pocket_residues()
    print(f"Canonical pocket residue set: {len(shared_pocket)} residues "
          f"(5 Å of ligand in {REFERENCE_RUN[:8]} frame 0)\n")

    with SessionLocal() as db:
        db_runs = db.query(MDRun).all()
        db_info = {str(m.id): m for m in db_runs}
    alive = {p.name for p in JOBS_DIR.iterdir()
             if p.is_dir() and p.name != "archive"
             and (p / "md" / "summary.json").is_file()}
    candidates = []
    for mid in sorted(set(alive) | set(db_info.keys())):
        sp = JOBS_DIR / mid / "md" / "summary.json"
        engine_fs = None
        if sp.is_file():
            try:
                engine_fs = (json.loads(sp.read_text()).get("engine") or {}).get("kind")
            except Exception:
                pass
        engine_db = db_info.get(mid).engine if mid in db_info else None
        if (engine_fs == "openmm_full") or (engine_db == "openmm_full"):
            candidates.append(mid)
    print(f"openmm_full candidates: {len(candidates)}")

    rows = []
    for i, mid in enumerate(candidates, 1):
        print(f"  [{i:>2}/{len(candidates)}] {mid[:8]} ...", end="", flush=True)
        result = analyze_one(mid, shared_pocket)
        m = db_info.get(mid)
        compound = compound_for(m.docking_job_id, m.ligand_name) if m else "?"
        result.update({
            "compound": compound,
            "ligand_name": m.ligand_name if m else None,
            "pose_rank": m.pose_rank if m else None,
            "docking_job_id": str(m.docking_job_id) if m and m.docking_job_id else None,
        })
        if "skipped" in result:
            print(f"  skipped ({result['skipped']})")
        else:
            tag = "QUALIFIES" if result["qualifies"] else "no"
            print(f"  {tag}")
        rows.append(result)

    def sort_key(r):
        qual = 0 if r.get("qualifies") else 1
        comp = r.get("compound") or ""
        comp_rank = {"taxol": 0, "Juliprosopine": 1, "Primaquine": 2}.get(comp, 5)
        if "wrong-pocket" in comp: comp_rank = 9
        pose = r.get("pose_rank") if r.get("pose_rank") is not None else 99
        return (qual, comp_rank, pose, r["md_job_id"])
    rows.sort(key=sort_key)

    # === CSV ===
    cols = [
        "md_job_id", "compound", "ligand_name", "pose_rank", "engine", "solvent",
        "verdict_ignored", "duration_ps", "n_frames",
        "eq_window_count", "eq_slice_count",
        "eq_begin_mean_a", "eq_end_mean_a", "eq_begin_end_delta_a", "c1_label",
        "in_pocket_frac_eq", "com_final_a",
        "top5_persistence_eq", "top5_residues",
        "c1_converged", "c3_in_pocket",
        "c4_interactions", "c5_eq_frames", "qualifies",
    ]
    csv_path = DOCS / "all_md_corrected_gate.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f); w.writerow(cols)
        for r in rows:
            row_vals = []
            for c in cols:
                v = r.get(c, "")
                if c == "top5_residues" and isinstance(v, list):
                    v = ";".join(f"{x['resname']} {x['chain']}:{x['resseq']}@"
                                 f"{x['trajectory_frac']:.2f}" for x in v)
                row_vals.append(v)
            w.writerow(row_vals)
    print(f"\nWrote {csv_path}")

    # === MARKDOWN ===
    def fmt(v):
        if v is None or v == "":
            return "—"
        if isinstance(v, bool):
            return "✅" if v else "❌"
        if isinstance(v, float):
            if abs(v) < 1e-4 and v != 0:
                return f"{v:.1e}"
            return f"{v:.2f}"
        return str(v)

    def yn(b):
        if b is None: return "—"
        return "✅" if b else "❌"

    def notes_for(r):
        parts = []
        if (r.get("solvent") or "").lower() == "implicit":
            parts.append("**implicit lineage** — declined strategy "
                         "(explicit TIP3P is the only supported FE-feeder; "
                         "kept for evidence, not an FE candidate)")
        if r.get("verdict_ignored"):
            parts.append(f"verdict.json said `{r['verdict_ignored']}` (ignored)")
        return " · ".join(parts) if parts else "—"

    evaluated = [r for r in rows if "skipped" not in r]
    qual_runs = [r for r in evaluated if r.get("qualifies")]

    lines: List[str] = []
    lines += [
        "# FE gate (final 2026-06-16) — every openmm_full MD run",
        "",
        "C1 collapses to a single eq-window begin-vs-end check — slope (a) "
        "and visual trend (c) are gone; both were measuring the same drift "
        "less directly than the begin/end means. The full gate is now 4 "
        "criteria. C2 (self-stable) remains removed (redundant). "
        "`summary.verdict` ignored; implicit-solvent runs flagged.",
        "",
        "## Criteria",
        "",
        f"**1. CONVERGED** — **|eq_end_mean − eq_begin_mean| < "
        f"{EQ_BEGIN_END_DELTA_A} Å**. eq_begin = mean pocket-aligned RMSD over "
        f"the first **{int(EQ_BEGIN_END_FRAC*100)} %** of the eq-window frames; "
        f"eq_end = mean over the last **{int(EQ_BEGIN_END_FRAC*100)} %**. "
        "For a 60-frame eq window that's the first and last 6 frames. "
        "Asks: once the ligand entered its settled window, did RMSD stop "
        "changing?",
        "",
        "**2. SELF-STABLE — REMOVED** (redundant once C1 + C3 + C4 hold).",
        "",
        f"**3. IN-POCKET** — in-pocket fraction in eq window ≥ "
        f"**{IN_POCKET_FRAC_MIN}** AND COM displacement at the final frame ≤ "
        f"**{COM_FINAL_MAX_A} Å**.",
        f"**4. INTERACTIONS** — **top-5 residue contact persistence** "
        "(≥ **65 %**). The top-5 are the 5 protein residues with the "
        "highest per-frame contact fraction across the whole trajectory "
        f"(any-heavy-atom ≤ **{CONTACT_CUTOFF_A} Å**). Per eq-window "
        "frame: count how many of those 5 are in contact. C4 metric = "
        "mean of (count ÷ 5). Captures “key interactions held” while "
        "tolerating the all-pair reshuffle.",
        f"**5. EQ-WINDOW FRAMES** — ≥ **{EQ_WINDOW_MIN_FRAMES}** frames in "
        f"the last {int(EQ_WINDOW_LAST_FRAC*100)} % of the trajectory.",
        "",
        "## Result",
        "",
        f"**{len(qual_runs)} of {len(evaluated)} openmm_full runs QUALIFY "
        "for FE computation.**",
        "",
        "**Solvent-strategy note (2026-06-16):** explicit TIP3P is the only "
        "supported FE-feeder. Implicit-solvent OBC2 runs are kept in this "
        "table as historical evidence — even if they qualify on the metric "
        "criteria they are NOT FE candidates; the **Notes** column flags "
        "them as `implicit lineage — declined strategy`.",
        "",
        "## Summary table",
        "",
        "| Job ID (8) | Compound | Sol | Δ(begin−end) Å | C3 in-pocket | "
        "C4 top-5 % | C5 eq-frames | **QUALIFIES** | Notes |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    for r in rows:
        if "skipped" in r:
            lines.append("| " + " | ".join([
                f"`{r['md_job_id'][:8]}`",
                str(r.get("compound") or "—"),
                "—", "—", "—", "—", "—", "—",
                f"skipped ({r['skipped']})",
            ]) + " |")
            continue
        delta = r.get("eq_begin_end_delta_a")
        delta_cell = (f"{yn(r['c1_converged'])} ({delta:.2f} Å)"
                      if delta is not None and np.isfinite(delta)
                      else f"{yn(r['c1_converged'])} (—)")
        c3_cell = f"{yn(r['c3_in_pocket'])} (frac {fmt(r['in_pocket_frac_eq'])}, COM {fmt(r['com_final_a'])} Å)"
        c4_cell = f"{yn(r['c4_interactions'])} ({fmt(r['top5_persistence_eq']*100 if r['top5_persistence_eq'] is not None and np.isfinite(r['top5_persistence_eq']) else float('nan'))}%)"
        c5_cell = f"{yn(r['c5_eq_frames'])} ({r['eq_window_count']}/{EQ_WINDOW_MIN_FRAMES})"
        lines.append("| " + " | ".join([
            f"`{r['md_job_id'][:8]}`",
            str(r.get("compound") or "—"),
            (r.get("solvent") or "—"),
            delta_cell, c3_cell, c4_cell, c5_cell,
            ("✅ **YES**" if r["qualifies"] else "❌ no"),
            notes_for(r),
        ]) + " |")

    # === Detailed per-run blocks ===
    lines += ["", "## Detailed per-run blocks", ""]
    for r in rows:
        if "skipped" in r:
            lines += [
                f"### `{r['md_job_id'][:8]}` — {r.get('compound') or '?'}  "
                f"_(skipped: {r['skipped']})_",
                "",
            ]
            continue
        eq_delta = r["eq_begin_end_delta_a"]
        eq_delta_str = (f"{eq_delta:.2f} Å" if eq_delta is not None
                        and np.isfinite(eq_delta) else "—")
        eq_begin_str = (f"{r['eq_begin_mean_a']:.2f}"
                        if np.isfinite(r['eq_begin_mean_a']) else "—")
        eq_end_str = (f"{r['eq_end_mean_a']:.2f}"
                      if np.isfinite(r['eq_end_mean_a']) else "—")
        ret_pct = (r["top5_persistence_eq"] * 100.0
                   if r["top5_persistence_eq"] is not None
                   and np.isfinite(r["top5_persistence_eq"]) else float("nan"))
        top5_str = "; ".join(
            f"{x['resname']} {x['chain']}:{x['resseq']} ({x['trajectory_frac']*100:.0f}%)"
            for x in (r.get("top5_residues") or [])
        ) or "—"
        verdict_note = (f"verdict.json said `{r['verdict_ignored']}` (ignored)"
                        if r.get("verdict_ignored") else "—")
        sol_note = (("**implicit lineage** — declined strategy. "
                     if (r.get("solvent") or "").lower() == "implicit" else ""))

        lines += [
            f"### `{r['md_job_id']}` — {r.get('compound') or '?'} "
            f"(pose {r.get('pose_rank')}, solvent={r.get('solvent')}, "
            f"{fmt(r['duration_ps'])} ps, {r['n_frames']} frames)",
            "",
            "```",
            f"Job: {r['md_job_id']} | Compound: {r.get('compound')} | Engine: openmm_full",
            "─" * 74,
            f"C1 Converged: eq-window begin vs end "
            f"({r['eq_slice_count']} frames each): {eq_begin_str} vs "
            f"{eq_end_str} Å, Δ(begin−end) = {eq_delta_str} "
            f"(threshold < {EQ_BEGIN_END_DELTA_A} Å) "
            f"→ {r['c1_label']}",
            "",
            f"In-pocket: fraction {r['in_pocket_frac_eq']:.2f}, "
            f"COM {r['com_final_a']:.2f} Å "
            f"→ {'PASS' if r['c3_in_pocket'] else 'FAIL'}",
            f"Interactions: top-5 residue persistence {ret_pct:.1f}% "
            f"(threshold {int(INTERACTIONS_PERSIST_MIN*100)}%) "
            f"→ {'PASS' if r['c4_interactions'] else 'FAIL'}",
            f"  Top-5 residues (whole-trajectory contact frac): {top5_str}",
            f"Equilibrated frames: {r['eq_window_count']} "
            f"(threshold {EQ_WINDOW_MIN_FRAMES}) "
            f"→ {'PASS' if r['c5_eq_frames'] else 'FAIL'}",
            "",
            f"**QUALIFIES: {'YES' if r['qualifies'] else 'NO'}** | Notes: "
            f"{sol_note}{verdict_note}",
            "```",
            "",
        ]

    # === Summary footer ===
    lines += ["## QUALIFYING runs", ""]
    if qual_runs:
        for r in qual_runs:
            pose = r.get("pose_rank")
            pose_str = f"pose {pose}" if pose is not None else ""
            lines.append(
                f"- `{r['md_job_id']}` — **{r['compound']}** {pose_str}, "
                f"solvent={r.get('solvent')}, {fmt(r['duration_ps'])} ps, "
                f"{r['n_frames']} frames"
            )
    else:
        lines.append("_None._ Per-criterion blocker counts:")
        for c, label in (
            ("c1_converged", "C1 CONVERGED"),
            ("c3_in_pocket", "C3 IN-POCKET"),
            ("c4_interactions", "C4 INTERACTIONS"),
            ("c5_eq_frames", "C5 EQ-FRAMES"),
        ):
            n_fail = sum(1 for r in evaluated if not r[c])
            lines.append(f"  - {label}: {n_fail} / {len(evaluated)} fail")

    md_path = DOCS / "all_md_corrected_gate.md"
    md_path.write_text("\n".join(lines))
    print(f"Wrote {md_path}")
    if qual_runs:
        print(f"\nQUALIFYING runs ({len(qual_runs)}):")
        for r in qual_runs:
            print(f"  {r['md_job_id'][:8]}  {r['compound']}  pose={r.get('pose_rank')}  "
                  f"solvent={r.get('solvent')}  frames={r['n_frames']}")
    else:
        print("\nNo run qualifies.")


if __name__ == "__main__":
    main()
