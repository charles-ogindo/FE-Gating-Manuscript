"""
Post-run analysis for the [B9] explicit-solvent production MD.

Reads frames from jobs/<md_id>/md/frames/, computes:
  - Backbone Cα RMSD over time vs frame 0 (global)
  - Ligand pose RMSD with GLOBAL backbone-Cα alignment
  - Ligand pose RMSD with POCKET-LOCAL backbone-Cα alignment
    (pocket residues = residues with any Cα within 10 Å of the
     taxane-site box center in frame 0)

Saves figures as PNGs under md/figures/ and a summary JSON to
md/pocket_rmsd_summary.json. Independent of the existing analyze.py
flow (which already wrote rmsd.csv with global alignment).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Same _WATERS / _NON_LIGAND_HETATMS rules analyze.py uses, copied here to
# keep this script self-contained for downstream re-runs without import
# coupling to the MD package layout.
_WATERS = {"HOH", "WAT", "H2O", "TIP3", "TIP", "SOL"}
_METAL_IONS = {
    "MG", "ZN", "CA", "MN", "FE", "FE2", "FE3", "CO", "NI", "CU", "CD",
    "HG", "NA", "K", "CL", "LI", "RB", "CS", "SR", "BA",
}
_COFACTORS = {"GTP", "GDP", "ATP", "ADP", "AMP"}
_NON_LIGAND_HETATMS = _WATERS | _METAL_IONS | _COFACTORS

# Taxane-site box centroid (verified for docking job 4a37bb0c).
BOX_CENTER_A = (-0.399, -16.403, 14.621)
POCKET_RADIUS_A = 10.0


def parse_pdb_frame(path: Path):
    """Return (receptor_ca, ligand_heavy, all_receptor_atoms) — each a list
    of dicts {chain, resseq, resname, name, xyz, element}."""
    receptor_ca = []
    receptor_all = []
    ligand = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        tag = line[:6]
        if tag not in ("ATOM  ", "HETATM"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            resseq = 0
        element = (line[76:78].strip() if len(line) >= 78 else "")
        atom = {"chain": chain, "resseq": resseq, "resname": resname,
                "name": name, "xyz": (x, y, z), "element": element}
        if tag == "ATOM  ":
            receptor_all.append(atom)
            if name == "CA":
                receptor_ca.append(atom)
        elif resname in _NON_LIGAND_HETATMS:
            continue
        else:
            if element.upper() != "H":
                ligand.append(atom)
    return receptor_ca, ligand, receptor_all


def kabsch_rotation(P, Q):
    """Return the optimal rotation matrix R that aligns P onto Q (both N×3
    arrays of corresponding atoms), in the least-squares sense."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def rmsd_after_alignment(p_align, q_align, p_apply, q_apply):
    """Compute the RMSD between p_apply (transformed via the Kabsch fit of
    p_align → q_align) and q_apply. p_align/q_align: N×3 arrays defining the
    alignment; p_apply/q_apply: M×3 arrays whose RMSD we want post-fit."""
    p_centroid = p_align.mean(axis=0)
    q_centroid = q_align.mean(axis=0)
    P = p_align - p_centroid
    Q = q_align - q_centroid
    R = kabsch_rotation(P, Q)
    # Apply: translate by -p_centroid, rotate, translate by +q_centroid.
    pa = (p_apply - p_centroid) @ R.T + q_centroid
    return float(np.sqrt(((pa - q_apply) ** 2).sum(axis=1).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("md_dir", help="Path to jobs/<md_id>/md/")
    args = ap.parse_args()
    md_dir = Path(args.md_dir)
    frames_dir = md_dir / "frames"
    figures_dir = md_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    frame_paths = sorted(frames_dir.glob("frame_*.pdb"))
    n_frames = len(frame_paths)
    print(f"Found {n_frames} frames in {frames_dir}")
    assert n_frames > 0, "no frames!"

    # Parse all frames once.
    snapshots = [parse_pdb_frame(p) for p in frame_paths]
    rec0_ca, lig0, _ = snapshots[0]

    # snapshot index: time in ps (snapshot_every_ps from summary.json)
    summary = json.loads((md_dir / "summary.json").read_text())
    snap_ps = float(summary["settings"]["snapshot_every_ps"])
    times_ps = [i * snap_ps for i in range(n_frames)]

    # ----- Pocket residue selection (frame 0 only) -----
    bx, by, bz = BOX_CENTER_A
    pocket_keys = set()
    for ca in rec0_ca:
        x, y, z = ca["xyz"]
        d2 = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
        if d2 <= POCKET_RADIUS_A ** 2:
            pocket_keys.add((ca["chain"], ca["resseq"]))
    print(f"Pocket: {len(pocket_keys)} Cα residues within "
          f"{POCKET_RADIUS_A:.1f} Å of box center")

    # ----- Reference arrays -----
    rec0_ca_pos = np.asarray([ca["xyz"] for ca in rec0_ca])
    lig0_heavy = np.asarray([a["xyz"] for a in lig0])
    # Pocket reference: a sub-array aligned to a stable atom-key ordering
    pocket_order = [(ca["chain"], ca["resseq"]) for ca in rec0_ca
                    if (ca["chain"], ca["resseq"]) in pocket_keys]
    pocket0_pos = np.asarray([
        next(ca["xyz"] for ca in rec0_ca
             if (ca["chain"], ca["resseq"]) == k)
        for k in pocket_order
    ])
    print(f"Reference shapes: rec0_ca={rec0_ca_pos.shape}, "
          f"lig0={lig0_heavy.shape}, pocket0={pocket0_pos.shape}")

    # ----- Per-frame RMSDs -----
    rmsd_bb = []          # global Cα RMSD vs frame 0 (no alignment — raw drift)
    rmsd_bb_aligned = []  # global Cα RMSD post-global-alignment
    rmsd_pose_global = []
    rmsd_pose_pocket = []
    rmsd_lig_internal = []

    for i, (rec_ca, lig, _) in enumerate(snapshots):
        rec_ca_pos = np.asarray([ca["xyz"] for ca in rec_ca])
        lig_heavy = np.asarray([a["xyz"] for a in lig])
        if rec_ca_pos.shape != rec0_ca_pos.shape or lig_heavy.shape != lig0_heavy.shape:
            # Skip frames whose atom-set doesn't match (truncation drift / read error).
            rmsd_bb.append(float("nan"))
            rmsd_bb_aligned.append(float("nan"))
            rmsd_pose_global.append(float("nan"))
            rmsd_pose_pocket.append(float("nan"))
            rmsd_lig_internal.append(float("nan"))
            continue
        # Raw backbone drift (no alignment — measures translation + rotation)
        rmsd_bb.append(float(np.sqrt(((rec_ca_pos - rec0_ca_pos) ** 2).sum(axis=1).mean())))
        # Backbone RMSD post-alignment (residual structural drift)
        rmsd_bb_aligned.append(
            rmsd_after_alignment(rec_ca_pos, rec0_ca_pos, rec_ca_pos, rec0_ca_pos)
        )
        # Pose RMSD with GLOBAL backbone alignment
        rmsd_pose_global.append(
            rmsd_after_alignment(rec_ca_pos, rec0_ca_pos, lig_heavy, lig0_heavy)
        )
        # Pose RMSD with POCKET-LOCAL backbone alignment
        pocket_pos_i = np.asarray([
            next(ca["xyz"] for ca in rec_ca
                 if (ca["chain"], ca["resseq"]) == k)
            for k in pocket_order
        ])
        rmsd_pose_pocket.append(
            rmsd_after_alignment(pocket_pos_i, pocket0_pos, lig_heavy, lig0_heavy)
        )
        # Ligand-internal Kabsch RMSD (frame i ligand aligned to frame 0 ligand)
        rmsd_lig_internal.append(
            rmsd_after_alignment(lig_heavy, lig0_heavy, lig_heavy, lig0_heavy)
        )

    # ----- Summary numbers (mirror Run B's reporting style) -----
    metrics = {
        "n_frames": n_frames,
        "snapshot_every_ps": snap_ps,
        "duration_ps": (n_frames - 1) * snap_ps,
        "pocket_residue_count": len(pocket_keys),
        "rmsd_backbone_final_a": rmsd_bb_aligned[-1],
        "rmsd_backbone_max_a": float(np.nanmax(rmsd_bb_aligned)),
        "rmsd_pose_global_final_a": rmsd_pose_global[-1],
        "rmsd_pose_global_max_a": float(np.nanmax(rmsd_pose_global)),
        "rmsd_pose_pocket_final_a": rmsd_pose_pocket[-1],
        "rmsd_pose_pocket_max_a": float(np.nanmax(rmsd_pose_pocket)),
        "rmsd_ligand_internal_final_a": rmsd_lig_internal[-1],
        "rmsd_ligand_internal_max_a": float(np.nanmax(rmsd_lig_internal)),
    }
    out_json = md_dir / "pocket_rmsd_summary.json"
    out_json.write_text(json.dumps(metrics, indent=2))
    print(f"Wrote {out_json.name}")

    # ----- Figures -----
    # 1) RMSD over time: backbone + ligand pose (global) + ligand pose (pocket)
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(times_ps, rmsd_bb_aligned, "-", color="C0",
            label=f"backbone Cα (global align)  final {rmsd_bb_aligned[-1]:.2f}")
    ax.plot(times_ps, rmsd_pose_global, "-", color="C1",
            label=f"ligand pose (global align)   final {rmsd_pose_global[-1]:.2f}  "
                  f"max {metrics['rmsd_pose_global_max_a']:.2f}")
    ax.plot(times_ps, rmsd_pose_pocket, "-", color="C2",
            label=f"ligand pose (pocket align)   final {rmsd_pose_pocket[-1]:.2f}  "
                  f"max {metrics['rmsd_pose_pocket_max_a']:.2f}")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title(f"[B9] explicit-solvent MD — RMSD vs frame 0  "
                 f"({n_frames} frames, {metrics['duration_ps']:.0f} ps)")
    ax.legend(loc="best", frameon=False)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig_path = figures_dir / "rmsd_global_vs_pocket.png"
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)
    print(f"Wrote {fig_path.relative_to(md_dir)}")

    # 2) Pose vs internal (the Q6b distinction)
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(times_ps, rmsd_pose_global, "-", color="C1",
            label=f"pose RMSD (global)  final {rmsd_pose_global[-1]:.2f}")
    ax.plot(times_ps, rmsd_pose_pocket, "-", color="C2",
            label=f"pose RMSD (pocket)  final {rmsd_pose_pocket[-1]:.2f}")
    ax.plot(times_ps, rmsd_lig_internal, "--", color="C3",
            label=f"internal RMSD       final {rmsd_lig_internal[-1]:.2f}")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Pose vs internal ligand RMSD")
    ax.legend(loc="best", frameon=False)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig_path = figures_dir / "rmsd_pose_vs_internal.png"
    fig.savefig(fig_path, dpi=140)
    plt.close(fig)
    print(f"Wrote {fig_path.relative_to(md_dir)}")

    # Top-line numbers for the user
    print("\n=== Headline numbers ===")
    print(f"  pose global  final/max:  {metrics['rmsd_pose_global_final_a']:.2f} / "
          f"{metrics['rmsd_pose_global_max_a']:.2f} Å")
    print(f"  pose pocket  final/max:  {metrics['rmsd_pose_pocket_final_a']:.2f} / "
          f"{metrics['rmsd_pose_pocket_max_a']:.2f} Å")
    print(f"  backbone     final/max:  {metrics['rmsd_backbone_final_a']:.2f} / "
          f"{metrics['rmsd_backbone_max_a']:.2f} Å")
    print(f"  internal     final/max:  {metrics['rmsd_ligand_internal_final_a']:.2f} / "
          f"{metrics['rmsd_ligand_internal_max_a']:.2f} Å")


if __name__ == "__main__":
    main()
