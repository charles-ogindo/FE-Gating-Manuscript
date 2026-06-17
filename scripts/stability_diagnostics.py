"""
Read-only diagnostic measurement of ligand binding-mode stability across MD
runs. Designed for the explicit-solvent vs Run A vs Run B comparison of
taxol at the taxane site (docking 4a37bb0c-7655-411e-a238-1ff5b7bab910),
but the metrics are run-agnostic; the only taxol-specific assumption is
the core/side-chain mask derived from the docked-ligand SMILES via RDKit
ring-system analysis (largest aliphatic fused-ring system = baccatin core).

Metrics computed per run:

  1. Whole-ligand pose RMSD vs frame 0:
       - global Cα alignment
       - pocket-local Cα alignment (pocket = residues with any heavy atom
         within 5 Å of any ligand heavy atom in frame 0)
     final / max / mean over the equilibrated window (last 30% of frames,
     matching the free-energy gate's DEFAULT_LAST_FRACTION).

  2. Core (baccatin) RMSD + flexible side-chain RMSD (both pocket-aligned).
     Same final/max/mean.

  3. Pocket-local backbone RMSD (N, CA, C, O of pocket residues) + per-
     residue Cα RMSF over the equilibrated window. Whole-protein backbone
     RMSD for comparison.

  4. Contact-network persistence: every (ligand heavy atom, protein heavy
     atom) pair within 4.0 Å in frame 0 is an initial contact. Per-contact
     occupancy across all frames; mean fraction of initial contacts
     retained per frame across the equilibrated window.

  5. Pocket occupancy: ligand center-of-mass displacement from frame 0
     (final/max, pocket-aligned coords); fraction of frames with COM
     displacement < 5.0 Å (in-pocket); mean fraction of ligand heavy
     atoms with any protein heavy atom within 4.5 Å (buried).

  6. Convergence: block-averaged pose RMSD over the first half vs the
     second half; linear-fit slope of pose RMSD over the second half.

Frame matching: pocket-local alignment is computed once per frame and
used for ALL pocket-aligned metrics (ligand pose, core, side-chain,
pocket backbone) — so the comparison between metric 3's pocket BB RMSD
and metric 1's pose RMSD lives in the same alignment frame.

Output:
  - One CSV table (rows = runs, columns = metrics) at
    docs/B9_stability_diagnostics.csv
  - A markdown summary at docs/B9_stability_diagnostics.md
  - Per-figure PNGs (whole-ligand / core / side-chain / pocket-BB RMSD
    over time, all runs overlaid) at docs/B9_stability_*.png
  - A per-run JSON dump at docs/B9_stability_<run_label>.json
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BOX_CENTER_A = (-0.399, -16.403, 14.621)  # taxane site for docking 4a37bb0c
POCKET_RADIUS_A = 5.0          # metric 3 def: residues within 5 Å of ligand in frame 0
CONTACT_CUTOFF_A = 4.0         # any heavy atom pair within this = contact
BURIED_CUTOFF_A = 4.5          # ligand heavy atom with protein within this = buried
IN_POCKET_COM_A = 5.0          # ligand COM within this of frame 0 = in pocket
EQ_WINDOW_LAST_FRAC = 0.30     # last 30% = equilibrated window (gates default)

# Same-as-analyze.py resname exclusion sets
_WATERS = {"HOH", "WAT", "H2O", "TIP3", "TIP", "SOL"}
_METAL_IONS = {
    "MG", "ZN", "CA", "MN", "FE", "FE2", "FE3", "CO", "NI", "CU", "CD",
    "HG", "NA", "K", "CL", "LI", "RB", "CS", "SR", "BA",
}
_COFACTORS = {"GTP", "GDP", "ATP", "ADP", "AMP", "NME"}
_NON_LIGAND_HETATMS = _WATERS | _METAL_IONS | _COFACTORS


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

@dataclass
class FrameAtoms:
    receptor_ca: List[dict]    # [{chain, resseq, resname, name, xyz}, ...]
    receptor_bb: List[dict]    # backbone N/CA/C/O of all protein residues
    ligand_heavy: List[dict]   # ligand heavy atoms, name like "C1", "O1", "N1", ...
    receptor_heavy: List[dict] # all protein heavy atoms (for contact + buried calcs)


def parse_pdb_frame(path: Path) -> FrameAtoms:
    rec_ca: List[dict] = []
    rec_bb: List[dict] = []
    rec_heavy: List[dict] = []
    lig: List[dict] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
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
                if element.upper() != "H":
                    rec_heavy.append(atom)
                if name == "CA":
                    rec_ca.append(atom)
                if name in ("N", "CA", "C", "O"):
                    rec_bb.append(atom)
            elif resname in _NON_LIGAND_HETATMS:
                continue
            else:
                if element.upper() != "H":
                    lig.append(atom)
    return FrameAtoms(receptor_ca=rec_ca, receptor_bb=rec_bb,
                      ligand_heavy=lig, receptor_heavy=rec_heavy)


# ---------------------------------------------------------------------------
# Kabsch alignment (P aligns to Q)
# ---------------------------------------------------------------------------

def kabsch_rotation(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def align_via(p_align: np.ndarray, q_align: np.ndarray,
              p_apply: np.ndarray) -> np.ndarray:
    """Return p_apply transformed to align p_align onto q_align (Kabsch)."""
    p_centroid = p_align.mean(axis=0)
    q_centroid = q_align.mean(axis=0)
    P = p_align - p_centroid
    Q = q_align - q_centroid
    R = kabsch_rotation(P, Q)
    return (p_apply - p_centroid) @ R.T + q_centroid


def rmsd_after_alignment(p_align: np.ndarray, q_align: np.ndarray,
                         p_apply: np.ndarray, q_apply: np.ndarray) -> float:
    pa = align_via(p_align, q_align, p_apply)
    return float(np.sqrt(((pa - q_apply) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------------------
# Core mask discovery from SMILES (taxol-specific via RDKit ring analysis)
# ---------------------------------------------------------------------------

def discover_core_atom_names_from_frame0(
    frame0_path: Path, ligand_smiles: str
) -> Tuple[set, set]:
    """Identify which LIG heavy-atom names belong to the rigid core (largest
    aliphatic fused-ring system) vs the flexible side chain.

    Strategy: parse the LIG heavy-atom block from frame 0, build an RDKit
    Mol via AssignBondOrdersFromTemplate against the docked-ligand SMILES,
    run ring analysis to find the largest aliphatic fused-ring component,
    map those atom indices back to LIG PDB atom names.

    Returns (core_names, sidechain_names) as sets of PDB atom name strings.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # 1) Read the LIG PDB block from frame 0 in HETATM form. Build an RDKit
    # Mol where atom order matches the PDB heavy-atom order.
    lig_block_lines = []
    pdb_atom_names_by_index: List[str] = []
    for line in frame0_path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[17:20].strip() != "LIG":
            continue
        elem = (line[76:78].strip() if len(line) >= 78 else "").upper()
        if elem == "H":
            continue
        lig_block_lines.append("HETATM" + line[6:])
        pdb_atom_names_by_index.append(line[12:16].strip())
    if not lig_block_lines:
        raise RuntimeError(f"no LIG heavy atoms in {frame0_path}")
    lig_block = "\n".join(lig_block_lines) + "\nEND\n"

    rdmol_pdb = Chem.MolFromPDBBlock(lig_block, removeHs=True, sanitize=False)
    if rdmol_pdb is None:
        raise RuntimeError("RDKit could not parse LIG PDB block")
    template = Chem.MolFromSmiles(ligand_smiles)
    if template is None:
        raise RuntimeError(f"invalid SMILES: {ligand_smiles[:60]}…")
    rdmol = AllChem.AssignBondOrdersFromTemplate(template, rdmol_pdb)
    # rdmol's atom order should mirror the PDB heavy-atom order. Sanity check:
    if rdmol.GetNumAtoms() != len(pdb_atom_names_by_index):
        raise RuntimeError(
            f"RDKit atom count {rdmol.GetNumAtoms()} != PDB heavy atom count "
            f"{len(pdb_atom_names_by_index)}"
        )

    # 2) Find largest aliphatic fused-ring system via union-find on ring bonds.
    ri = rdmol.GetRingInfo()
    ring_atoms_all = set()
    for ring in ri.AtomRings():
        # Only count aliphatic rings (the side-chain phenyls are aromatic
        # rings, NOT part of the rigid taxane core).
        is_arom = all(rdmol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_arom:
            ring_atoms_all.update(ring)
    if not ring_atoms_all:
        # No aliphatic rings — fall back to treating ALL ring atoms as core.
        for ring in ri.AtomRings():
            ring_atoms_all.update(ring)
    parent = {a: a for a in ring_atoms_all}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
    for b in rdmol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in ring_atoms_all and j in ring_atoms_all and b.IsInRing():
            union(i, j)
    from collections import defaultdict
    groups = defaultdict(list)
    for a in ring_atoms_all:
        groups[find(a)].append(a)
    largest = max(groups.values(), key=len)
    core_indices = set(largest)

    # 3) Map indices → PDB atom names.
    core_names = {pdb_atom_names_by_index[i] for i in core_indices}
    all_names = set(pdb_atom_names_by_index)
    sidechain_names = all_names - core_names
    return core_names, sidechain_names


# ---------------------------------------------------------------------------
# Per-run analysis pipeline
# ---------------------------------------------------------------------------

@dataclass
class RunSpec:
    label: str               # "Run A", "Run B", "Explicit", "Juliprosopine", "Primaquine"
    md_id: str               # full UUID
    md_dir: Path
    solvent: str             # implicit | explicit
    duration_ps: float       # production length
    snapshot_every_ps: float # ps per snapshot
    construct_note: str      # construct/loop notes for the table
    ligand_label: str        # "lig_000000" etc.
    ligand_smiles: str       # for RDKit core/sidechain mask derivation
    ligand_common_name: str  # "taxol" / "Juliprosopine" / "Primaquine"


def analyze_run(spec: RunSpec, shared_pocket_residues: set) -> dict:
    """Compute the 6 metrics for one run. Returns a dict ready for table emit.

    The core / side-chain mask is derived per-run from the ligand's SMILES
    (largest fused-ring component; aliphatic-preferred, falls back to all
    ring atoms when no aliphatic rings exist — e.g. Primaquine, whose only
    rings are the aromatic quinoline pair).
    """
    # Per-run core mask
    core_names, sidechain_names = discover_core_atom_names_from_frame0(
        spec.md_dir / "frames" / "frame_000.pdb", spec.ligand_smiles,
    )
    frames_dir = spec.md_dir / "frames"
    frame_paths = sorted(frames_dir.glob("frame_*.pdb"))
    n_frames = len(frame_paths)
    print(f"\n[{spec.label}] {spec.md_id[:8]}  frames={n_frames}  "
          f"solvent={spec.solvent}  duration={spec.duration_ps} ps")
    snapshots = [parse_pdb_frame(p) for p in frame_paths]
    times_ps = [i * spec.snapshot_every_ps for i in range(n_frames)]
    eq_start_idx = int(round(n_frames * (1.0 - EQ_WINDOW_LAST_FRAC)))

    # Frame 0 references
    f0 = snapshots[0]
    rec0_ca = f0.receptor_ca
    lig0 = f0.ligand_heavy

    rec0_ca_pos = np.asarray([a["xyz"] for a in rec0_ca])
    lig0_heavy = np.asarray([a["xyz"] for a in lig0])
    lig0_names = [a["name"] for a in lig0]
    lig0_name_to_idx = {n: i for i, n in enumerate(lig0_names)}
    core_mask = np.array([n in core_names for n in lig0_names], dtype=bool)
    side_mask = np.array([n in sidechain_names for n in lig0_names], dtype=bool)
    n_core_present = int(core_mask.sum())
    n_side_present = int(side_mask.sum())
    print(f"  ligand heavy atoms: {len(lig0_names)}  "
          f"core mask: {n_core_present}  side mask: {n_side_present}")

    # Pocket residues for THIS run = intersection of the run's available
    # residues with the shared pocket-residue keys.
    rec_keys = {(a["chain"], a["resseq"]) for a in rec0_ca}
    pocket_keys_present = shared_pocket_residues & rec_keys
    pocket_order = [(a["chain"], a["resseq"]) for a in rec0_ca
                    if (a["chain"], a["resseq"]) in pocket_keys_present]
    pocket0_pos = np.asarray([
        next(a["xyz"] for a in rec0_ca
             if (a["chain"], a["resseq"]) == k)
        for k in pocket_order
    ])
    print(f"  pocket Cα residues: {len(pocket_order)} (shared set: "
          f"{len(shared_pocket_residues)})")

    # Pocket backbone (N,CA,C,O of pocket residues), in frame 0
    # Build pocket-bb atom selector by (chain, resseq, name) for the run's frame 0
    pocket_bb_selector = []
    for a in f0.receptor_bb:
        if (a["chain"], a["resseq"]) in pocket_keys_present:
            pocket_bb_selector.append((a["chain"], a["resseq"], a["name"]))
    pocket_bb0_pos = np.asarray([
        next(a["xyz"] for a in f0.receptor_bb
             if (a["chain"], a["resseq"], a["name"]) == k)
        for k in pocket_bb_selector
    ])

    # Whole-protein backbone
    rec0_bb_pos = np.asarray([a["xyz"] for a in f0.receptor_bb])

    # ----- Initial contact set (metric 4) -----
    # Contacts at frame 0: ligand heavy ↔ protein heavy within 4 Å.
    lig0_arr = np.asarray([a["xyz"] for a in lig0])
    rec0_heavy = np.asarray([a["xyz"] for a in f0.receptor_heavy])
    rec0_heavy_keys = [(a["chain"], a["resseq"], a["resname"], a["name"])
                       for a in f0.receptor_heavy]
    # Compute distance matrix; iterate via numpy for speed.
    dmat = np.sqrt(((lig0_arr[:, None, :] - rec0_heavy[None, :, :]) ** 2).sum(axis=2))
    contact_pairs = np.argwhere(dmat < CONTACT_CUTOFF_A)
    # Each entry: (lig_idx, rec_heavy_idx)
    initial_contacts: List[Tuple[int, Tuple[str, int, str, str]]] = []
    for li, ri in contact_pairs:
        initial_contacts.append((int(li), rec0_heavy_keys[int(ri)]))
    n_initial_contacts = len(initial_contacts)
    print(f"  initial contacts (ligand−protein ≤ {CONTACT_CUTOFF_A:.1f} Å): "
          f"{n_initial_contacts}")

    # ----- Per-frame metrics -----
    pose_global = []
    pose_pocket = []
    core_pocket = []
    side_pocket = []
    pocket_bb_rmsd = []
    wholeprot_bb_rmsd = []
    ca_per_residue_aligned: List[np.ndarray] = []  # for RMSF
    lig_com_displ = []
    in_pocket_flag = []
    buried_frac_series = []
    contacts_retained_frac = []
    # Per-contact persistence count
    contact_persist_counts = [0] * n_initial_contacts

    for idx, snap in enumerate(snapshots):
        rec_ca_pos = np.asarray([a["xyz"] for a in snap.receptor_ca])
        lig_heavy = np.asarray([a["xyz"] for a in snap.ligand_heavy])
        if (rec_ca_pos.shape != rec0_ca_pos.shape
                or lig_heavy.shape != lig0_heavy.shape):
            # Skip mismatched-shape frames (defensive); produce NaNs
            pose_global.append(float("nan"))
            pose_pocket.append(float("nan"))
            core_pocket.append(float("nan"))
            side_pocket.append(float("nan"))
            pocket_bb_rmsd.append(float("nan"))
            wholeprot_bb_rmsd.append(float("nan"))
            lig_com_displ.append(float("nan"))
            in_pocket_flag.append(False)
            buried_frac_series.append(float("nan"))
            contacts_retained_frac.append(float("nan"))
            ca_per_residue_aligned.append(np.full_like(rec0_ca_pos, np.nan))
            continue
        # GLOBAL Cα alignment
        pose_global.append(
            rmsd_after_alignment(rec_ca_pos, rec0_ca_pos, lig_heavy, lig0_heavy)
        )
        # POCKET Cα alignment — reuse one Kabsch for all pocket-aligned metrics
        pocket_pos_i = np.asarray([
            next(a["xyz"] for a in snap.receptor_ca
                 if (a["chain"], a["resseq"]) == k)
            for k in pocket_order
        ])
        # Align using pocket Cα; transform the rest from the same fit
        p_centroid = pocket_pos_i.mean(axis=0)
        q_centroid = pocket0_pos.mean(axis=0)
        P = pocket_pos_i - p_centroid
        Q = pocket0_pos - q_centroid
        R = kabsch_rotation(P, Q)
        def xfm(arr_xyz):
            return (arr_xyz - p_centroid) @ R.T + q_centroid

        lig_xfm = xfm(lig_heavy)
        pose_pocket.append(
            float(np.sqrt(((lig_xfm - lig0_heavy) ** 2).sum(axis=1).mean()))
        )
        if n_core_present > 0:
            core_pocket.append(
                float(np.sqrt(((lig_xfm[core_mask] - lig0_heavy[core_mask]) ** 2)
                              .sum(axis=1).mean()))
            )
        else:
            core_pocket.append(float("nan"))
        if n_side_present > 0:
            side_pocket.append(
                float(np.sqrt(((lig_xfm[side_mask] - lig0_heavy[side_mask]) ** 2)
                              .sum(axis=1).mean()))
            )
        else:
            side_pocket.append(float("nan"))
        # Pocket BB RMSD — align pocket Cα then RMSD on pocket BB (N,CA,C,O)
        bb_pos_i = np.asarray([
            next(a["xyz"] for a in snap.receptor_bb
                 if (a["chain"], a["resseq"], a["name"]) == k)
            for k in pocket_bb_selector
        ])
        bb_pos_i_xfm = xfm(bb_pos_i)
        pocket_bb_rmsd.append(
            float(np.sqrt(((bb_pos_i_xfm - pocket_bb0_pos) ** 2)
                          .sum(axis=1).mean()))
        )
        # Whole-protein BB — global Cα align THEN BB RMSD
        rec_bb_pos = np.asarray([a["xyz"] for a in snap.receptor_bb])
        wholeprot_bb_rmsd.append(
            rmsd_after_alignment(rec_ca_pos, rec0_ca_pos, rec_bb_pos, rec0_bb_pos)
        )
        # Per-residue Cα aligned (for RMSF over equilibrated window)
        rec_ca_xfm_global = align_via(rec_ca_pos, rec0_ca_pos, rec_ca_pos)
        ca_per_residue_aligned.append(rec_ca_xfm_global)
        # Ligand center of mass displacement (POCKET-aligned, since we want
        # binding-mode drift, not global translation/rotation)
        com_i = lig_xfm.mean(axis=0)
        com_0 = lig0_heavy.mean(axis=0)
        com_d = float(np.linalg.norm(com_i - com_0))
        lig_com_displ.append(com_d)
        in_pocket_flag.append(com_d < IN_POCKET_COM_A)
        # Buried fraction: per ligand heavy atom, any protein heavy atom within 4.5 Å
        rec_heavy = np.asarray([a["xyz"] for a in snap.receptor_heavy])
        dmat_i = np.sqrt(((lig_heavy[:, None, :] - rec_heavy[None, :, :]) ** 2)
                         .sum(axis=2))
        buried = (dmat_i.min(axis=1) < BURIED_CUTOFF_A).mean()
        buried_frac_series.append(float(buried))
        # Contacts: re-check each initial pair against the run's CURRENT
        # protein heavy atoms (selected by (chain, resseq, resname, name)).
        # Build a quick lookup once per frame for speed.
        rec_lookup = {(a["chain"], a["resseq"], a["resname"], a["name"]): a["xyz"]
                      for a in snap.receptor_heavy}
        retained = 0
        for ci, (li, key) in enumerate(initial_contacts):
            ri_xyz = rec_lookup.get(key)
            if ri_xyz is None:
                continue
            l_xyz = lig_heavy[li]
            d = ((ri_xyz[0] - l_xyz[0]) ** 2 + (ri_xyz[1] - l_xyz[1]) ** 2
                 + (ri_xyz[2] - l_xyz[2]) ** 2) ** 0.5
            if d < CONTACT_CUTOFF_A:
                retained += 1
                contact_persist_counts[ci] += 1
        contacts_retained_frac.append(
            retained / n_initial_contacts if n_initial_contacts > 0 else float("nan")
        )

    # --- RMSF over equilibrated window (per-residue Cα fluctuation) ---
    ca_arr = np.array(ca_per_residue_aligned[eq_start_idx:])  # (T, N, 3)
    if ca_arr.size > 0:
        ca_mean = np.nanmean(ca_arr, axis=0)
        diffs = ca_arr - ca_mean[None, :, :]
        sqd = (diffs ** 2).sum(axis=2)
        rmsf_per_residue = np.sqrt(np.nanmean(sqd, axis=0))
        # Restrict to pocket Cα RMSF
        pocket_resi_idx = [i for i, a in enumerate(rec0_ca)
                           if (a["chain"], a["resseq"]) in pocket_keys_present]
        pocket_rmsf_mean = float(np.nanmean(rmsf_per_residue[pocket_resi_idx]))
    else:
        pocket_rmsf_mean = float("nan")

    # Helpers
    def nm_mean(arr, lo=None, hi=None):
        a = np.asarray(arr)
        if lo is None and hi is None:
            return float(np.nanmean(a))
        return float(np.nanmean(a[lo:hi]))
    def nm_final(arr):
        return float(arr[-1])
    def nm_max(arr):
        return float(np.nanmax(np.asarray(arr)))

    # Block averages for convergence (metric 6)
    half = n_frames // 2
    pose_pocket_arr = np.asarray(pose_pocket)
    first_half_mean = nm_mean(pose_pocket, 0, half)
    last_half_mean = nm_mean(pose_pocket, half, n_frames)
    # Slope of pose_pocket over the last half via least-squares fit
    if half < n_frames - 1:
        x = np.asarray(times_ps[half:])
        y = pose_pocket_arr[half:]
        # Filter NaN
        mask = ~np.isnan(y)
        if mask.sum() >= 2:
            slope, _ = np.polyfit(x[mask], y[mask], 1)
        else:
            slope = float("nan")
    else:
        slope = float("nan")

    summary = {
        "label": spec.label,
        "md_id": spec.md_id,
        "solvent": spec.solvent,
        "duration_ps": spec.duration_ps,
        "snapshot_every_ps": spec.snapshot_every_ps,
        "construct_note": spec.construct_note,
        "ligand_label": spec.ligand_label,
        "ligand_common_name": spec.ligand_common_name,
        "ligand_smiles": spec.ligand_smiles,
        "n_frames": n_frames,
        "n_eq_window_frames": n_frames - eq_start_idx,
        "n_pocket_residues_shared_intersection": len(pocket_order),
        "n_core_atoms": n_core_present,
        "n_sidechain_atoms": n_side_present,

        # Metric 1
        "pose_global_final": nm_final(pose_global),
        "pose_global_max":   nm_max(pose_global),
        "pose_global_mean_eq": nm_mean(pose_global, eq_start_idx, None),
        "pose_pocket_final": nm_final(pose_pocket),
        "pose_pocket_max":   nm_max(pose_pocket),
        "pose_pocket_mean_eq": nm_mean(pose_pocket, eq_start_idx, None),

        # Metric 2
        "core_pocket_final": nm_final(core_pocket),
        "core_pocket_max":   nm_max(core_pocket),
        "core_pocket_mean_eq": nm_mean(core_pocket, eq_start_idx, None),
        "sidechain_pocket_final": nm_final(side_pocket),
        "sidechain_pocket_max":   nm_max(side_pocket),
        "sidechain_pocket_mean_eq": nm_mean(side_pocket, eq_start_idx, None),

        # Metric 3
        "pocket_bb_final":  nm_final(pocket_bb_rmsd),
        "pocket_bb_max":    nm_max(pocket_bb_rmsd),
        "pocket_bb_mean_eq": nm_mean(pocket_bb_rmsd, eq_start_idx, None),
        "pocket_rmsf_mean_eq": pocket_rmsf_mean,
        "wholeprot_bb_final": nm_final(wholeprot_bb_rmsd),
        "wholeprot_bb_max":   nm_max(wholeprot_bb_rmsd),

        # Metric 4
        "n_initial_contacts": n_initial_contacts,
        "contacts_retained_frac_mean_eq":
            nm_mean(contacts_retained_frac, eq_start_idx, None),
        "per_contact_occupancy": [
            {
                "lig_atom_name": lig0_names[li],
                "protein_residue": f"{key[0]}{key[1]}{key[2]}",
                "protein_atom_name": key[3],
                "occupancy": (cnt / n_frames) if n_frames > 0 else 0.0,
            }
            for (li, key), cnt in zip(initial_contacts, contact_persist_counts)
        ],

        # Metric 5
        "lig_COM_displ_final": nm_final(lig_com_displ),
        "lig_COM_displ_max":   nm_max(lig_com_displ),
        "in_pocket_frac_eq": float(
            np.mean(in_pocket_flag[eq_start_idx:])
        ) if eq_start_idx < n_frames else float("nan"),
        "buried_frac_mean_eq": nm_mean(buried_frac_series, eq_start_idx, None),

        # Metric 6
        "conv_first_half_mean_pose_pocket": first_half_mean,
        "conv_last_half_mean_pose_pocket": last_half_mean,
        "conv_last_half_slope_pose_pocket_a_per_ps": slope,

        # Time series for figures (clipped)
        "_series_times_ps": times_ps,
        "_series_pose_global": pose_global,
        "_series_pose_pocket": pose_pocket,
        "_series_core_pocket": core_pocket,
        "_series_sidechain_pocket": side_pocket,
        "_series_pocket_bb": pocket_bb_rmsd,
        "_series_lig_com_displ": lig_com_displ,
    }
    return summary


# ---------------------------------------------------------------------------
# Pocket-residue shared set derivation
# ---------------------------------------------------------------------------

def derive_shared_pocket_residues(specs: List[RunSpec]) -> set:
    """Pocket = union over all runs of residues with ANY heavy atom within
    POCKET_RADIUS_A of any ligand heavy atom in frame 0. Using the union
    (not intersection) preserves the binding-site definition even when one
    run truncates a flanking residue. Each run's per-frame computation
    then restricts to the residues IT actually has.
    """
    all_keys = set()
    for spec in specs:
        frame0 = spec.md_dir / "frames" / "frame_000.pdb"
        f0 = parse_pdb_frame(frame0)
        lig = np.asarray([a["xyz"] for a in f0.ligand_heavy])
        rec_heavy_keys = []
        rec_heavy_pos = []
        for a in f0.receptor_heavy:
            rec_heavy_keys.append((a["chain"], a["resseq"]))
            rec_heavy_pos.append(a["xyz"])
        rec_heavy_pos = np.asarray(rec_heavy_pos)
        # any heavy atom of residue within 5 Å of any ligand heavy atom
        dmat = np.sqrt(((lig[:, None, :] - rec_heavy_pos[None, :, :]) ** 2)
                       .sum(axis=2))
        close = (dmat.min(axis=0) < POCKET_RADIUS_A)
        for i, is_close in enumerate(close):
            if is_close:
                all_keys.add(rec_heavy_keys[i])
    return all_keys


# ---------------------------------------------------------------------------
# Figure rendering
# ---------------------------------------------------------------------------

def plot_overlay(series_per_run: List[Tuple[str, List[float], List[float]]],
                 ylabel: str, title: str, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    for label, times_ps, values in series_per_run:
        ax.plot(times_ps, values, "-", label=label, linewidth=1.2)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"  Wrote {out_path}")


# ---------------------------------------------------------------------------
# Cross-run consensus RMSD (pose-selection robustness)
# ---------------------------------------------------------------------------

def compute_pairwise_consensus_rmsds(
    specs: List[RunSpec], shared_pocket_residues: set,
) -> Dict:
    """Pairwise RMSD between equilibrated-window MEAN ligand poses of the
    given runs, computed in a single common reference frame.

    Procedure:
      1. Use the FIRST spec's frame 0 as the common reference. Capture its
         pocket-Cα atoms and ligand heavy atoms.
      2. For each run, walk the equilibrated window; pocket-align every
         frame to the common reference; emit the transformed ligand
         heavy-atom positions.
      3. Take the per-atom MEAN of those positions across the window — the
         "consensus" ligand pose for that run, in the common frame.
      4. Pairwise: for each (i, j) pair, RMSD between consensus pose i and
         consensus pose j (no further alignment — they live in the same
         frame already).

    The result tests "do the alternative-pose MDs relax toward the same
    bound mode as pose 0?". Small pairwise RMSDs → consensus (poses
    converged to one bound mode). Large pairwise RMSDs → distinct
    persistent modes.

    Requires all input runs to share the same ligand SMILES (same atom
    set and ordering) — otherwise the per-atom mean is meaningless.
    """
    n = len(specs)
    print(f"\nComputing pairwise consensus RMSDs across {n} runs ...")

    # --- common reference: first run's frame 0
    ref_spec = specs[0]
    ref_frame = parse_pdb_frame(ref_spec.md_dir / "frames" / "frame_000.pdb")
    rec0_ca = ref_frame.receptor_ca
    rec_keys = {(a["chain"], a["resseq"]) for a in rec0_ca}
    pocket_keys = shared_pocket_residues & rec_keys
    pocket_order = [(a["chain"], a["resseq"]) for a in rec0_ca
                    if (a["chain"], a["resseq"]) in pocket_keys]
    pocket0_pos = np.asarray([
        next(a["xyz"] for a in rec0_ca
             if (a["chain"], a["resseq"]) == k)
        for k in pocket_order
    ])
    ref_lig_names = [a["name"] for a in ref_frame.ligand_heavy]
    n_lig = len(ref_lig_names)
    print(f"  common reference = {ref_spec.label} frame 0; "
          f"{len(pocket_order)} pocket Cα, {n_lig} ligand heavy atoms")

    consensus_poses: List[np.ndarray] = []
    labels: List[str] = []
    for spec in specs:
        frames_dir = spec.md_dir / "frames"
        frame_paths = sorted(frames_dir.glob("frame_*.pdb"))
        n_frames = len(frame_paths)
        eq_start = int(round(n_frames * (1.0 - EQ_WINDOW_LAST_FRAC)))
        print(f"  {spec.label}: n_frames={n_frames}, "
              f"equilibrated window = frames {eq_start}..{n_frames-1}")

        # Verify ligand atom ordering matches the reference
        snap0 = parse_pdb_frame(frame_paths[0])
        names_this = [a["name"] for a in snap0.ligand_heavy]
        if names_this != ref_lig_names:
            raise RuntimeError(
                f"{spec.label} ligand atom names diverge from reference — "
                f"cannot compute consensus across heterogeneous topologies"
            )

        # Walk the equilibrated window; transform ligand atoms into the
        # common reference frame via pocket-Cα alignment.
        lig_transformed_per_frame = []
        for fp in frame_paths[eq_start:]:
            snap = parse_pdb_frame(fp)
            pocket_pos_i = np.asarray([
                next(a["xyz"] for a in snap.receptor_ca
                     if (a["chain"], a["resseq"]) == k)
                for k in pocket_order
            ])
            lig_pos_i = np.asarray([a["xyz"] for a in snap.ligand_heavy])
            lig_xfm = align_via(pocket_pos_i, pocket0_pos, lig_pos_i)
            lig_transformed_per_frame.append(lig_xfm)
        # MEAN consensus pose for this run in the common frame
        consensus = np.mean(np.array(lig_transformed_per_frame), axis=0)
        consensus_poses.append(consensus)
        labels.append(spec.label)

    # Pairwise RMSD matrix
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(0.0)
            else:
                rmsd = float(np.sqrt(((consensus_poses[i] - consensus_poses[j])
                                      ** 2).sum(axis=1).mean()))
                row.append(rmsd)
        matrix.append(row)

    return {
        "labels": labels,
        "matrix": matrix,
        "reference": ref_spec.label,
        "n_pocket_residues": len(pocket_order),
        "n_ligand_heavy": n_lig,
        "equilibrated_window_last_fraction": EQ_WINDOW_LAST_FRAC,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

JOBS_ROOT = Path(__file__).resolve().parents[1] / "jobs"
DOCS_ROOT = Path(__file__).resolve().parents[1] / "docs"
JOBS_ROOT.mkdir(parents=True, exist_ok=True)
DOCS_ROOT.mkdir(parents=True, exist_ok=True)
TAXOL_SMILES = ("O[C@@]12[C@@H](OC(=O)c3ccccc3)[C@@H]3[C@]4(OC(=O)C)CO[C@@H]4"
                "C[C@H](O)[C@@]3(C)C(=O)[C@H](OC(=O)C)C(=C([C@@H](OC(=O)[C@H]"
                "(O)[C@@H](NC(=O)c3ccccc3)c3ccccc3)C1)C)C2(C)C")
JULIPROSOPINE_SMILES = ("C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@]"
                        "(O3)(OO4)C)C")
PRIMAQUINE_SMILES = "COc1cc2cccnc2c(NC(C)CCCN)c1"

# Same explicit-mode protocol as the [B9] Explicit taxol run: truncated
# octahedron + α E-hook truncated + OXT cap + Zn stripped + 1 ns NPT.
EXPLICIT_NOTE = ("explicit TIP3P, NPT, truncated octahedron (1.0 nm padding, "
                 "0.15 M); α C-terminal E-hook 439–451 truncated (OXT cap @ "
                 "A438); Mg retained, Zn stripped (artifact)")

RUNS: List[RunSpec] = [
    # Taxol — bona-fide binder, three runs at different protocols.
    RunSpec(
        label="Run A",
        md_id="1f01da83-ac18-4103-8631-2053873967c0",
        md_dir=JOBS_ROOT / "1f01da83-ac18-4103-8631-2053873967c0" / "md",
        solvent="implicit",
        duration_ps=500.0, snapshot_every_ps=5.0,
        construct_note=("taxol implicit OBC2; full receptor (chain A 1–451 + "
                       "chain B 2–437); Mg + Zn both retained"),
        ligand_label="lig_000000",
        ligand_smiles=TAXOL_SMILES,
        ligand_common_name="taxol",
    ),
    RunSpec(
        label="Run B",
        md_id="5817afff-8e6b-4130-9997-46c9965379b2",
        md_dir=JOBS_ROOT / "5817afff-8e6b-4130-9997-46c9965379b2" / "md",
        solvent="implicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=("taxol implicit OBC2; full receptor (chain A 1–451 + "
                       "chain B 2–437); Mg + Zn both retained"),
        ligand_label="lig_000000",
        ligand_smiles=TAXOL_SMILES,
        ligand_common_name="taxol",
    ),
    RunSpec(
        label="Explicit (pose 0)",
        md_id="5bc61f59-834f-4e71-a492-d32ddfdc7326",
        md_dir=JOBS_ROOT / "5bc61f59-834f-4e71-a492-d32ddfdc7326" / "md",
        solvent="explicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=f"taxol pose_rank 0 (vina −8.68, cluster 1); {EXPLICIT_NOTE}",
        ligand_label="lig_000000",
        ligand_smiles=TAXOL_SMILES,
        ligand_common_name="taxol",
    ),
    # [B11] pose-selection robustness: same construct/protocol as pose 0,
    # only the starting pose differs. pose 1 and pose 2 sit in DIFFERENT
    # docking clusters from pose 0 (rmsd_to_best 6.05 / 5.99 Å respectively)
    # so this is a real alternative-pose test, not a refinement.
    RunSpec(
        label="Explicit (pose 1)",
        md_id="34840aa1-cbf5-4c6b-a665-ea4f52110f5d",
        md_dir=JOBS_ROOT / "34840aa1-cbf5-4c6b-a665-ea4f52110f5d" / "md",
        solvent="explicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=f"taxol pose_rank 1 (vina −8.48, cluster 0); {EXPLICIT_NOTE}",
        ligand_label="lig_000000",
        ligand_smiles=TAXOL_SMILES,
        ligand_common_name="taxol",
    ),
    RunSpec(
        label="Explicit (pose 2)",
        md_id="80e53d8a-1926-4525-b2e1-55cb1e30eedd",
        md_dir=JOBS_ROOT / "80e53d8a-1926-4525-b2e1-55cb1e30eedd" / "md",
        solvent="explicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=f"taxol pose_rank 2 (vina −8.34, cluster 0); {EXPLICIT_NOTE}",
        ligand_label="lig_000000",
        ligand_smiles=TAXOL_SMILES,
        ligand_common_name="taxol",
    ),
    # Controls — repurposed compounds docked into the taxane box, NOT
    # known binders. Identical explicit-mode protocol so any deviation in
    # behaviour reflects ligand identity, not protocol drift.
    RunSpec(
        label="Juliprosopine",
        md_id="a0b04941-e4b0-40ef-9459-67becac4a61c",
        md_dir=JOBS_ROOT / "a0b04941-e4b0-40ef-9459-67becac4a61c" / "md",
        solvent="explicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=f"Juliprosopine (repurposed; non-binder candidate); {EXPLICIT_NOTE}",
        ligand_label="lig_000001",
        ligand_smiles=JULIPROSOPINE_SMILES,
        ligand_common_name="Juliprosopine",
    ),
    RunSpec(
        label="Primaquine",
        md_id="00dda37f-e218-49ce-8ac0-a69bfff6851d",
        md_dir=JOBS_ROOT / "00dda37f-e218-49ce-8ac0-a69bfff6851d" / "md",
        solvent="explicit",
        duration_ps=1000.0, snapshot_every_ps=5.0,
        construct_note=f"Primaquine (repurposed; non-binder candidate); {EXPLICIT_NOTE}",
        ligand_label="lig_000002",
        ligand_smiles=PRIMAQUINE_SMILES,
        ligand_common_name="Primaquine",
    ),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(DOCS_ROOT),
                    help="Directory to write CSV/MD/PNG output (default: docs/)")
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== Stability diagnostics ===")
    print(f"Output dir: {out_dir}")
    print(f"Runs: {[r.label for r in RUNS]}")

    # 1) Derive shared pocket residue set (union across all runs' frame 0).
    # The core/sidechain mask is derived PER RUN from its ligand's SMILES,
    # so different ligands get different masks (taxol's 17-atom baccatin
    # vs Primaquine's ~10-atom quinoline, etc.).
    print(f"\nDeriving pocket residue set ({POCKET_RADIUS_A:.1f} Å of any "
          f"ligand atom, union across runs) ...")
    shared_pocket = derive_shared_pocket_residues(RUNS)
    print(f"  Pocket residues (union): {len(shared_pocket)}")
    print(f"  Example: {sorted(shared_pocket)[:10]}")

    # 2) Run analysis per run
    results = []
    for spec in RUNS:
        result = analyze_run(spec, shared_pocket)
        results.append(result)
        # JSON dump
        json_path = out_dir / f"B9_stability_{spec.label.replace(' ', '_').lower()}.json"
        # Strip the long series for the per-run JSON (still keep them for figures)
        json_payload = {k: v for k, v in result.items() if not k.startswith("_series_")}
        json_payload["per_contact_occupancy"] = result["per_contact_occupancy"]
        json_path.write_text(json.dumps(json_payload, indent=2))
        print(f"  Wrote {json_path}")

    # 4) Figures (overlays)
    print("\nWriting overlay figures ...")
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_pose_pocket"])
         for r in results],
        ylabel="ligand pose RMSD (Å)  (pocket-aligned)",
        title="Whole-ligand pose RMSD — pocket-aligned",
        out_path=out_dir / "B9_stability_pose_pocket_overlay.png",
    )
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_pose_global"])
         for r in results],
        ylabel="ligand pose RMSD (Å)  (global Cα align)",
        title="Whole-ligand pose RMSD — global align",
        out_path=out_dir / "B9_stability_pose_global_overlay.png",
    )
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_core_pocket"])
         for r in results],
        ylabel="core RMSD (Å)  (pocket-aligned)",
        title=("Core (baccatin/taxane scaffold) RMSD — pocket-aligned"),
        out_path=out_dir / "B9_stability_core_overlay.png",
    )
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_sidechain_pocket"])
         for r in results],
        ylabel="side-chain RMSD (Å)  (pocket-aligned)",
        title="Side-chain (C13 phenylisoserine etc.) RMSD — pocket-aligned",
        out_path=out_dir / "B9_stability_sidechain_overlay.png",
    )
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_pocket_bb"])
         for r in results],
        ylabel="pocket backbone RMSD (Å)  (pocket Cα align)",
        title="Pocket-local backbone RMSD",
        out_path=out_dir / "B9_stability_pocket_bb_overlay.png",
    )
    plot_overlay(
        [(r["label"], r["_series_times_ps"], r["_series_lig_com_displ"])
         for r in results],
        ylabel="ligand COM displacement (Å) (pocket-aligned)",
        title="Ligand center-of-mass drift from frame 0",
        out_path=out_dir / "B9_stability_lig_com_displ_overlay.png",
    )

    # 5) Write the table — CSV + markdown
    print("\nWriting tables ...")
    cols = [
        "label", "md_id", "solvent", "duration_ps", "n_frames",
        "construct_note",
        "pose_global_final", "pose_global_max", "pose_global_mean_eq",
        "pose_pocket_final", "pose_pocket_max", "pose_pocket_mean_eq",
        "core_pocket_final", "core_pocket_max", "core_pocket_mean_eq",
        "sidechain_pocket_final", "sidechain_pocket_max", "sidechain_pocket_mean_eq",
        "pocket_bb_final", "pocket_bb_max", "pocket_bb_mean_eq",
        "pocket_rmsf_mean_eq",
        "wholeprot_bb_final", "wholeprot_bb_max",
        "n_initial_contacts", "contacts_retained_frac_mean_eq",
        "lig_COM_displ_final", "lig_COM_displ_max",
        "in_pocket_frac_eq", "buried_frac_mean_eq",
        "conv_first_half_mean_pose_pocket",
        "conv_last_half_mean_pose_pocket",
        "conv_last_half_slope_pose_pocket_a_per_ps",
    ]
    csv_path = out_dir / "B9_stability_diagnostics.csv"
    import csv as _csv
    with csv_path.open("w", newline="") as f:
        w = _csv.writer(f)
        w.writerow(cols)
        for r in results:
            w.writerow([r.get(c, "") for c in cols])
    print(f"  Wrote {csv_path}")

    # Markdown table — fewer columns, more readable
    md_cols = [
        ("Run",                "label"),
        ("solvent",            "solvent"),
        ("ps",                 "duration_ps"),
        ("n_frames",           "n_frames"),
        ("pose-glob fin/max",  None),
        ("pose-pkt fin/max",   None),
        ("core-pkt fin/max",   None),
        ("side-pkt fin/max",   None),
        ("pkt-BB fin/max",     None),
        ("pkt-RMSF mean",      "pocket_rmsf_mean_eq"),
        ("whole-BB fin/max",   None),
        ("n_init_contacts",    "n_initial_contacts"),
        ("retained_frac (eq)", "contacts_retained_frac_mean_eq"),
        ("COM displ fin/max",  None),
        ("in_pocket (eq)",     "in_pocket_frac_eq"),
        ("buried (eq)",        "buried_frac_mean_eq"),
        ("1st½ / 2nd½ pose_pkt", None),
        ("2nd½ slope (Å/ps)",  "conv_last_half_slope_pose_pocket_a_per_ps"),
    ]

    def fmt(v):
        if v is None:
            return "—"
        if isinstance(v, float):
            if abs(v) < 1e-4:
                return f"{v:.2e}"
            return f"{v:.2f}"
        return str(v)

    def fmt_pair(r, k_a, k_b):
        return f"{fmt(r[k_a])} / {fmt(r[k_b])}"

    md_path = out_dir / "B9_stability_diagnostics.md"
    n_taxol = sum(1 for r in results if "taxol" in r["construct_note"].lower())
    n_control = len(results) - n_taxol
    lines = [
        "# B9 stability diagnostics — measure-only table (READ-ONLY)",
        "",
        f"Read-only diagnostic across {len(results)} long-form OpenMM runs "
        f"at the verified taxane pocket of docking "
        f"4a37bb0c-7655-411e-a238-1ff5b7bab910 — {n_taxol} taxol "
        f"(lig_000000, pose 0; bona-fide binder) plus {n_control} "
        f"control compounds docked into the same pocket. **No thresholds "
        f"are imposed; no pass/fail is assigned.** The table is the raw "
        f"material for later criterion design.",
        "",
        "**Controls** (lig_000001 = Juliprosopine, lig_000002 = Primaquine): "
        "repurposed-compound docking poses, run through the SAME explicit-"
        "mode protocol as the taxol Explicit run so any deviation in "
        "behaviour reflects ligand identity, not protocol drift.",
        "",
        "**Equilibrated window**: last 30% of frames "
        "(matches `DEFAULT_LAST_FRACTION` in the free-energy gating module).",
        "",
        "**Pocket definition**: protein residues with any heavy atom within "
        f"{POCKET_RADIUS_A:.1f} Å of any ligand heavy atom in frame 0, "
        f"unioned across all {len(results)} runs.",
        "",
        "**Core / side-chain mask** is derived PER RUN from the ligand's "
        "SMILES via RDKit ring-system analysis (largest aliphatic fused-ring "
        "component when present; falls back to the largest ring system overall "
        "for ligands like Primaquine whose only rings are aromatic). "
        "Per-run atom counts (n_core_atoms / n_sidechain_atoms heavy):",
        "",
    ]
    for r in results:
        # Pick a short "what is this run" tag for the per-run mask list —
        # ligand_common_name when available, otherwise the first phrase of
        # the construct note (strip after the first comma to avoid the
        # long full-protocol blurb).
        short = r.get("ligand_common_name") or (
            r.get("construct_note", "").split(",")[0].rstrip(";")
        )
        lines.append(
            f"  - **{r['label']}** ({short}): "
            f"core={r['n_core_atoms']}, side-chain={r['n_sidechain_atoms']} "
            f"(of {r['n_core_atoms'] + r['n_sidechain_atoms']} ligand heavy atoms)"
        )
    lines += [
        "",
        "## Per-run differences (do not normalize away)",
        "",
        *[f"- **{r['label']}** ({r['md_id']}): {r['construct_note']}"
          for r in results],
        "",
        "## Table",
        "",
    ]

    header = "| " + " | ".join(c for c, _ in md_cols) + " |"
    sep = "|" + "|".join(["---"] * len(md_cols)) + "|"
    lines.append(header)
    lines.append(sep)
    for r in results:
        row = []
        for label_md, k in md_cols:
            if k is not None:
                row.append(fmt(r.get(k)))
            elif label_md.startswith("pose-glob"):
                row.append(fmt_pair(r, "pose_global_final", "pose_global_max"))
            elif label_md.startswith("pose-pkt"):
                row.append(fmt_pair(r, "pose_pocket_final", "pose_pocket_max"))
            elif label_md.startswith("core-pkt"):
                row.append(fmt_pair(r, "core_pocket_final", "core_pocket_max"))
            elif label_md.startswith("side-pkt"):
                row.append(fmt_pair(r, "sidechain_pocket_final", "sidechain_pocket_max"))
            elif label_md.startswith("pkt-BB"):
                row.append(fmt_pair(r, "pocket_bb_final", "pocket_bb_max"))
            elif label_md.startswith("whole-BB"):
                row.append(fmt_pair(r, "wholeprot_bb_final", "wholeprot_bb_max"))
            elif label_md.startswith("COM displ"):
                row.append(fmt_pair(r, "lig_COM_displ_final", "lig_COM_displ_max"))
            elif label_md.startswith("1st½"):
                row.append(fmt_pair(r, "conv_first_half_mean_pose_pocket",
                                       "conv_last_half_mean_pose_pocket"))
            else:
                row.append("—")
        lines.append("| " + " | ".join(row) + " |")

    lines += [
        "",
        "Units: lengths in Å; durations in ps. Slope is least-squares fit "
        "of pose-pocket RMSD over the second half. `retained_frac (eq)` = "
        "mean over the equilibrated window of the fraction of frame-0 "
        "contacts retained per frame.",
        "",
        "## Notes per metric",
        "",
        "1. **Pose RMSD** — `pose-glob` aligns on whole-protein Cα; `pose-pkt` "
        f"aligns on the {POCKET_RADIUS_A:.1f}-Å pocket Cα set. The same pocket "
        "alignment frame is reused for `core-pkt`, `side-pkt`, and `pkt-BB` "
        "(see the script: one Kabsch per frame, four downstream RMSDs).",
        "2. **Core vs side-chain** — `core` is the largest rigid fused-ring "
        "component of the ligand (per-ligand mask derived from SMILES via "
        "RDKit ring-system analysis); `side-chain` is everything else. "
        "For taxol the core is the baccatin/oxetane scaffold (17 heavy "
        "atoms) with the C13 phenylisoserine + the C2 OBz / C10 / C4 OAc "
        "decorations as side-chain (45 atoms). For Juliprosopine the core "
        "is the polycyclic peroxide scaffold (16 atoms, the whole molecule "
        "minus 4 methyls). For Primaquine the core is the quinoline ring "
        "(10 atoms) with the flexible aminopentyl chain as side-chain "
        "(9 atoms).",
        "3. **Pocket BB RMSD + RMSF** — pocket BB is N/CA/C/O of the pocket "
        "residue set; RMSF is per-residue Cα fluctuation in the same "
        "global-Cα-aligned frame, reduced over the equilibrated window.",
        "4. **Contacts** — every (ligand heavy atom × protein heavy atom) "
        f"pair within {CONTACT_CUTOFF_A:.1f} Å of each other in frame 0 is an "
        "initial contact. Per-frame retention is checked against the same "
        f"{CONTACT_CUTOFF_A:.1f}-Å threshold; the per-contact occupancy table "
        "lives in the per-run JSON dumps.",
        "5. **Pocket occupancy** — `COM displ` is computed AFTER pocket-Cα "
        f"alignment (binding-mode drift, not bulk translation); `in_pocket "
        f"(eq)` = fraction of equilibrated-window frames with COM displ < "
        f"{IN_POCKET_COM_A:.1f} Å; `buried (eq)` = mean fraction of ligand "
        f"heavy atoms with any protein heavy atom within {BURIED_CUTOFF_A:.1f} "
        "Å, over the equilibrated window.",
        "6. **Convergence** — `1st½ / 2nd½` are block-averaged "
        "pose-pocket RMSD; `2nd½ slope` is the linear-fit slope over the "
        "second half. Plateaued → similar block averages + slope near 0; "
        "still rising → 2nd½ noticeably higher than 1st½ + positive slope.",
        "",
        "## Figures",
        "",
        f"All saved as overlay PNGs ({len(results)} runs, color-coded):",
        "",
        "- `B9_stability_pose_pocket_overlay.png` — whole-ligand pose RMSD (pocket-aligned)",
        "- `B9_stability_pose_global_overlay.png` — whole-ligand pose RMSD (global align)",
        "- `B9_stability_core_overlay.png` — baccatin-core RMSD (pocket-aligned)",
        "- `B9_stability_sidechain_overlay.png` — side-chain RMSD (pocket-aligned)",
        "- `B9_stability_pocket_bb_overlay.png` — pocket backbone RMSD",
        "- `B9_stability_lig_com_displ_overlay.png` — ligand COM displacement",
        "",
    ]
    md_path.write_text("\n".join(lines))
    print(f"  Wrote {md_path}")

    # 6) Pose-selection consensus across the taxol runs (B11).
    # Only run when there are >= 2 taxol runs; computing it requires identical
    # ligand atom orderings, which is only guaranteed within a single
    # ligand SMILES.
    taxol_specs = [r for r in RUNS if r.ligand_common_name == "taxol"]
    if len(taxol_specs) >= 2:
        consensus = compute_pairwise_consensus_rmsds(taxol_specs, shared_pocket)
        consensus_path = out_dir / "B11_consensus_rmsd.json"
        consensus_path.write_text(json.dumps(consensus, indent=2))
        print(f"  Wrote {consensus_path}")
        # Append a summary section to the existing markdown
        lines = md_path.read_text().splitlines()
        lines += [
            "",
            "## Pose-selection consensus across taxol runs (B11)",
            "",
            f"Pairwise RMSD (Å) between equilibrated-window MEAN ligand "
            f"poses, in the common reference frame "
            f"({consensus['reference']} frame 0, pocket-Cα aligned). "
            f"Equilibrated window = last "
            f"{int(consensus['equilibrated_window_last_fraction']*100)} % of "
            f"frames; {consensus['n_pocket_residues']} pocket Cα atoms drive "
            f"the alignment; {consensus['n_ligand_heavy']} heavy atoms in "
            "each consensus pose.",
            "",
        ]
        header = "| | " + " | ".join(consensus["labels"]) + " |"
        sep = "|" + "|".join(["---"] * (len(consensus["labels"]) + 1)) + "|"
        lines.append(header)
        lines.append(sep)
        for i, lbl in enumerate(consensus["labels"]):
            row = ["**" + lbl + "**"]
            for j in range(len(consensus["labels"])):
                v = consensus["matrix"][i][j]
                row.append(f"{v:.2f}")
            lines.append("| " + " | ".join(row) + " |")
        lines += [
            "",
            "Interpretation:",
            "  - Small off-diagonals (≲ 2 Å) → the alternative-pose MDs "
            "have relaxed toward the same bound mode as the reference "
            "(consensus achieved; pose selection is robust to docking-rank "
            "uncertainty).",
            "  - Large off-diagonals (close to the starting RMSD between "
            "docking poses — ~6 Å for pose 0 vs poses 1/2 here) → the "
            "alternative poses stayed in their own basins (pose selection "
            "is NOT robust at this timescale).",
            "",
        ]
        md_path.write_text("\n".join(lines))
        print(f"  Appended consensus block to {md_path.name}")

    print("\nDone.")


if __name__ == "__main__":
    main()
