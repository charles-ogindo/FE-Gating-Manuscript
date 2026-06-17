"""
MD analysis: takes a list of snapshot PDB files (one per saved frame), each
containing receptor + ligand atoms in the same coordinate frame, and emits:

  - rmsd.csv     : backbone-Cα RMSD and ligand heavy-atom RMSD vs the t=0 frame
  - hbonds.csv   : per-frame protein-ligand H-bond count + persistence summary
  - contacts.csv : per-residue contact frequency (fraction of frames within 4 Å)

The verdict logic (`classify_stability`) operates on the final values of those
time series; thresholds are calibrated to standard MD conventions and live in
one place so future calibration is a single-line change.

HETATM classification (parse_pdb): Distinguishes the docked drug ligand from
bound cofactors and ions. Without this distinction, RMSD metrics would lump
the drug + GTP + GDP + metal ions into a single 'ligand' cloud and produce
meaningless combined-cloud Kabsch RMSD values (each component is on an
independent rigid-body trajectory). The exclusion list shares `_METAL_IONS`
and `_COFACTOR_SMILES` with the receptor-prep pipeline
(`backend/app/docking/md_receptor_prep`) to maintain a single source of
truth for what counts as a cofactor vs the docked compound.

Per diagnostic 2026-06-03: prior versions reported ligand RMSD = 4.94 Å on a
taxol back-dock that, after correct classification, has true taxol-only
Kabsch internal RMSD = 0.65 Å and receptor-aligned pose displacement = 1.5 Å
— both consistent with a stable in-pocket pose.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from backend.app.utils.rmsd import kabsch_fit, kabsch_rmsd
from backend.app.md import (
    VERDICT_STABLE,
    VERDICT_DRIFTING,
    VERDICT_UNSTABLE,
)
from backend.app.docking.md_receptor_prep import _METAL_IONS, _COFACTOR_SMILES

HBOND_DONOR_ELEMENTS = {"N", "O", "S"}
HBOND_ACCEPTOR_ELEMENTS = {"N", "O", "S"}
HBOND_DISTANCE_A = 3.5
CONTACT_DISTANCE_A = 4.0

# Waters carry various resnames depending on the upstream pipeline; track all
# common ones so a model converted via any tool gets handled.
_WATERS = {"HOH", "WAT", "H2O", "TIP3", "TIP", "SOL"}

# Composite set of HETATM resnames to EXCLUDE from both receptor and ligand
# atom lists in parse_pdb. Cofactors (parameterized via GAFF at MD time but
# not the docked drug), metal ions (single-atom residues coordinated to the
# protein), and waters. The docked drug's HETATM (typically `LIG`, `UNL`, or
# `MOL` depending on the engine) is whatever's left over after this exclusion.
_NON_LIGAND_HETATMS = (
    set(_METAL_IONS) |
    set(_COFACTOR_SMILES.keys()) |
    _WATERS
)


@dataclass
class Snapshot:
    receptor_atoms: List[dict]   # {chain, resseq, resname, name, element, xyz}
    ligand_atoms: List[dict]


def parse_pdb(path: Path) -> Snapshot:
    receptor: List[dict] = []
    ligand: List[dict] = []
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
            chain = line[21:22]
            try:
                resseq = int(line[22:26])
            except ValueError:
                resseq = 0
            element = line[76:78].strip()
            if not element:
                # Fall back to the leading non-digit char of the atom name.
                element = "".join(c for c in name if not c.isdigit())[:1].upper()
            row = {
                "chain": chain,
                "resseq": resseq,
                "resname": resname,
                "name": name,
                "element": element,
                "xyz": (x, y, z),
            }
            if tag == "ATOM  ":
                receptor.append(row)
            elif resname in _NON_LIGAND_HETATMS:
                # Cofactor / metal / water — neither receptor nor docked ligand
                # for the purposes of RMSD or contact analysis. Drops these
                # from both lists. See module docstring for the full rationale.
                continue
            else:
                ligand.append(row)
    return Snapshot(receptor_atoms=receptor, ligand_atoms=ligand)


# ---------------------------------------------------------------------------
# RMSD time series
# ---------------------------------------------------------------------------
def _backbone_coords(snap: Snapshot) -> np.ndarray:
    pts = [a["xyz"] for a in snap.receptor_atoms if a["name"] == "CA"]
    return np.asarray(pts, dtype=float) if pts else np.empty((0, 3))


def _ligand_heavy_coords(snap: Snapshot) -> np.ndarray:
    pts = [a["xyz"] for a in snap.ligand_atoms if a["element"] != "H"]
    return np.asarray(pts, dtype=float) if pts else np.empty((0, 3))


def compute_rmsd_series(
    snapshots: Sequence[Snapshot],
    times_ps: Sequence[float],
) -> List[Tuple[float, float, float, float]]:
    """Returns [(time_ps, bb_ca_rmsd, lig_internal_rmsd, lig_pose_rmsd), ...]
    against snapshot 0.

    Q6b reframing:
      - `lig_internal_rmsd` is the pre-Q6b ligand-on-ligand Kabsch RMSD — it
        measures conformational change of the ligand itself (sliding bond
        rotations, torsion noise) but is blind to whether the ligand stayed
        in the pocket. Kept as a secondary diagnostic.
      - `lig_pose_rmsd` is the post-Q6b POSE RMSD: align each frame's
        receptor backbone Cα onto the reference receptor backbone Cα to get
        the rigid-body (R, t), apply that SAME transform to the frame's
        ligand heavy atoms, then compute plain RMSD vs the reference ligand.
        No second alignment of the ligand. This is the metric that answers
        "did the ligand stay in the pocket?" — and is what classify_stability
        keys the verdict on from Q6b onwards.

    Atom count drift between snapshot 0 and any later frame propagates NaN
    for the affected metric — same convention as the pre-Q6b code.
    """
    if not snapshots:
        return []
    ref_bb = _backbone_coords(snapshots[0])
    ref_lig = _ligand_heavy_coords(snapshots[0])
    out: List[Tuple[float, float, float, float]] = []
    for t, snap in zip(times_ps, snapshots):
        bb = _backbone_coords(snap)
        lig = _ligand_heavy_coords(snap)

        # Backbone RMSD (independent receptor superposition; reported as-is).
        bb_rmsd = (
            kabsch_rmsd(ref_bb, bb)
            if ref_bb.shape == bb.shape and ref_bb.size > 0
            else float("nan")
        )

        # Ligand internal RMSD (ligand-on-ligand). Same value the pre-Q6b
        # code reported; kept verbatim under the new name.
        lig_internal = (
            kabsch_rmsd(ref_lig, lig)
            if ref_lig.shape == lig.shape and ref_lig.size > 0
            else float("nan")
        )

        # Ligand POSE RMSD: receptor-frame superposition, then plain RMSD.
        # Requires both the backbone and the ligand to be shape-consistent
        # with the reference; otherwise NaN.
        lig_pose = float("nan")
        if (
            ref_bb.shape == bb.shape
            and ref_bb.size > 0
            and ref_lig.shape == lig.shape
            and ref_lig.size > 0
        ):
            R, tvec, _ = kabsch_fit(ref_bb, bb)
            lig_aligned = lig @ R + tvec
            diff = ref_lig - lig_aligned
            lig_pose = float(
                np.sqrt((diff * diff).sum() / ref_lig.shape[0])
            )

        out.append((float(t), float(bb_rmsd), float(lig_internal), lig_pose))
    return out


def write_rmsd_csv(
    path: Path,
    series: Sequence[Tuple[float, float, float, float]],
) -> None:
    """Q6b: adds the `rmsd_ligand_pose_A` column. CSV column order documents
    the post-Q6b primary/secondary distinction — pose first, then internal."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "time_ps",
            "rmsd_backbone_ca_A",
            "rmsd_ligand_pose_A",
            "rmsd_ligand_internal_A",
        ])
        for t, bb, lig_internal, lig_pose in series:
            w.writerow([
                f"{t:.3f}",
                f"{bb:.3f}" if math.isfinite(bb) else "",
                f"{lig_pose:.3f}" if math.isfinite(lig_pose) else "",
                f"{lig_internal:.3f}" if math.isfinite(lig_internal) else "",
            ])


# ---------------------------------------------------------------------------
# H-bonds (count per snapshot; donor/acceptor heavy-atom distance heuristic)
# ---------------------------------------------------------------------------
def _polar_atoms(atoms: List[dict]) -> List[Tuple[Tuple[float, float, float], str, dict]]:
    out = []
    for a in atoms:
        if a["element"] in HBOND_DONOR_ELEMENTS or a["element"] in HBOND_ACCEPTOR_ELEMENTS:
            out.append((a["xyz"], a["element"], a))
    return out


def count_hbonds(snap: Snapshot) -> int:
    """Permissive heavy-atom distance heuristic (no angle check, since trajectories
    rarely retain hydrogens in PDB snapshots). Conservative threshold to avoid
    over-counting."""
    rec = _polar_atoms(snap.receptor_atoms)
    lig = _polar_atoms(snap.ligand_atoms)
    if not rec or not lig:
        return 0
    rec_xyz = np.asarray([p[0] for p in rec], dtype=float)
    lig_xyz = np.asarray([p[0] for p in lig], dtype=float)
    diff = rec_xyz[:, None, :] - lig_xyz[None, :, :]
    d2 = (diff * diff).sum(axis=2)
    return int(np.sum(d2 < (HBOND_DISTANCE_A ** 2)))


def write_hbonds_csv(
    path: Path,
    counts: Sequence[Tuple[float, int]],
    persistence: Dict[str, float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_ps", "hbond_count"])
        for t, c in counts:
            w.writerow([f"{t:.3f}", c])
        w.writerow([])
        w.writerow(["#summary"])
        w.writerow(["mean_count", f"{persistence.get('mean', 0.0):.3f}"])
        w.writerow(["min_count", str(int(persistence.get("min", 0)))])
        w.writerow(["max_count", str(int(persistence.get("max", 0)))])
        w.writerow(["frac_frames_with_any", f"{persistence.get('frac_with_any', 0.0):.3f}"])


# ---------------------------------------------------------------------------
# Per-residue contact frequency
# ---------------------------------------------------------------------------
def per_residue_contact_frequency(
    snapshots: Sequence[Snapshot],
    cutoff_a: float = CONTACT_DISTANCE_A,
) -> List[Tuple[str, int, str, float]]:
    """Returns (chain, resseq, resname, fraction_of_frames_with_contact)."""
    if not snapshots:
        return []
    n = len(snapshots)
    counts: Dict[Tuple[str, int, str], int] = {}
    for snap in snapshots:
        lig_xyz = np.asarray(
            [a["xyz"] for a in snap.ligand_atoms if a["element"] != "H"],
            dtype=float,
        )
        if lig_xyz.size == 0:
            continue
        seen: set = set()
        # Bucket receptor atoms by residue key.
        for a in snap.receptor_atoms:
            if a["element"] == "H":
                continue
            key = (a["chain"], a["resseq"], a["resname"])
            xyz = np.asarray(a["xyz"], dtype=float)
            d2 = ((lig_xyz - xyz) ** 2).sum(axis=1)
            if np.any(d2 < cutoff_a * cutoff_a):
                seen.add(key)
        for key in seen:
            counts[key] = counts.get(key, 0) + 1
    out = [(k[0], k[1], k[2], v / n) for k, v in counts.items()]
    out.sort(key=lambda x: -x[3])
    return out


def write_contacts_csv(
    path: Path,
    rows: Sequence[Tuple[str, int, str, float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chain", "resseq", "resname", "frac_frames_in_contact"])
        for ch, rs, rn, frac in rows:
            w.writerow([ch, rs, rn, f"{frac:.3f}"])


# ---------------------------------------------------------------------------
# Q6b PART 2 — residue-number relabel (MD → docking author numbering)
# ---------------------------------------------------------------------------
def _chain_ca_sequence(path: Path) -> Dict[str, List[Tuple[int, str]]]:
    """Walk a PDB; return {chain_id: [(resseq, resname), ...]} in file order,
    keeping only the first Cα per residue (PDB residues normally have one
    Cα). Used by the relabel map builder."""
    out: Dict[str, List[Tuple[int, str]]] = {}
    seen: set = set()
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21:22]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            resname = line[17:20].strip()
            key = (chain, resseq)
            if key in seen:
                continue
            seen.add(key)
            out.setdefault(chain, []).append((resseq, resname))
    return out


def build_md_to_docking_resseq_map(
    md_pdb: Path,
    docking_pdb: Path,
) -> Optional[Dict[Tuple[str, int], int]]:
    """Build a per-chain (chain, md_resseq) → docking_resseq map.

    The MD pipeline renumbers chains to be contiguous from 1; the docking
    pipeline preserves author numbering (with explicit gaps where the source
    PDB had missing-density residues). The contacts pipeline labels
    interactions in AUTHOR numbering — which is what the literature uses
    (PHE272, LEU275, GLY370 on β-tubulin for the taxane site).

    Map construction: per-chain, we expect the SAME number of Cα atoms in
    both files (the prep dropped resseq markers, not residues), and the
    residue NAMES at corresponding positions must match. The MD i-th
    chain-B residue (1-indexed resseq i) maps to the docking i-th chain-B
    residue (author resseq).

    PER-CHAIN SAFETY: a length / residue-name disagreement on ONE chain
    drops that chain from the map but does NOT abort the whole map. Other
    chains whose alignment is clean still get relabeled. Rows for the
    disagreement chain pass through `relabel_contacts` unchanged
    (preserving MD numbering — honest, not silent-wrong).

    Returns None only when neither chain can be aligned (e.g. the two
    files describe entirely different receptors) — in which case the
    caller skips the relabel and surfaces `receptor_renumbered=False` in
    the artifact.
    """
    try:
        md_seq = _chain_ca_sequence(md_pdb)
        dock_seq = _chain_ca_sequence(docking_pdb)
    except Exception:
        return None

    mapping: Dict[Tuple[str, int], int] = {}
    aligned_chains: List[str] = []
    for chain, md_rows in md_seq.items():
        dock_rows = dock_seq.get(chain)
        if dock_rows is None:
            # Chain only present in MD — skip this chain's relabel.
            continue
        if len(dock_rows) != len(md_rows):
            # Per-chain length mismatch — sequence alignment would be
            # needed; outside Q6b's clean-relabel guard. Skip this
            # chain; other chains may still align cleanly.
            continue
        chain_ok = True
        partial: Dict[Tuple[str, int], int] = {}
        for (md_rs, md_name), (dock_rs, dock_name) in zip(md_rows, dock_rows):
            if md_name != dock_name:
                # Residue identity mismatch at this position → the
                # two files describe different chemistry on this
                # chain; skip the whole chain (don't half-relabel).
                chain_ok = False
                break
            partial[(chain, md_rs)] = dock_rs
        if chain_ok:
            mapping.update(partial)
            aligned_chains.append(chain)

    # If NO chain aligned, surface as None so callers can flag
    # receptor_renumbered=False. If at least one chain aligned, return
    # the partial map — the unaligned chains' contacts pass through
    # unchanged in relabel_contacts.
    return mapping if aligned_chains else None


def relabel_contacts(
    rows: Sequence[Tuple[str, int, str, float]],
    resseq_map: Optional[Dict[Tuple[str, int], int]],
) -> List[Tuple[str, int, str, float]]:
    """Apply a (chain, md_resseq) → docking_resseq map to a contacts list.
    Rows with no mapping are passed through unchanged (preserves the MD
    resseq for any chain not covered by the map). Idempotent on rows that
    are already author-numbered AND coincidentally collide with an
    md_resseq — that collision would only happen on a chain whose MD
    resseq == docking resseq for that residue, in which case the relabel
    is the identity for that row. Callers should NOT call relabel twice
    on the same rows in general."""
    if not resseq_map:
        return list(rows)
    out: List[Tuple[str, int, str, float]] = []
    for ch, rs, rn, frac in rows:
        new_rs = resseq_map.get((ch, rs), rs)
        out.append((ch, new_rs, rn, frac))
    return out


# ---------------------------------------------------------------------------
# Stability verdict
# ---------------------------------------------------------------------------
@dataclass
class StabilitySummary:
    """Q6b: verdict is now keyed on the POSE RMSD (receptor-frame ligand
    displacement) — the metric that actually answers "did the ligand stay
    in the pocket?". The pre-Q6b ligand-on-ligand RMSD is preserved
    verbatim as the *_internal_a fields, surfaced as a diagnostic so
    consumers can see ligand conformational change separately from pose
    drift. Thresholds (≤2 stable / 2-4 drifting / >4 unstable) are
    unchanged — the metric they read is now strictly stronger."""
    verdict: str
    rmsd_backbone_final_a: Optional[float]
    # Primary ligand metrics (Q6b): receptor-aligned pose displacement.
    rmsd_ligand_pose_final_a: Optional[float]
    rmsd_ligand_pose_max_a: Optional[float]
    # Secondary diagnostic (pre-Q6b semantics, kept under new name):
    # ligand-on-ligand Kabsch RMSD — internal conformational change.
    rmsd_ligand_internal_final_a: Optional[float]
    rmsd_ligand_internal_max_a: Optional[float]
    hbond_persistence_frac: Optional[float]
    top_contacts: List[Tuple[str, int, str, float]]
    rationale: str


def classify_stability(
    rmsd_series: Sequence[Tuple[float, float, float, float]],
    hbond_counts: Sequence[Tuple[float, int]],
    contacts: Sequence[Tuple[str, int, str, float]],
) -> StabilitySummary:
    """Q6b: keys the verdict on `rmsd_ligand_pose_*` (receptor-frame
    superposition), not on the legacy ligand-on-ligand value. Some pre-Q6b
    "stable" jobs will reclassify under the stronger metric — that is the
    metric gaining discriminating power, not a regression."""
    if not rmsd_series:
        return StabilitySummary(
            verdict=VERDICT_UNSTABLE,
            rmsd_backbone_final_a=None,
            rmsd_ligand_pose_final_a=None,
            rmsd_ligand_pose_max_a=None,
            rmsd_ligand_internal_final_a=None,
            rmsd_ligand_internal_max_a=None,
            hbond_persistence_frac=None,
            top_contacts=list(contacts[:10]),
            rationale="No frames produced.",
        )
    # Series rows: (time_ps, bb_rmsd, lig_internal, lig_pose).
    bb_final = rmsd_series[-1][1]
    lig_internal_final = rmsd_series[-1][2]
    lig_pose_final = rmsd_series[-1][3]
    lig_internal_max = max(
        (row[2] for row in rmsd_series if math.isfinite(row[2])),
        default=float("nan"),
    )
    lig_pose_max = max(
        (row[3] for row in rmsd_series if math.isfinite(row[3])),
        default=float("nan"),
    )
    persist = (
        sum(1 for _, c in hbond_counts if c > 0) / len(hbond_counts)
        if hbond_counts else 0.0
    )

    # Thresholds (standard MD convention; ligand-driven verdict).
    # Backbone RMSD up to ~5 Å is normal equilibrium fluctuation for a large
    # protein in implicit solvent (e.g. 3.1 Å for an 877-residue tubulin
    # dimer); it informs the rationale but does not independently downgrade
    # an otherwise-stable pose. Q6b: the ligand metric below is now pose
    # RMSD (receptor-frame displacement), not internal RMSD.
    if math.isfinite(lig_pose_final) and lig_pose_final > 4.0:
        verdict = VERDICT_UNSTABLE
        rationale = (
            f"Ligand pose RMSD {lig_pose_final:.2f} Å — likely unbinding "
            f"(internal {lig_internal_final:.2f} Å)."
        )
    elif math.isfinite(lig_pose_final) and lig_pose_final > 2.0:
        verdict = VERDICT_DRIFTING
        rationale = (
            f"Ligand pose RMSD {lig_pose_final:.2f} Å — pose is sliding but "
            f"still bound (internal {lig_internal_final:.2f} Å)."
        )
    elif math.isfinite(bb_final) and bb_final > 5.0:
        verdict = VERDICT_DRIFTING
        rationale = (
            f"Backbone moved {bb_final:.2f} Å — protein not equilibrated; "
            f"ligand pose RMSD {lig_pose_final:.2f} Å is acceptable but the "
            f"system may need longer equilibration."
        )
    else:
        verdict = VERDICT_STABLE
        rationale = (
            f"Backbone RMSD {bb_final:.2f} Å, ligand pose RMSD "
            f"{lig_pose_final:.2f} Å (internal {lig_internal_final:.2f} Å); "
            f"H-bonds present in {persist*100:.0f}% of frames."
        )

    return StabilitySummary(
        verdict=verdict,
        rmsd_backbone_final_a=bb_final if math.isfinite(bb_final) else None,
        rmsd_ligand_pose_final_a=lig_pose_final if math.isfinite(lig_pose_final) else None,
        rmsd_ligand_pose_max_a=lig_pose_max if math.isfinite(lig_pose_max) else None,
        rmsd_ligand_internal_final_a=lig_internal_final if math.isfinite(lig_internal_final) else None,
        rmsd_ligand_internal_max_a=lig_internal_max if math.isfinite(lig_internal_max) else None,
        hbond_persistence_frac=persist,
        top_contacts=list(contacts[:10]),
        rationale=rationale,
    )
