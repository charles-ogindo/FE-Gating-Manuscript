"""Single-trajectory MM-GBSA free-energy estimator on OpenMM.

Computes ΔG_bind ≈ <E_complex> − <E_receptor> − <E_ligand>
on a single MD trajectory where receptor and ligand subsystems are
extracted from each complex frame (no separate simulations). Entropy
is OMITTED — configurational entropy (normal-mode or quasi-harmonic)
is out of scope for this estimator; the reported number is the
enthalpic + solvation part of the binding free energy, which is the
standard MM-GBSA quantity in the literature.

Force-field reuse policy
------------------------
The estimator REUSES the exact MD parameterization (amber14-all +
amber14/tip3p ion templates + implicit/obc2 + gaff-2.11) so the
single-trajectory subtraction is consistent. The MD itself was already
implicit-solvent OBC2 — no explicit waters/ions to strip.

Per-frame components
--------------------
The estimator separates by OpenMM force-group dispatch:
  * `bonded`     = HarmonicBond + HarmonicAngle + PeriodicTorsion
  * `nonbonded`  = NonbondedForce (vacuum LJ + Coulomb under
                   CutoffNonPeriodic; the SAME force class as the MD)
  * `solvation`  = GBSAOBCForce (polar GB + SA nonpolar; the OBC2
                   implicit-solvent contribution)
ΔG component = <E_component_complex − E_component_receptor − E_component_ligand>.

Sampling-adequacy
-----------------
The estimator is run BEHIND the existing free_energy/gating.py gate.
It NEVER lowers MIN_TOTAL_FRAMES to make a job "pass" — that threshold
is the integrity of the gating discipline. When the gate blocks (e.g.
taxol's 101 < 200), the artifact is marked `sampling_adequate=false`
and `preliminary=true`, the gate reason is surfaced verbatim, and the
ΔG is still computed (an honestly-flagged number is the goal). When
the gate passes, `sampling_adequate=true` and the number stands as
the project's standard MM-GBSA estimate.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

# Schema bump should accompany any backwards-incompatible artifact change.
FE_ARTIFACT_SCHEMA_VERSION = "1.0.0"


# ---------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------

@dataclass
class FrameEnergy:
    """Per-frame breakdown. All values in kcal/mol after conversion
    from OpenMM's native kJ/mol."""
    frame_index: int
    e_complex_total: float
    e_receptor_total: float
    e_ligand_total: float
    delta_g_total: float
    components: Dict[str, float] = field(default_factory=dict)
    # components keys: bonded, nonbonded, solvation  (each is the Δ value
    # = complex − receptor − ligand for that force group)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


@dataclass
class MMGBSAResult:
    """Aggregate output. mean ± SEM in kcal/mol; n is the number of
    frames that actually evaluated (skipped frames excluded)."""
    n_frames_used: int
    n_frames_skipped: int
    delta_g_mean: float
    delta_g_sem: float
    delta_g_stddev: float
    components_mean: Dict[str, float]
    components_sem: Dict[str, float]
    per_frame: List[FrameEnergy]
    method: Dict[str, Any]
    sampling_adequate: bool
    gate_reason: Optional[str]
    preliminary: bool

    def to_dict(self) -> dict:
        return {
            "n_frames_used": self.n_frames_used,
            "n_frames_skipped": self.n_frames_skipped,
            "delta_g_mean_kcal_per_mol": self.delta_g_mean,
            "delta_g_sem_kcal_per_mol":  self.delta_g_sem,
            "delta_g_stddev_kcal_per_mol": self.delta_g_stddev,
            "components_mean_kcal_per_mol": self.components_mean,
            "components_sem_kcal_per_mol":  self.components_sem,
            "method": self.method,
            "sampling_adequate": self.sampling_adequate,
            "gate_reason": self.gate_reason,
            "preliminary": self.preliminary,
            "per_frame": [fe.to_dict() for fe in self.per_frame],
        }


# ---------------------------------------------------------------------
# Force-group dispatch — split E_total into bonded / nonbonded / solvation.
# ---------------------------------------------------------------------

# Map each OpenMM force class name to a string component bucket. The
# class names come from `force.__class__.__name__`; OpenMM's stable
# Python wrappers keep these stable across 8.x.
_FORCE_GROUP_BUCKET = {
    "HarmonicBondForce":      "bonded",
    "HarmonicAngleForce":     "bonded",
    "PeriodicTorsionForce":   "bonded",
    "RBTorsionForce":         "bonded",
    "CMAPTorsionForce":       "bonded",
    "NonbondedForce":         "nonbonded",
    "CustomNonbondedForce":   "nonbonded",
    "GBSAOBCForce":           "solvation",
    "CustomGBForce":          "solvation",
    # Misc. constraint/center-of-mass forces report 0 in single-point
    # energies under typical settings; bucket as "other" so an
    # unexpected class doesn't silently corrupt a known component.
}

_DEFAULT_GROUP_FOR_BUCKET = {"bonded": 1, "nonbonded": 2, "solvation": 3, "other": 4}


def _assign_force_groups(system) -> Dict[int, str]:
    """Assign each force in `system` to a force group keyed by component
    bucket. Returns {group_int -> bucket_name} so the caller can read
    per-component energies via getState(groups={N})."""
    assigned: Dict[int, str] = {}
    for force in system.getForces():
        bucket = _FORCE_GROUP_BUCKET.get(force.__class__.__name__, "other")
        group = _DEFAULT_GROUP_FOR_BUCKET[bucket]
        force.setForceGroup(group)
        assigned[group] = bucket
    return assigned


# ---------------------------------------------------------------------
# Subsystem extraction — given a complex Modeller, derive receptor and
# ligand Modellers by topology subset. Positions slice along the same
# subset. This is the single-trajectory contract: same positions, three
# topology slices.
# ---------------------------------------------------------------------

def _split_modeller(complex_modeller, ligand_residue_name: str = "LIG"):
    """Return (receptor_modeller, ligand_modeller). The split is by
    residue name — `ligand_residue_name` (default "LIG") goes to the
    ligand; everything else goes to the receptor (cofactors + ions +
    protein). The engine writes the docked ligand as residue "LIG"
    in `_openmm_ligand_with_h.pdb` (`engine_openmm.py:133`), so this
    default matches the MD's own labeling.
    """
    from openmm import app
    import numpy as np

    top = complex_modeller.topology
    pos = complex_modeller.positions

    # Build atom-index masks.
    receptor_idx: List[int] = []
    ligand_idx: List[int] = []
    for atom in top.atoms():
        if atom.residue.name == ligand_residue_name:
            ligand_idx.append(atom.index)
        else:
            receptor_idx.append(atom.index)

    if not ligand_idx:
        raise ValueError(
            f"No atoms found with residue name {ligand_residue_name!r} — "
            "cannot extract ligand subsystem"
        )
    if not receptor_idx:
        raise ValueError("No receptor atoms — entire complex is the ligand")

    # Modeller has delete(atoms) — pass the atoms to KEEP-NOT, i.e.,
    # we delete the complement. Do this on copies so the original is
    # untouched.
    receptor_modeller = app.Modeller(top, pos)
    receptor_modeller.delete([
        a for a in receptor_modeller.topology.atoms()
        if a.residue.name == ligand_residue_name
    ])

    ligand_modeller = app.Modeller(top, pos)
    ligand_modeller.delete([
        a for a in ligand_modeller.topology.atoms()
        if a.residue.name != ligand_residue_name
    ])

    return receptor_modeller, ligand_modeller, receptor_idx, ligand_idx


# ---------------------------------------------------------------------
# Energy evaluation
# ---------------------------------------------------------------------

def _build_context(system, topology, *, platform_name: Optional[str] = None):
    """Build ONE OpenMM Context for a System. The caller drives many
    single-point evals against the same Context via setPositions() —
    Context construction is the expensive step for a ~13k-atom OBC2
    system on CPU (it allocates per-platform buffers, builds the
    neighborlist, compiles the GBSA kernel, validates the topology
    against the system). Reusing one Context across all frames cuts
    303 builds → 3 on the taxol smoke (~3.7h → ~minutes).

    Returns the openmm.Context. The integrator is wired to the Context
    but never stepped — we only ever call setPositions + getState.

    The caller MUST have already called `_assign_force_groups(system)`
    before this; once a Context exists, force-group assignments on the
    underlying System aren't re-read on subsequent getState(groups=...)
    calls.
    """
    from openmm import LangevinMiddleIntegrator, Platform, unit, app

    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        2.0 * unit.femtosecond,
    )
    if platform_name is None:
        try:
            platform = Platform.getPlatformByName("CPU")
        except Exception:
            platform = None
    else:
        platform = Platform.getPlatformByName(platform_name)

    if platform is not None:
        sim = app.Simulation(topology, system, integrator, platform)
    else:
        sim = app.Simulation(topology, system, integrator)
    # Stash the System-derived group→bucket map on the Context so the
    # per-frame eval helper doesn't have to re-walk the forces 101 times.
    groups_to_bucket: Dict[int, str] = {}
    for force in system.getForces():
        g = force.getForceGroup()
        bucket = _FORCE_GROUP_BUCKET.get(force.__class__.__name__, "other")
        groups_to_bucket[g] = bucket
    sim.context._mmgbsa_groups_to_bucket = groups_to_bucket
    return sim.context


def _evaluate_context(context, positions) -> Tuple[float, Dict[str, float]]:
    """Set positions on an existing Context, return
    (total_kj, per-bucket_kj_dict). No Context construction — that's
    why this is fast."""
    from openmm import unit

    context.setPositions(positions)

    state_all = context.getState(getEnergy=True)
    total_kj = state_all.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    per_bucket: Dict[str, float] = {"bonded": 0.0, "nonbonded": 0.0, "solvation": 0.0, "other": 0.0}
    groups_to_bucket: Dict[int, str] = getattr(context, "_mmgbsa_groups_to_bucket", {})
    for group, bucket in groups_to_bucket.items():
        st = context.getState(getEnergy=True, groups={group})
        per_bucket[bucket] += st.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )

    return total_kj, per_bucket


_KJ_PER_KCAL = 4.184


def _kj_to_kcal(x: float) -> float:
    return x / _KJ_PER_KCAL


# ---------------------------------------------------------------------
# Estimator — public entry point
# ---------------------------------------------------------------------

def estimate_mmgbsa(
    *,
    complex_system,
    complex_topology,
    receptor_system,
    receptor_topology,
    ligand_system,
    ligand_topology,
    frame_positions: Sequence,  # iterable of OpenMM positions or numpy arrays
    receptor_idx: Sequence[int],
    ligand_idx: Sequence[int],
    gate_can_run: bool,
    gate_reason: Optional[str],
    method_meta: Dict[str, Any],
) -> MMGBSAResult:
    """Run the single-trajectory MM-GBSA estimator.

    All three systems must be built consistently (same FF, same OBC2,
    same nonbonded settings). Each frame's positions are sliced by
    `receptor_idx` and `ligand_idx` for the two subsystem evaluations.

    `gate_can_run` and `gate_reason` flow from free_energy/gating.py's
    validate() — when False, the result is marked preliminary +
    sampling_adequate=False but the numbers are STILL computed.
    """
    from openmm import unit
    import numpy as np

    # Force-group assignment must happen BEFORE Context construction.
    # OpenMM caches the group integer at Context build time; mutating
    # it post-hoc on the System would not propagate to the live
    # Contexts. The slow path called _single_point_energy → which
    # built a fresh Context each frame → which re-picked up the
    # mutation; the fast path builds Contexts once and depends on
    # this ordering.
    _assign_force_groups(complex_system)
    _assign_force_groups(receptor_system)
    _assign_force_groups(ligand_system)

    # Build ONE Context per subsystem. 303 → 3 in the per-frame loop.
    # This is the entire fast-path win; the actual energy math is
    # unchanged from the slow path.
    cpx_ctx = _build_context(complex_system, complex_topology)
    rec_ctx = _build_context(receptor_system, receptor_topology)
    lig_ctx = _build_context(ligand_system, ligand_topology)

    receptor_idx = list(receptor_idx)
    ligand_idx = list(ligand_idx)

    frames: List[FrameEnergy] = []
    n_skipped = 0
    for i, pos in enumerate(frame_positions):
        try:
            # Slice positions for the subsystems. `pos` may be an
            # OpenMM Quantity-wrapped array or a plain numpy array
            # (nm). We handle both.
            if hasattr(pos, "value_in_unit"):
                pos_nm = pos.value_in_unit(unit.nanometer)
                pos_array = np.asarray(pos_nm)
            else:
                pos_array = np.asarray(pos)

            recv_pos = pos_array[receptor_idx, :] * unit.nanometer
            lig_pos = pos_array[ligand_idx, :] * unit.nanometer
            cpx_pos = pos_array * unit.nanometer

            e_cpx_total, e_cpx_buckets = _evaluate_context(cpx_ctx, cpx_pos)
            e_rec_total, e_rec_buckets = _evaluate_context(rec_ctx, recv_pos)
            e_lig_total, e_lig_buckets = _evaluate_context(lig_ctx, lig_pos)

            delta_g_total_kj = e_cpx_total - e_rec_total - e_lig_total
            components_kcal = {}
            for bucket in ("bonded", "nonbonded", "solvation"):
                d_kj = e_cpx_buckets.get(bucket, 0.0) \
                    - e_rec_buckets.get(bucket, 0.0) \
                    - e_lig_buckets.get(bucket, 0.0)
                components_kcal[bucket] = _kj_to_kcal(d_kj)

            frames.append(FrameEnergy(
                frame_index=i,
                e_complex_total=_kj_to_kcal(e_cpx_total),
                e_receptor_total=_kj_to_kcal(e_rec_total),
                e_ligand_total=_kj_to_kcal(e_lig_total),
                delta_g_total=_kj_to_kcal(delta_g_total_kj),
                components=components_kcal,
            ))
        except Exception as e:
            logger.warning("frame %d failed energy evaluation: %s", i, e)
            n_skipped += 1
            continue

    if not frames:
        raise RuntimeError("No frames evaluated successfully")

    # Aggregate
    deltas = [f.delta_g_total for f in frames]
    n = len(deltas)
    mean_dg = sum(deltas) / n
    var = sum((d - mean_dg) ** 2 for d in deltas) / max(n - 1, 1)
    stddev = math.sqrt(var)
    sem = stddev / math.sqrt(n) if n > 1 else float("nan")

    comp_means: Dict[str, float] = {}
    comp_sems: Dict[str, float] = {}
    for key in ("bonded", "nonbonded", "solvation"):
        vals = [f.components.get(key, 0.0) for f in frames]
        m = sum(vals) / n
        v = sum((x - m) ** 2 for x in vals) / max(n - 1, 1)
        s = math.sqrt(v)
        comp_means[key] = m
        comp_sems[key] = (s / math.sqrt(n)) if n > 1 else float("nan")

    return MMGBSAResult(
        n_frames_used=n,
        n_frames_skipped=n_skipped,
        delta_g_mean=mean_dg,
        delta_g_sem=sem,
        delta_g_stddev=stddev,
        components_mean=comp_means,
        components_sem=comp_sems,
        per_frame=frames,
        method=method_meta,
        sampling_adequate=bool(gate_can_run),
        gate_reason=gate_reason if not gate_can_run else None,
        preliminary=not bool(gate_can_run),
    )
