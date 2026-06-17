"""End-to-end MM-GBSA runner for an MD job: gate → systems → frames → estimator
→ self-describing free_energy artifact.

Two orchestration entry points:

  - `prepare_mmgbsa_systems(md_dir, ...)` — builds (or loads from cache) the
    three OpenMM Systems consistent with the MD parameterization. The 3-System
    build runs `sqm` (AmberTools AM1-BCC) on the small molecules, which takes
    1-4 hours per ligand cold — so the result is cached under
    `<md_dir>/_mmgbsa_cache/` as serialized System XML + a topology PDB. The
    cache key is a hash of (ligand_smiles, cofactor SMILES list, force-field
    set) so a SMILES change OR a force-field bump invalidates it.

  - `compute_md_fe(md_id, ...)` — runs the gate + systems + estimator and
    returns the self-describing `free_energy` artifact block, ready to be
    merged into `summary.json`. Honors `gate_can_run=False` (still computes
    but marks the result preliminary; does NOT lower MIN_TOTAL_FRAMES — that
    threshold is the integrity of the gating discipline).

Single-trajectory contract: receptor and ligand subsystems share positions
with the complex frame; no separate sims; entropy omitted (normal-mode
out of scope — flagged in method.configurational_entropy).
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from backend.app.core.config import JOBS_DIR
from backend.app.free_energy.gating import validate
from backend.app.free_energy.mmgbsa import (
    FE_ARTIFACT_SCHEMA_VERSION,
    MMGBSAResult,
    _split_modeller,
    estimate_mmgbsa,
)
from backend.app.md.io import read_ligand_smiles

logger = logging.getLogger(__name__)

# Forcefield identifiers — kept in one place so the cache key stays in sync
# with whatever the MD engine used. If `engine_openmm.py` is ever upgraded
# (e.g., GAFF version bump), bump these too OR the cache will return a stale
# parameterization for a newly-rebuilt MD.
_FF_PROTEIN_AND_IONS = ["amber14-all.xml", "amber14/tip3p.xml"]
_FF_IMPLICIT = "implicit/obc2.xml"
_FF_SMALL_MOL = "gaff-2.11"


# ---------------------------------------------------------------------
# System cache — XML-serialized Systems + topology PDB under
# <md_dir>/_mmgbsa_cache/. Cache key = SHA256(JSON(ligand_smiles +
# cofactor SMILES list + FF set)).
# ---------------------------------------------------------------------

def _cache_dir(md_dir: Path) -> Path:
    return md_dir / "_mmgbsa_cache"


def _cache_key(ligand_smiles: str, cofactor_smiles: List[str]) -> str:
    payload = {
        "ligand_smiles": ligand_smiles,
        "cofactor_smiles": sorted(cofactor_smiles),
        "ff_protein_ions": _FF_PROTEIN_AND_IONS,
        "ff_implicit": _FF_IMPLICIT,
        "ff_small_mol": _FF_SMALL_MOL,
    }
    blob = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def _write_cache(
    cache: Path,
    cache_key: str,
    cpx_sys, rec_sys, lig_sys,
    cpx_topology, rec_topology, lig_topology,
    receptor_idx: Sequence[int], ligand_idx: Sequence[int],
) -> None:
    """Serialize 3 Systems + 3 topologies + atom-index map under cache/."""
    from openmm import XmlSerializer, app

    cache.mkdir(parents=True, exist_ok=True)

    # System XMLs
    for name, sys_obj in (
        ("complex_system.xml", cpx_sys),
        ("receptor_system.xml", rec_sys),
        ("ligand_system.xml", lig_sys),
    ):
        (cache / name).write_text(XmlSerializer.serialize(sys_obj))

    # Topologies — PDBFile.writeFile needs positions; write zero coords
    # since we only need the topology shape on reload.
    import numpy as np
    from openmm import unit

    def _write_topology_pdb(top, path: Path):
        n = sum(1 for _ in top.atoms())
        zero_pos = np.zeros((n, 3)) * unit.nanometer
        with path.open("w") as f:
            app.PDBFile.writeFile(top, zero_pos, f)

    _write_topology_pdb(cpx_topology, cache / "complex_topology.pdb")
    _write_topology_pdb(rec_topology, cache / "receptor_topology.pdb")
    _write_topology_pdb(lig_topology, cache / "ligand_topology.pdb")

    (cache / "cache_meta.json").write_text(json.dumps({
        "cache_key": cache_key,
        "schema_version": "1.0.0",
        "force_fields": _FF_PROTEIN_AND_IONS + [_FF_IMPLICIT],
        "small_molecule_forcefield": _FF_SMALL_MOL,
        "receptor_idx": list(receptor_idx),
        "ligand_idx": list(ligand_idx),
    }, indent=2))


def _try_load_cache(cache: Path, cache_key: str):
    """Return (cpx_sys, rec_sys, lig_sys, cpx_top, rec_top, lig_top,
    receptor_idx, ligand_idx) if cache matches key, else None."""
    meta_path = cache / "cache_meta.json"
    if not meta_path.is_file():
        return None
    try:
        meta = json.loads(meta_path.read_text())
        if meta.get("cache_key") != cache_key:
            return None
        from openmm import XmlSerializer, app

        def _load_sys(name):
            return XmlSerializer.deserializeSystem((cache / name).read_text())

        def _load_top(name):
            return app.PDBFile(str(cache / name)).topology

        return (
            _load_sys("complex_system.xml"),
            _load_sys("receptor_system.xml"),
            _load_sys("ligand_system.xml"),
            _load_top("complex_topology.pdb"),
            _load_top("receptor_topology.pdb"),
            _load_top("ligand_topology.pdb"),
            list(meta["receptor_idx"]),
            list(meta["ligand_idx"]),
        )
    except Exception:
        logger.exception("mmgbsa cache load failed; will rebuild")
        return None


# ---------------------------------------------------------------------
# System construction — mirrors engine_openmm.py's parameterization
# choices so the estimator is consistent with the MD it consumes.
# ---------------------------------------------------------------------

def _build_complex_modeller(md_dir: Path, docking_dir: Path, docking_job_id: str,
                            ligand_name: str):
    """Build the GBSA complex Modeller from the MD trajectory's frame_000.pdb.

    Topology is taken DIRECTLY from frame_000 — its (chain, resname, resseq,
    atom_name) tuples are the only ones the trajectory frames write, so the
    strip mask in `_collect_frame_positions` then matches 1:1 by construction.

    Previous implementation reconstructed the topology from
    `docking/viewer/receptor.pdb` + cofactor sidecars + the ligand-with-H
    sidecar and added them via `Modeller.add(...)`. That route drifted from
    what `engine_openmm` actually wrote into the frames:

      • N-terminal MET H atom: receptor.pdb has H1; frame writer renames
        to H.
      • Cofactor + ligand chain IDs: sidecars carry GTP on chain A:500,
        GDP on chain B:600, ligand on chain X:1 (Modeller.add may further
        reassign chain letters). The frame writer always lands MG on C,
        GTP on D, GDP on E, LIG on F regardless of sidecar origin.

    Net result was an atom-tuple mismatch that broke MM-GBSA on every real
    trajectory we tried (5/5 in the 2026-06-16 sweep). Using frame_000 as
    the topology source removes the entire class of writer-vs-builder
    convention drift — the strip mask cannot mismatch frame N if it was
    built from frame 0.

    `small_mols` are still built from SMILES (via the docking metadata +
    `read_ligand_smiles`) because openmmforcefields' GAFFTemplateGenerator
    matches topology residues to Molecule objects by graph isomorphism on
    the SMILES side — independent of which atoms-on-which-chain layout the
    topology uses.

    Returns (modeller, small_mols, ligand_residue_name, ligand_smiles,
    cofactor_smiles).
    """
    import openmm.app as app
    from openff.toolkit import Molecule
    from backend.app.docking.md_receptor_prep import _build_cofactor_molecule

    frame0 = md_dir / "frames" / "frame_000.pdb"
    if not frame0.is_file():
        raise FileNotFoundError(
            "MM-GBSA topology source frame_000.pdb missing; cannot build "
            f"complex modeller (looked at {frame0})"
        )
    pdb = app.PDBFile(str(frame0))
    modeller = app.Modeller(pdb.topology, pdb.positions)

    # small_mols feed SystemGenerator's GAFFTemplateGenerator. It maps each
    # HETATM residue in the topology to one of these Molecule objects via
    # graph isomorphism, so as long as the SMILES correspond to GTP / GDP /
    # ligand the chain-ID layout above doesn't matter.
    meta = json.loads((docking_dir / "metadata.json").read_text())
    cofactor_info = (meta.get("receptor_prep") or {}).get("parameterized_cofactors") or []
    small_mols = []
    cofactor_smiles: List[str] = []
    for cof in cofactor_info:
        smi = cof.get("smiles")
        if not smi:
            continue   # MG / ZN — handled by amber14 ion templates
        small_mols.append(Molecule.from_smiles(smi, allow_undefined_stereo=True))
        cofactor_smiles.append(smi)

    ligand_smiles = read_ligand_smiles(docking_job_id, ligand_name)
    if ligand_smiles is None:
        raise RuntimeError(
            f"no ligand SMILES for {ligand_name} in docking job {docking_job_id}"
        )
    # Use _build_cofactor_molecule (PDB + SMILES → stereo-correct Molecule)
    # so the GAFF cache key matches what engine_openmm wrote at MD time.
    ligand_pdb = md_dir / "_openmm_input_ligand.pdb"
    ligand_mol = _build_cofactor_molecule(
        ligand_pdb.read_text(), ligand_smiles, "docked_ligand",
    )
    small_mols.append(ligand_mol)

    return modeller, small_mols, "LIG", ligand_smiles, cofactor_smiles


def _make_system_generator(small_mols, cache_path: Optional[str] = None):
    """Build one SystemGenerator shared across the complex / receptor / ligand
    subsystems. Same parameterization knobs as engine_openmm.py — keep in sync.

    Why one shared SystemGenerator instead of three: the per-docking-job GAFF
    template cache written by the MD engine is DELTA-ENCODED — entry 1 (e.g.
    GTP) defines the full GAFF atom-type set; later entries (GDP, taxol, …)
    only define types that hadn't been registered yet. A fresh ligand-only
    SystemGenerator hits a ForceField with no prior atom types loaded, so its
    cached delta-ffxml fails with KeyError: 'c3'. A single SystemGenerator
    threaded through complex → receptor → ligand accumulates all atom types
    during the complex build; the receptor + ligand calls then reuse the
    fully-populated ForceField.

    `cache_path` (optional) — points at the MD engine's per-docking-job
    `_gaff_template_cache.json`. Cache HIT avoids the AM1-BCC backend
    entirely (which is env-dependent and not always present); MISS falls
    through to AM1-BCC via the toolkit registry and will fail loudly if no
    backend is installed.
    """
    from openmm import unit, app
    from openmmforcefields.generators import SystemGenerator

    kwargs = dict(
        forcefields=_FF_PROTEIN_AND_IONS + [_FF_IMPLICIT],
        small_molecule_forcefield=_FF_SMALL_MOL,
        molecules=small_mols,
        forcefield_kwargs={"constraints": app.HBonds},
        nonperiodic_forcefield_kwargs={
            "nonbondedMethod": app.CutoffNonPeriodic,
            "nonbondedCutoff": 1.0 * unit.nanometer,
        },
    )
    if cache_path:
        kwargs["cache"] = cache_path
    return SystemGenerator(**kwargs)


def prepare_mmgbsa_systems(
    md_id: str,
    *,
    docking_job_id: Optional[str] = None,
    ligand_name: Optional[str] = None,
    jobs_dir: Optional[Path] = None,
    use_cache: bool = True,
):
    """Return (cpx_sys, rec_sys, lig_sys, cpx_top, rec_top, lig_top,
    receptor_idx, ligand_idx). Loads from `<md_dir>/_mmgbsa_cache/` if
    available; otherwise builds via SystemGenerator (slow, sqm-bound)
    and writes the cache so subsequent calls are instant.

    `docking_job_id` and `ligand_name` are read from
    `<md_dir>/summary.json` when omitted — matches the existing MD job
    conventions.
    """
    base = jobs_dir if jobs_dir is not None else JOBS_DIR
    md_dir = base / md_id / "md"

    summary = json.loads((md_dir / "summary.json").read_text())
    if docking_job_id is None:
        docking_job_id = summary["docking_job_id"]
    if ligand_name is None:
        ligand_name = summary["ligand"]
    docking_dir = base / docking_job_id

    # Resolve the cache key BEFORE the (potentially expensive) build, so
    # a cache hit short-circuits the whole pipeline.
    ligand_smiles = read_ligand_smiles(docking_job_id, ligand_name)
    cofactor_info = (json.loads((docking_dir / "metadata.json").read_text())
                     .get("receptor_prep") or {}).get("parameterized_cofactors") or []
    cofactor_smiles = [c["smiles"] for c in cofactor_info if c.get("smiles")]
    key = _cache_key(ligand_smiles or "", cofactor_smiles)

    cache = _cache_dir(md_dir)
    if use_cache:
        hit = _try_load_cache(cache, key)
        if hit is not None:
            logger.info("mmgbsa cache HIT for %s (key=%s)", md_id, key)
            return hit

    logger.info("mmgbsa cache MISS for %s (key=%s) — rebuilding systems", md_id, key)
    modeller, small_mols, lig_resname, _, _ = _build_complex_modeller(
        md_dir, docking_dir, docking_job_id, ligand_name,
    )
    recm, ligm, ridx, lidx = _split_modeller(modeller, ligand_residue_name=lig_resname)
    # Reuse the MD pipeline's persistent GAFF template cache (per docking
    # job) AND share ONE SystemGenerator across the three subsystem builds.
    # The cache is delta-encoded — atom types are only re-written when a
    # later molecule introduces a new type — so a fresh SystemGenerator
    # for the ligand-only subsystem hits a ForceField with no prior atom
    # types loaded and KeyError-explodes on the first delta-entry it loads
    # (e.g. 'c3'). One shared generator threads through complex → receptor
    # → ligand, accumulating atom types on the complex build.
    gaff_cache = docking_dir / "_gaff_template_cache.json"
    gaff_cache_path = str(gaff_cache) if gaff_cache.is_file() else None
    if gaff_cache_path:
        logger.info("reusing GAFF template cache from %s", gaff_cache_path)
    sg = _make_system_generator(small_mols, gaff_cache_path)
    cpx_sys = sg.create_system(modeller.topology)
    rec_sys = sg.create_system(recm.topology)
    lig_sys = sg.create_system(ligm.topology)

    if use_cache:
        _write_cache(
            cache, key,
            cpx_sys, rec_sys, lig_sys,
            modeller.topology, recm.topology, ligm.topology,
            ridx, lidx,
        )

    return (cpx_sys, rec_sys, lig_sys,
            modeller.topology, recm.topology, ligm.topology,
            ridx, lidx)


def _collect_frame_positions(
    md_dir: Path,
    *,
    expected_topology=None,
    solvent_mode: str = "implicit",
) -> List:
    """Yield positions (nm) for each frame in md/frames/frame_*.pdb.

    When the trajectory comes from an explicit-solvent MD run, frames may
    contain waters + neutralization ions that the MM-GBSA complex system
    does NOT contain (GBSA is always implicit; waters + ions are stripped
    at the GBSA-system construction stage by `_build_complex_modeller`).
    Feeding waterful positions into the GBSA context would either
    out-of-range crash (too many atoms) or — worse — silently bind the
    wrong atoms.

    Stripping rule, defense-in-depth in both modes:

      1. Build an atom-key mask from `expected_topology` —
         {(chain.id, residue.name, residue.id, atom.name) → expected index}.
         This is the GBSA complex system's atom set; any frame atom whose
         key isn't here is by definition not part of the GBSA system.

      2. For each frame, walk PDB ATOM/HETATM lines in file order; for
         each line, parse its key; if the key is in the mask, record its
         (expected_index, x, y, z). Atoms whose key isn't in the mask —
         waters (HOH/WAT/TIP3), neutralization ions (NA/CL/K added by
         addSolvent), and any stray HETATM that the GBSA build dropped —
         are skipped silently.

      3. Emit positions sorted by expected_index, producing a length-N
         array that aligns 1:1 with the GBSA topology's atom order.

      4. Assert the matched atom count equals `expected_topology.getNumAtoms()`.
         An undermatch means the trajectory is missing solute atoms the
         GBSA system needs (e.g. legacy implicit run where receptor.pdb
         had ZN A:900 that the post-[A1] GBSA build no longer expects, or
         a renamed ligand residue); we raise rather than guess.

    When `expected_topology` is None the function falls back to the
    legacy behavior (read positions as-is, no mask, no validation), so
    callers that don't have the topology handy can still use this — but
    every production call site should pass `expected_topology`. The
    `solvent_mode` argument is surfaced in error messages and used by
    callers to log the trigger.
    """
    import openmm.app as app
    import numpy as np
    from openmm import unit

    expected_keys: Optional[Dict[Tuple[str, str, str, str], int]] = None
    if expected_topology is not None:
        expected_keys = {}
        for atom in expected_topology.atoms():
            ch = (atom.residue.chain.id or " ")
            key = (ch, atom.residue.name, str(atom.residue.id), atom.name)
            expected_keys[key] = atom.index

    out: List = []
    for p in sorted(md_dir.glob("frames/frame_*.pdb")):
        pdb = app.PDBFile(str(p))
        all_positions_nm = np.asarray(
            pdb.positions.value_in_unit(unit.nanometer)
        )
        if expected_keys is None:
            out.append(all_positions_nm)
            continue

        # Walk the topology of the loaded PDB (NOT the raw ATOM lines —
        # OpenMM's parser already normalizes element/coords/residue IDs,
        # so the key construction matches expected_keys exactly).
        matched: List[Tuple[int, int]] = []  # (expected_index, src_index)
        for atom in pdb.topology.atoms():
            ch = (atom.residue.chain.id or " ")
            key = (ch, atom.residue.name, str(atom.residue.id), atom.name)
            tgt_idx = expected_keys.get(key)
            if tgt_idx is None:
                continue  # solvent / ion / unmasked HETATM — skip silently
            matched.append((tgt_idx, atom.index))

        n_expected = expected_topology.getNumAtoms()
        if len(matched) != n_expected:
            missing = n_expected - len({m[0] for m in matched})
            raise RuntimeError(
                f"MM-GBSA frame strip mismatch on {p.name} "
                f"(solvent_mode={solvent_mode}): expected {n_expected} "
                f"GBSA atoms, matched {len(matched)} ({missing} missing). "
                "Trajectory atoms (chain, resname, resseq, atom_name) don't "
                "align 1:1 with the GBSA complex topology. Likely cause: "
                "the MD receptor was prepped with a different artifact-strip "
                "rule or atom-name convention than the GBSA build sees now."
            )

        # Emit positions in GBSA-topology order.
        matched.sort(key=lambda m: m[0])
        filtered = np.empty((n_expected, 3), dtype=all_positions_nm.dtype)
        for tgt_idx, src_idx in matched:
            filtered[tgt_idx] = all_positions_nm[src_idx]
        out.append(filtered)
    return out


# ---------------------------------------------------------------------
# Top-level orchestrator — self-describing FE artifact
# ---------------------------------------------------------------------

def _git_sha(repo_root: Optional[Path] = None) -> Optional[str]:
    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[3]
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(repo_root), capture_output=True, text=True, check=False, timeout=5,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    except (FileNotFoundError, subprocess.SubprocessError, OSError):
        pass
    return None


def compute_md_fe(
    md_id: str,
    *,
    jobs_dir: Optional[Path] = None,
    use_cache: bool = True,
    equilibration_discard_frames: Optional[int] = None,
) -> Dict[str, Any]:
    """Run the full pipeline on an MD job. Returns a `free_energy` block
    suitable for merging into the MD `summary.json`.

    `equilibration_discard_frames` drops the first N frames as
    equilibration. None (default) → auto-derive from the MD summary's
    settings:

      - If `settings.equilibration_discard_ps` is set, use that.
      - Else if `settings.solvent == 'explicit'`, use the default
        20 ps restraint-release-transient window (production NPT starts
        the instant the position restraints are released; waters need
        a few ps to settle around the now-unrestrained solute, and
        those frames distort the GBSA energy if included).
      - Else 0 (implicit runs come out of the heat phase already
        equilibrated; no production-side discard needed).

    Pass 0 to force-disable the discard (e.g. in tests). The threshold
    gates (≥200 total / ≥50 sampled) intentionally count ALL trajectory
    frames, not post-discard — the gates are about trajectory length,
    not analyzability; the discard is about which subset feeds the
    estimator. This deliberately keeps the gates honest under explicit
    mode (a shorter explicit run does NOT auto-pass).
    """
    base = jobs_dir if jobs_dir is not None else JOBS_DIR
    md_dir = base / md_id / "md"

    t0 = time.perf_counter()

    # Read the MD summary's solvent mode to drive both the [B4] trajectory
    # strip AND the auto equilibration-discard window. An explicit-solvent MD
    # writes solute-only frames per [B2], so the strip is a no-op in the
    # happy path; but the topology-mask filter is the safety net that
    # catches a future variant (or a buggy implicit MD whose frame atom set
    # doesn't match the GBSA system).
    md_solvent_mode = "implicit"  # default for pre-[B1] summaries
    md_settings: Dict[str, Any] = {}
    try:
        md_summary = json.loads((md_dir / "summary.json").read_text())
        md_settings = md_summary.get("settings") or {}
        md_solvent_mode = md_settings.get("solvent", "implicit")
    except Exception:
        pass  # missing/unreadable summary is gated elsewhere

    # Auto-derive the discard window from settings when the caller didn't
    # supply one. `equilibration_discard_frames=0` is a deliberate override
    # (skip discard); only `None` triggers auto-sizing.
    discard_source = "caller"
    if equilibration_discard_frames is None:
        snap_ps = float(md_settings.get("snapshot_every_ps") or 5.0)
        settings_discard_ps = md_settings.get("equilibration_discard_ps")
        if settings_discard_ps is not None:
            discard_ps = float(settings_discard_ps)
            discard_source = "settings.equilibration_discard_ps"
        elif md_solvent_mode == "explicit":
            discard_ps = 20.0
            discard_source = "auto: 20 ps for explicit restraint-release transient"
        else:
            discard_ps = 0.0
            discard_source = "auto: 0 ps for implicit (heat is pre-snapshot)"
        equilibration_discard_frames = max(0, int(round(discard_ps / snap_ps)))
    else:
        # Caller-supplied integer — derive a ps value for provenance.
        snap_ps = float(md_settings.get("snapshot_every_ps") or 5.0)
        discard_ps = equilibration_discard_frames * snap_ps

    # Live gate decision — same call the dashboard + report bundle use.
    gate = validate(md_id, jobs_dir=base).to_dict()

    # Systems (cache-aware).
    (cpx_sys, rec_sys, lig_sys,
     cpx_top, rec_top, lig_top,
     ridx, lidx) = prepare_mmgbsa_systems(
        md_id, jobs_dir=base, use_cache=use_cache,
    )

    # Frames — masked against the GBSA complex topology so explicit-solvent
    # trajectories have their waters + neutralization ions filtered out and
    # any atom-ordering mismatch raises loudly instead of corrupting energies.
    positions = _collect_frame_positions(
        md_dir,
        expected_topology=cpx_top,
        solvent_mode=md_solvent_mode,
    )
    if equilibration_discard_frames > 0:
        positions = positions[equilibration_discard_frames:]
    # Last-30 % equilibrated window — matches the corrected-gate
    # convention in docs/all_md_corrected_gate.md, which qualifies a
    # run on the last 30 % of frames. MM-GBSA was previously
    # integrating over ALL post-discard frames, including the
    # unequilibrated head; the ΔG would then average in the early
    # portion the gate explicitly excludes.
    n_post_discard = len(positions)
    eq_window_start = int(n_post_discard * 0.7)
    n_pre_window = eq_window_start
    positions = positions[eq_window_start:]
    if not positions:
        raise RuntimeError(
            f"no frames to evaluate (after discarding {equilibration_discard_frames} "
            f"head + skipping {n_pre_window} pre-equilibrated)"
        )

    method_meta = {
        "name": "single-trajectory MM-GBSA",
        "implicit_solvent": "OBC2",
        "force_fields": _FF_PROTEIN_AND_IONS + [_FF_IMPLICIT],
        "small_molecule_forcefield": _FF_SMALL_MOL,
        "nonbonded_method": "CutoffNonPeriodic",
        "nonbonded_cutoff_nm": 1.0,
        "configurational_entropy": "omitted (normal-mode out of scope)",
        "components_basis": "force-group dispatch (bonded / nonbonded / solvation)",
        "estimator": "backend.app.free_energy.mmgbsa:estimate_mmgbsa",
        # MD-trajectory solvent context. The GBSA energy eval itself is
        # always implicit OBC2 (above); this records the SOURCE trajectory's
        # solvent mode so the report bundle can disclose it. "explicit"
        # trajectories have their waters + ions filtered at frame load.
        "md_trajectory_solvent": md_solvent_mode,
    }

    result: MMGBSAResult = estimate_mmgbsa(
        complex_system=cpx_sys, complex_topology=cpx_top,
        receptor_system=rec_sys, receptor_topology=rec_top,
        ligand_system=lig_sys, ligand_topology=lig_top,
        frame_positions=positions,
        receptor_idx=ridx, ligand_idx=lidx,
        gate_can_run=gate["can_run"],
        gate_reason=(gate["reasons"][0]["message"] if gate["reasons"] else None),
        method_meta=method_meta,
    )

    wall_s = time.perf_counter() - t0

    # Self-describing free_energy block.
    return {
        "status": "completed",
        "schema_version": FE_ARTIFACT_SCHEMA_VERSION,
        "method": method_meta,
        "criteria": {
            "gate": {
                "can_run": gate["can_run"],
                "reasons": gate.get("reasons") or [],
                "warnings": gate.get("warnings") or [],
                "thresholds": gate.get("thresholds") or {},
            },
        },
        "provenance": {
            "git_sha": _git_sha(),
            "computed_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_frames_used": result.n_frames_used,
            "n_frames_skipped": result.n_frames_skipped,
            "equilibration_discard_frames": equilibration_discard_frames,
            "equilibration_discard_ps": discard_ps,
            "equilibration_discard_source": discard_source,
            "eq_window_last_fraction": 0.30,
            "eq_window_skipped_frames": n_pre_window,
            "md_trajectory_solvent": md_solvent_mode,
            "wall_seconds": round(wall_s, 2),
            "estimator_module": "backend.app.free_energy.mmgbsa_runner",
        },
        "result": {
            "sampling_adequate": result.sampling_adequate,
            "preliminary": result.preliminary,
            "gate_reason": result.gate_reason,
            "delta_g_mean_kcal_per_mol":   result.delta_g_mean,
            "delta_g_sem_kcal_per_mol":    result.delta_g_sem,
            "delta_g_stddev_kcal_per_mol": result.delta_g_stddev,
            "components_mean_kcal_per_mol": result.components_mean,
            "components_sem_kcal_per_mol":  result.components_sem,
        },
        "per_frame": [fe.to_dict() for fe in result.per_frame],
    }


def write_md_fe(md_id: str, *, jobs_dir: Optional[Path] = None,
                use_cache: bool = True,
                equilibration_discard_frames: Optional[int] = None,
                ) -> Dict[str, Any]:
    """Run compute_md_fe and persist the result into the MD summary.json's
    free_energy block, flipping status from `planned` → `completed`.

    Backs up the prior summary as `summary.json.pre-fe` once so the
    original `free_energy.status: "planned"` row is recoverable.
    """
    base = jobs_dir if jobs_dir is not None else JOBS_DIR
    md_dir = base / md_id / "md"
    summary_path = md_dir / "summary.json"

    fe_block = compute_md_fe(md_id, jobs_dir=base, use_cache=use_cache,
                             equilibration_discard_frames=equilibration_discard_frames)

    summary = json.loads(summary_path.read_text())
    backup = md_dir / "summary.json.pre-fe"
    if not backup.exists():
        backup.write_text(json.dumps(summary, indent=2))

    summary["free_energy"] = fe_block
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return fe_block
