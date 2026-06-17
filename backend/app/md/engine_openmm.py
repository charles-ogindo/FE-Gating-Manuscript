"""
OpenMM engine — preferred MD backend.

Lazy-imports openmm so that importing this module is always cheap. Returns
either a list of snapshot PDB paths + log path, or an `OpenMMNotAvailable`
exception that the orchestrator catches and falls back from.

Force-field stack: amber14-all + amber14/tip3p (for ion templates) + implicit
OBC2 solvent, with GAFF-2.11 (via openmmforcefields SystemGenerator) for the
docked ligand and any parameterized cofactors. openmmforcefields + openff
are hard runtime requirements.

Checkpoint / restart
--------------------
Every successful run serializes its final simulation State (positions +
velocities + periodic-box vectors + time) to ``<out_dir>/final_state.xml`` via
``Simulation.saveState`` — the portable XML form, not the binary,
platform-locked ``saveCheckpoint`` — plus a ``final_state_meta.json`` sidecar
recording the run conditions (force-field stack, temperature, timestep,
friction, snapshot interval, final frame index, timestamp).

Passing ``restart_from_checkpoint=<parent>/final_state.xml`` resumes a run:
the engine rebuilds the *identical* System/topology from the same inputs,
``loadState``s the parent's positions AND velocities, and continues production
with **no re-minimization, no re-heating, and no velocity reseed** — so the
extended leg is thermodynamically continuous with the parent's final step
(unlike a fresh run, which Maxwell-Boltzmann-seeds velocities and would
introduce a re-equilibration discontinuity). Compatibility (temperature,
timestep, friction, snapshot interval, force-field stack) is validated against
the checkpoint metadata and a mismatch raises rather than silently drifting.
New frames continue the parent's numbering (parent ended at frame_{offset};
the restart writes frame_{offset+1}..), so the two legs concatenate into one
logical trajectory. A continuous restart is only possible for runs produced
*after* this feature landed — older runs (e.g. 1f01da83) saved coordinate-only
PDB frames with no velocities and cannot be continuously extended.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# Checkpoint format (bump when the meta schema changes).
CHECKPOINT_SCHEMA_VERSION = "1.1.0"  # 1.1.0 adds explicit-solvent fields
CHECKPOINT_STATE_FILENAME = "final_state.xml"
CHECKPOINT_META_FILENAME = "final_state_meta.json"
# Explicit-mode only — the System and the solvated topology cannot be
# rebuilt deterministically on restart (addSolvent generates fresh waters
# each call). [B3] adds the load path; [B2] writes the files for any
# future restart-ready run.
CHECKPOINT_SYSTEM_FILENAME = "system.xml"
CHECKPOINT_TOPOLOGY_FILENAME = "topology.pdb"

# Force-field stack — single source of truth shared by create_system and the
# restart-compatibility validator so the two can never silently diverge.
# Implicit (OBC2) stack — kept verbatim for the default path.
_FORCE_FIELDS: Tuple[str, ...] = (
    "amber14-all.xml",
    "amber14/tip3p.xml",
    "implicit/obc2.xml",
)
_SMALL_MOLECULE_FF = "gaff-2.11"


def _explicit_forcefields(water_model: str) -> Tuple[str, ...]:
    """Force-field XML stack for explicit-solvent runs.

    amber14-all for protein/ions, amber14/<water_model>.xml for the water
    template. The implicit/obc2 term is dropped (PME handles solvation
    electrostatics; OBC2 must not be summed on top).
    """
    return ("amber14-all.xml", f"amber14/{water_model}.xml")


class OpenMMNotAvailable(RuntimeError):
    pass


@dataclass
class OpenMMResult:
    snapshots: List[Path]
    times_ps: List[float]
    log_path: Path
    engine_kind: str  # always "openmm_full"
    # Restart support (populated on every successful run; restart_info is None
    # for fresh runs and a provenance dict for continuations).
    checkpoint_path: Optional[Path] = None
    checkpoint_meta: Optional[dict] = None
    restart_info: Optional[dict] = None
    # Engine-side artifact-metal strip: provenance for any crystallographic
    # packing ions the engine removed from the receptor at load time. None
    # when the receptor was already filtered at receptor-prep time (post
    # [A1]) or when no monoatomic metals were present.
    receptor_strip: Optional[List[dict]] = None
    # Tail-truncation provenance. None when no MdRunSettings.truncate_chain_ranges
    # was configured. Each entry: {chain, from_resseq, to_resseq, n_residues,
    # cap, cap_residue_id, anchor_residue_id, ...}.
    receptor_truncation: Optional[List[dict]] = None


def run_openmm(
    receptor_pdb: Path,
    ligand_pdb: Path,
    out_dir: Path,
    *,
    ligand_smiles: str,
    cofactor_info: Optional[List[dict]] = None,
    docking_job_dir: Optional[Path] = None,
    production_ps: float = 500.0,
    timestep_fs: float = 2.0,
    snapshot_every_ps: float = 5.0,
    temperature_k: float = 300.0,
    friction_per_ps: float = 1.0,
    minimize_max_iter: int = 200,
    restart_from_checkpoint: Optional[Path] = None,
    parent_md_job_id: Optional[str] = None,
    frame_offset_hint: Optional[int] = None,
    # --- Explicit-solvent mode (additive; defaults preserve implicit path) ---
    solvent: str = "implicit",
    water_model: str = "tip3p",
    water_padding_nm: float = 1.0,
    ionic_strength_molar: float = 0.15,
    pressure_bar: float = 1.0,
    barostat_frequency_steps: int = 25,
    npt_equilibration_ps: float = 100.0,
    position_restraint_k_kj_per_mol_per_nm2: float = 1000.0,
    box_shape: str = "cube",
    truncate_chain_ranges: Optional[List[Dict[str, Any]]] = None,
    # --- Reproducible-replicate plumbing ---
    # `random_seed` applies to BOTH the Langevin integrator's stochastic
    # noise (`integrator.setRandomNumberSeed`) AND the initial Maxwell-
    # Boltzmann velocity draw (`setVelocitiesToTemperature(T, seed)`).
    # When None (default) OpenMM's clock-based defaults are used and the
    # run is non-deterministic; when set, the trajectory is reproducible
    # given the same starting structure + force-field stack. The seed
    # is also persisted into the engine result so the replicate that
    # produced a given trajectory is recoverable from disk.
    random_seed: Optional[int] = None,
) -> OpenMMResult:
    """Run a short implicit-solvent MD of the receptor + ligand complex.

    Raises OpenMMNotAvailable if openmm is not importable.

    Fresh run (``restart_from_checkpoint=None``): minimize → 50 ps NVT heat
    (Maxwell-Boltzmann velocities) → production; frame_000 is the post-heating
    t=0 snapshot.

    Restart (``restart_from_checkpoint`` set): rebuild the identical System
    from the same receptor/ligand inputs, ``loadState`` the parent's positions
    + velocities + box (NO minimize/heat/reseed), and continue production. The
    parent's serialized State must sit next to its ``final_state_meta.json``
    sidecar; run conditions are validated against it and a mismatch raises
    ``ValueError``. New frames are numbered from ``final_frame_index + 1`` so
    they concatenate onto the parent trajectory. ``frame_offset_hint``, when
    given, is cross-checked against the checkpoint's recorded final frame index.

    Either way, the run ends by serializing a fresh checkpoint (State + meta)
    into ``out_dir`` so this run can itself be extended later.
    """
    try:
        # Lazy import — this raises ImportError if openmm wasn't installed.
        import openmm as mm
        from openmm import app, unit
    except Exception as e:
        raise OpenMMNotAvailable(
            "openmm is not importable in the active Python environment "
            "(install with `pip install openmm`)"
        ) from e

    # Solvent-mode dispatch. The explicit branch lands in [B2]; until then
    # this engine implements only the implicit-OBC2 path. A caller asking
    # for "explicit" should fail loudly rather than silently get implicit
    # results — that's exactly the kind of mode-mismatch the orchestrator
    # cannot recover from.
    if solvent not in ("implicit", "explicit"):
        raise ValueError(
            f"solvent must be 'implicit' or 'explicit', got {solvent!r}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = out_dir / "frames"
    frames_dir.mkdir(exist_ok=True)
    log_path = out_dir / "log.txt"
    log = log_path.open("w", encoding="utf-8")

    def say(msg: str) -> None:
        log.write(msg + "\n")
        log.flush()
        logger.info("[md/openmm] %s", msg)

    try:
        say(f"OpenMM version {mm.__version__}")
        say(f"Loading receptor {receptor_pdb.name} and ligand {ligand_pdb.name}")

        # Vina pose PDBs arrive as multi-model files whose obabel-written
        # short-form `MODEL n` records crash OpenMM's strict PDB parser
        # (int('') on the model-serial column). Sanitize the ligand pose to
        # a single clean model; receptor is already MD-grade from
        # md_receptor_prep (cofactors live in sidecar PDBs, H is built in).
        clean_ligand = out_dir / "_openmm_input_ligand.pdb"
        _write_openmm_input_pdb(ligand_pdb, clean_ligand, uniquify_atom_names=True)

        # Build modeller: protein-with-ions receptor + cofactor sidecars +
        # docked ligand. SystemGenerator parameterizes peptides via amber14,
        # ions via amber14/tip3p ion templates, cofactors + ligand via GAFF.
        # SystemGenerator's molecules list is graph-only — Molecule.from_smiles
        # is sufficient for template registration; coords come from the
        # modeller (sidecar PDBs already have H from md_receptor_prep).
        from openmmforcefields.generators import SystemGenerator
        from openff.toolkit import Molecule
        from backend.app.docking.md_receptor_prep import _build_cofactor_molecule

        receptor_struct = app.PDBFile(str(receptor_pdb))
        modeller = app.Modeller(receptor_struct.topology, receptor_struct.positions)

        # Engine-side artifact-metal strip. Safety net for receptors that
        # pre-date the receptor-prep strip ([A1]). Idempotent: returns []
        # when nothing is flagged. Runs BEFORE cofactor add so the classifier
        # sees the structure exactly as it landed from receptor.pdb (and so
        # cofactors don't accidentally count as their own near-neighbors).
        receptor_strip = _strip_artifact_metals_from_modeller(
            modeller,
            docking_job_dir=docking_job_dir,
            cofactor_info=cofactor_info,
            keep_metals=None,
            say=say,
        )

        # Tail truncation (e.g., remove the α-tubulin C-terminal E-hook so
        # the explicit-solvent box doesn't have to encase residues that
        # bloat the water shell without affecting the binding site).
        # Mutates the modeller in place; adds NME amide caps at each new
        # C-terminus. Empty default → no-op. After truncation, the modeller's
        # NME residues lack hydrogens (only N + CH3 placed geometrically);
        # we'll add them via Modeller.addHydrogens against amber14 below
        # — but ONLY when a truncation actually occurred. Skipping the
        # addHydrogens call when no truncation ran preserves the implicit
        # path's atom counts bit-for-bit (md_receptor_prep already added
        # hydrogens at receptor-prep time; rerunning would be a no-op for
        # existing residues but would still emit a log line).
        receptor_truncation = _truncate_chain_residues_with_oxt_cap(
            modeller,
            truncate_chain_ranges or [],
            mm=mm, app=app, unit=unit, say=say,
        )
        if receptor_truncation:
            # After OXT insertion + PDB-roundtrip, the new C-terminal residue
            # matches amber14's C-terminal residue template (CASP/CGLN/...
            # all carry OXT). addHydrogens fills in any C-terminal-specific
            # H atoms. ForceField needs protein (amber14-all.xml) + ion
            # templates (amber14/tip3p.xml) since md_receptor_prep retained
            # the Mg ion in the receptor.
            ff_h = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
            say(f"Adding hydrogens to truncated termini + newly-exposed atoms "
                f"({len(receptor_truncation)} truncation(s))")
            modeller.addHydrogens(ff_h)

        small_mols = []
        for cof in (cofactor_info or []):
            sidecar_rel = cof.get("sidecar_pdb")
            if not sidecar_rel:
                continue  # metal ion — already in receptor PDB (or stripped above)
            sidecar = (docking_job_dir / sidecar_rel) if docking_job_dir else Path(sidecar_rel)
            if not sidecar.exists():
                say(f"cofactor sidecar missing: {sidecar} — skipping")
                continue
            pdbf = app.PDBFile(str(sidecar))
            modeller.add(pdbf.topology, pdbf.positions)
            mol = Molecule.from_smiles(cof["smiles"], allow_undefined_stereo=True)
            small_mols.append(mol)
            say(f"cofactor added: {cof['residue_name']} {cof['chain']}{cof['residue_number']} ({mol.n_atoms} atoms)")

        # Docked ligand: build Molecule with H from SMILES + cleaned pose,
        # write a normalized H-complete sidecar so PDBFile gives one residue.
        ligand_mol = _build_cofactor_molecule(
            clean_ligand.read_text(), ligand_smiles, "docked_ligand",
        )
        ligand_h_sidecar = out_dir / "_openmm_ligand_with_h.pdb"
        ligand_mol.to_file(str(ligand_h_sidecar), file_format="PDB")
        # Uniquify atom names: openff's PDB writer assigns inline polar Hs (H1..H4)
        # AND a separate bulk-AddHs sequence (H1..H47), causing duplicate names
        # within the LIG residue that OpenMM's PDBFile parser silently drops.
        _lig_norm = []
        _lig_counts: dict = {}
        for _ln in ligand_h_sidecar.read_text().splitlines():
            if _ln.startswith(("ATOM  ", "HETATM")):
                _ln = "HETATM" + _ln[6:]
                _ln = _rename_atom(_ln, _lig_counts)
                _ln = _ln[:17] + "LIG" + _ln[20:]
                _ln = _ln[:21] + "X" + _ln[22:]
                _ln = _ln[:22] + "   1" + _ln[26:]
            _lig_norm.append(_ln)
        ligand_h_sidecar.write_text("\n".join(_lig_norm) + "\n")
        ligand_pdbf = app.PDBFile(str(ligand_h_sidecar))
        modeller.add(ligand_pdbf.topology, ligand_pdbf.positions)
        small_mols.append(ligand_mol)
        say(f"docked ligand added: {ligand_mol.n_atoms} atoms")

        engine_kind = "openmm_full"

        # ---------- Restart meta pre-load + validation ----------
        # Read the checkpoint's metadata upfront so the System-build dispatch
        # below can take the LOAD path for an explicit restart (which must
        # NOT re-run addSolvent — the water shell would be different and the
        # continuation would not be thermodynamically continuous).
        ckpt_meta_for_restart: Optional[dict] = None
        if restart_from_checkpoint is not None:
            ckpt = Path(restart_from_checkpoint)
            if not ckpt.exists():
                raise ValueError(f"restart checkpoint not found: {ckpt}")
            ckpt_meta_for_restart = _load_checkpoint_meta(ckpt)
            _validate_restart_compatibility(
                ckpt_meta_for_restart,
                temperature_k=temperature_k,
                timestep_fs=timestep_fs,
                friction_per_ps=friction_per_ps,
                snapshot_every_ps=snapshot_every_ps,
                solvent=solvent,
                water_model=water_model,
                water_padding_nm=water_padding_nm,
                ionic_strength_molar=ionic_strength_molar,
                pressure_bar=pressure_bar,
                barostat_frequency_steps=barostat_frequency_steps,
                box_shape=box_shape,
                say=say,
            )

        # ---------- System-build dispatch ----------
        # Four mini-branches: (implicit|explicit) × (fresh|restart). Implicit
        # fresh == implicit restart at the build level (same SystemGenerator);
        # explicit fresh runs addSolvent + builds restraints; explicit restart
        # loads the parent's serialized system + topology + solute template.
        # All four converge on `system`, `sim_topology`, `sim_positions`,
        # `solute_topology`, `n_solute_atoms`, with sim_topology being what
        # Simulation gets bound to.
        solute_topology = None
        n_solute_atoms = 0
        sim_topology = None
        sim_positions = None

        if solvent == "explicit" and restart_from_checkpoint is not None:
            # ---- explicit RESTART: load saved System + topology ----
            # ckpt was already validated above; load the three companion files.
            sys_xml_path = ckpt.parent / CHECKPOINT_SYSTEM_FILENAME
            top_pdb_path = ckpt.parent / CHECKPOINT_TOPOLOGY_FILENAME
            solute_template_path = ckpt.parent / "_explicit_solute_template.pdb"
            for p, label in (
                (sys_xml_path, "system.xml"),
                (top_pdb_path, "topology.pdb"),
                (solute_template_path, "_explicit_solute_template.pdb"),
            ):
                if not p.exists():
                    raise FileNotFoundError(
                        f"explicit checkpoint missing {label} (looked at {p}); "
                        "only runs produced after [B2] can be continuously extended."
                    )
            say(f"Restart-explicit: loading saved system + topology + solute "
                f"template from {ckpt.parent}")
            system = mm.XmlSerializer.deserialize(
                sys_xml_path.read_text(encoding="utf-8")
            )
            top_pdb = app.PDBFile(str(top_pdb_path))
            sim_topology = top_pdb.topology
            sim_positions = top_pdb.positions  # unused — loadState supplies positions
            solute_template = app.PDBFile(str(solute_template_path))
            solute_topology = solute_template.topology
            meta_n_solute = ckpt_meta_for_restart.get("n_solute_atoms")
            n_solute_atoms = (int(meta_n_solute) if meta_n_solute is not None
                              else solute_topology.getNumAtoms())
            if n_solute_atoms != solute_topology.getNumAtoms():
                raise ValueError(
                    f"checkpoint meta.n_solute_atoms={meta_n_solute} "
                    f"!= solute template atom count "
                    f"{solute_topology.getNumAtoms()}"
                )
            say(f"Restart-explicit: n_solute_atoms={n_solute_atoms}; "
                f"restraint force lives in serialized system at k=0 "
                f"(no re-add)")

        elif solvent == "explicit":
            # ---- explicit FRESH: capture solute template, addSolvent, build ----
            solute_template_pdb = out_dir / "_explicit_solute_template.pdb"
            with solute_template_pdb.open("w", encoding="utf-8") as _f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, _f)
            solute_template = app.PDBFile(str(solute_template_pdb))
            solute_topology = solute_template.topology
            n_solute_atoms = solute_topology.getNumAtoms()

            # Build the SystemGenerator BEFORE addSolvent. addSolvent uses the
            # ForceField it's given to identify the solute (so it can place
            # waters around it); a forcefield without GAFF templates for the
            # cofactors + ligand fails with "No template found for residue
            # NNN (GTP)" the moment it walks into a small-molecule residue.
            #
            # SystemGenerator's underlying `sg.forcefield` IS the right
            # forcefield — it carries amber14 protein/ion templates AND the
            # GAFF small-molecule template generator. But template
            # registration is lazy: SystemGenerator only registers a small
            # molecule's GAFF template the first time forcefield is asked to
            # build a system for that molecule. addSolvent's forcefield
            # walk happens BEFORE the registration trigger, so we need to
            # force-trigger it ourselves with a throwaway create_system on
            # the pre-solvation modeller. openmmforcefields caches the sqm
            # AM1-BCC parameterization internally, so the post-solvation
            # create_system below pays no extra sqm cost.
            # Persistent GAFF cache. openmmforcefields' GAFFTemplateGenerator
            # accepts a `cache=<path>` argument and writes a JSON registry of
            # parameterized small-molecule templates keyed by molecular hash
            # (effectively SMILES). The first MD run against a docking job
            # pays the full sqm AM1-BCC cost (~25-90 min for taxol-sized
            # ligands); every subsequent MD run of the SAME docking_job_id
            # (same ligand + same cofactors) hits the cache and skips sqm
            # entirely — dropping MD wall from ~3 h to <5 min for the
            # non-sqm phases.
            sg_cache = (str(docking_job_dir / "_gaff_template_cache.json")
                        if docking_job_dir else None)
            if sg_cache:
                say(f"GAFF template cache: {sg_cache}")
            sg = SystemGenerator(
                forcefields=list(_explicit_forcefields(water_model)),
                small_molecule_forcefield=_SMALL_MOLECULE_FF,
                molecules=small_mols,
                forcefield_kwargs={
                    "constraints": app.HBonds,
                    "rigidWater": True,
                },
                periodic_forcefield_kwargs={
                    "nonbondedMethod": app.PME,
                    "nonbondedCutoff": 1.0 * unit.nanometer,
                },
                cache=sg_cache,
            )
            say(f"Pre-registering GAFF templates for {len(small_mols)} small "
                f"molecules (sqm AM1-BCC, cached on disk for reuse)")
            _ = sg.create_system(modeller.topology)

            n_atoms_before = modeller.topology.getNumAtoms()  # == n_solute_atoms
            say(f"Solvating: model={water_model}, padding={water_padding_nm} nm, "
                f"ionic_strength={ionic_strength_molar} M, neutralize=True, "
                f"boxShape={box_shape}")
            modeller.addSolvent(
                sg.forcefield,
                model=water_model,
                padding=water_padding_nm * unit.nanometer,
                ionicStrength=ionic_strength_molar * unit.molar,
                neutralize=True,
                boxShape=box_shape,
            )
            n_atoms_after = modeller.topology.getNumAtoms()
            say(f"Solvent added: {n_atoms_after - n_atoms_before} solvent atoms "
                f"({n_atoms_after} total, {n_solute_atoms} solute)")

            say(f"Building system: amber14 + {water_model} + PME + gaff-2.11 "
                f"({len(small_mols)} small mols, explicit)")
            system = sg.create_system(modeller.topology)

            # Position restraint on solute heavy atoms only (no H, no waters,
            # no neutralization ions). Released at production start by setting
            # the global parameter `k` to 0; the force stays in the System but
            # contributes 0 to energy and forces.
            restraint = mm.CustomExternalForce(
                "k*periodicdistance(x, y, z, x0, y0, z0)^2"
            )
            restraint.addGlobalParameter(
                "k", position_restraint_k_kj_per_mol_per_nm2
            )
            restraint.addPerParticleParameter("x0")
            restraint.addPerParticleParameter("y0")
            restraint.addPerParticleParameter("z0")
            n_restrained = 0
            for atom in modeller.topology.atoms():
                if atom.index >= n_solute_atoms:
                    break  # waters + neutralization ions appended after solute
                if atom.element is None or atom.element.symbol == "H":
                    continue
                p = modeller.positions[atom.index].value_in_unit(unit.nanometer)
                restraint.addParticle(atom.index, [p.x, p.y, p.z])
                n_restrained += 1
            system.addForce(restraint)
            say(f"Position restraints: k={position_restraint_k_kj_per_mol_per_nm2:.0f} "
                f"kJ/mol/nm² on {n_restrained} solute heavy atoms "
                f"(released to k=0 at production start)")
            sim_topology = modeller.topology
            sim_positions = modeller.positions
        else:
            # ---- implicit (fresh OR restart) ----
            sg_cache = (str(docking_job_dir / "_gaff_template_cache.json")
                        if docking_job_dir else None)
            sg = SystemGenerator(
                forcefields=list(_FORCE_FIELDS),
                small_molecule_forcefield=_SMALL_MOLECULE_FF,
                molecules=small_mols,
                forcefield_kwargs={"constraints": app.HBonds},
                nonperiodic_forcefield_kwargs={
                    "nonbondedMethod": app.CutoffNonPeriodic,
                    "nonbondedCutoff": 1.0 * unit.nanometer,
                },
                cache=sg_cache,
            )
            say(f"Building system: amber14 + tip3p ions + obc2 + gaff-2.11 "
                f"({len(small_mols)} small mols, implicit)")
            system = sg.create_system(modeller.topology)
            sim_topology = modeller.topology
            sim_positions = modeller.positions

        # ---------- Integrator + platform + Simulation ----------
        integrator = mm.LangevinMiddleIntegrator(
            temperature_k * unit.kelvin,
            friction_per_ps / unit.picosecond,
            timestep_fs * unit.femtosecond,
        )
        if random_seed is not None:
            # Reproducible Langevin noise — distinct seeds per replicate
            # decorrelate the stochastic-force draws, which (together with
            # the velocity-seed below) is what makes independent replicates
            # statistically independent rather than re-runs of the same
            # trajectory.
            integrator.setRandomNumberSeed(int(random_seed))
            say(f"Integrator seeded: {int(random_seed)}")
        platform = _pick_platform(mm, say)
        simulation = app.Simulation(sim_topology, system, integrator, platform)

        if restart_from_checkpoint is not None:
            # loadState restores positions + velocities + box vectors + time, so
            # the continuation is thermodynamically continuous with the parent's
            # final step. Deliberately NO minimize / NO heat / NO velocity reseed.
            frame_offset = int(ckpt_meta_for_restart.get("final_frame_index", 0))
            if (frame_offset_hint is not None
                    and int(frame_offset_hint) != frame_offset):
                raise ValueError(
                    f"frame_offset_hint ({frame_offset_hint}) disagrees with the "
                    f"checkpoint's final_frame_index ({frame_offset}); refusing "
                    "to guess trajectory numbering."
                )
            simulation.loadState(str(ckpt))
            say(f"Restart: loaded State {ckpt} (parent {parent_md_job_id or '?'}, "
                f"final_frame_index={frame_offset}); continuing with no "
                f"re-equilibration.")
        else:
            simulation.context.setPositions(sim_positions)
            # For explicit solvent, minimize needs to resolve water-solute
            # clashes that addSolvent may have produced (waters placed inside
            # van der Waals radii of solute atoms). 200 iters is fine for the
            # implicit OBC2 system (~13k atoms, no clashing waters); explicit
            # at ~50k atoms with PME needs an order of magnitude more, and
            # rigidWater + HBonds constraints sometimes still trip on initial
            # geometry — apply constraints first to clean those up.
            if solvent == "explicit":
                effective_minimize_iter = max(minimize_max_iter, 5000)
                simulation.context.applyConstraints(1e-6)
                say(f"Minimizing (explicit, max_iter={effective_minimize_iter}, "
                    f"constraints applied to initial positions)")
                simulation.minimizeEnergy(maxIterations=effective_minimize_iter)
            else:
                say(f"Minimizing (max_iter={minimize_max_iter})")
                simulation.minimizeEnergy(maxIterations=minimize_max_iter)
            heat_steps = int((50.0 * 1000.0) / timestep_fs)
            if solvent == "explicit":
                say(f"Heating: 50 ps restrained NVT @ {temperature_k} K "
                    f"({heat_steps} steps)")
            else:
                say(f"Heating: 50 ps NVT @ {temperature_k} K ({heat_steps} steps)")
            if random_seed is not None:
                # Distinct velocity draws per replicate — paired with the
                # integrator seed above, this is what makes independent
                # replicates statistically independent (PRD requirement (a)
                # + (b) for the replicate sweep).
                simulation.context.setVelocitiesToTemperature(
                    temperature_k * unit.kelvin, int(random_seed),
                )
            else:
                simulation.context.setVelocitiesToTemperature(
                    temperature_k * unit.kelvin,
                )
            simulation.step(heat_steps)

            if solvent == "explicit":
                # Switch ensemble NVT → NPT by adding MonteCarloBarostat to
                # the System and reinitializing the context (preserveState=True
                # keeps positions + velocities + box + time intact). Then
                # restrained NPT equilibration; then release restraints and
                # enter production NPT.
                barostat = mm.MonteCarloBarostat(
                    pressure_bar * unit.bar,
                    temperature_k * unit.kelvin,
                    barostat_frequency_steps,
                )
                system.addForce(barostat)
                simulation.context.reinitialize(preserveState=True)
                npt_steps = int((npt_equilibration_ps * 1000.0) / timestep_fs)
                say(f"NPT equilibration: {npt_equilibration_ps} ps restrained @ "
                    f"{pressure_bar} bar, {temperature_k} K ({npt_steps} steps, "
                    f"barostat every {barostat_frequency_steps})")
                simulation.step(npt_steps)
                # Release restraints — k=0 zeroes the CustomExternalForce
                # contribution without removing the force itself (cheap).
                simulation.context.setParameter("k", 0.0)
                say("Restraints released (k=0); entering production NPT")

            frame_offset = 0

        # Production with snapshots.
        prod_steps = int((production_ps * 1000.0) / timestep_fs)
        snap_steps = max(1, int((snapshot_every_ps * 1000.0) / timestep_fs))
        say(f"Production: {production_ps} ps ({prod_steps} steps, snapshot every "
            f"{snap_steps} steps = {snapshot_every_ps} ps)")

        def _write_frame(_path: Path) -> None:
            if solvent == "explicit":
                _dump_solute_pdb(
                    simulation, _path,
                    solute_topology=solute_topology,
                    n_solute_atoms=n_solute_atoms,
                )
            else:
                _dump_pdb(simulation, mm, unit, _path)

        snapshots: List[Path] = []
        times: List[float] = []
        if restart_from_checkpoint is None:
            # Frame 000 = post-equilibration starting point (this is "t=0"
            # for analysis). In explicit mode it is the moment AFTER
            # restraint release, so the analyzer's window starts in a
            # density-converged, unrestrained state.
            path0 = frames_dir / "frame_000.pdb"
            _write_frame(path0)
            snapshots.append(path0); times.append(0.0)
        else:
            # The loaded State == the parent's final frame (frame_{offset}),
            # which already exists in the parent run — do NOT re-dump it (that
            # would duplicate a frame). New frames continue at offset+1.
            say(f"Appending new frames starting at frame_{frame_offset + 1:03d}.")

        n_dumps = max(1, prod_steps // snap_steps)
        wall_start = time.perf_counter()
        for i in range(n_dumps):
            simulation.step(snap_steps)
            idx = frame_offset + i + 1
            path = frames_dir / f"frame_{idx:03d}.pdb"
            _write_frame(path)
            snapshots.append(path); times.append(idx * snapshot_every_ps)

        # Serialize a fresh checkpoint so THIS run can itself be extended later.
        # Explicit mode also writes system.xml + topology.pdb so [B3]'s restart
        # path can rebuild the identical solvated system without re-running
        # addSolvent (which would generate a different water shell each call).
        last_frame_index = frame_offset + n_dumps
        checkpoint_path, checkpoint_meta = _write_checkpoint(
            simulation, out_dir,
            engine_kind=engine_kind,
            final_frame_index=last_frame_index,
            n_frames_this_leg=len(snapshots),
            temperature_k=temperature_k,
            timestep_fs=timestep_fs,
            friction_per_ps=friction_per_ps,
            snapshot_every_ps=snapshot_every_ps,
            production_ps=production_ps,
            solvent=solvent,
            water_model=water_model,
            water_padding_nm=water_padding_nm,
            ionic_strength_molar=ionic_strength_molar,
            pressure_bar=pressure_bar,
            barostat_frequency_steps=barostat_frequency_steps,
            npt_equilibration_ps=npt_equilibration_ps,
            position_restraint_k_kj_per_mol_per_nm2=
                position_restraint_k_kj_per_mol_per_nm2,
            n_solute_atoms=n_solute_atoms,
            box_shape=box_shape,
            modeller_topology=modeller.topology,
            modeller_positions=modeller.positions,
            system=system,
            mm=mm,
            app=app,
            say=say,
        )

        restart_info: Optional[dict] = None
        if restart_from_checkpoint is not None:
            restart_info = {
                "is_restart": True,
                "parent_md_job_id": parent_md_job_id,
                "restart_from_checkpoint": str(restart_from_checkpoint),
                "restart_from_frame": frame_offset,
                "frame_offset": frame_offset,
                "first_new_frame_index": frame_offset + 1,
                "last_new_frame_index": last_frame_index,
                "frames_written_this_leg": len(snapshots),
            }

        say(f"Done. {len(snapshots)} new frames "
            f"(frame_{frame_offset + 1:03d}..frame_{last_frame_index:03d}) in "
            f"{time.perf_counter() - wall_start:.1f}s wall ({engine_kind}).")
        return OpenMMResult(
            snapshots=snapshots,
            times_ps=times,
            log_path=log_path,
            engine_kind=engine_kind,
            checkpoint_path=checkpoint_path,
            checkpoint_meta=checkpoint_meta,
            restart_info=restart_info,
            receptor_strip=receptor_strip or None,
            receptor_truncation=receptor_truncation or None,
        )
    finally:
        try:
            log.close()
        except Exception:
            pass


def _write_openmm_input_pdb(
    src: Path, dst: Path, *, uniquify_atom_names: bool = False
) -> None:
    """Write a cleaned copy of `src` that OpenMM's strict PDB parser accepts.

    Keeps only the first model's coordinate/bond records (ATOM/HETATM/TER/
    CONECT); drops MODEL/ENDMDL/REMARK/COMPND/etc. Vina pose PDBs arrive here
    as multi-model files whose obabel-written short-form `MODEL n` records
    crash OpenMM (`int('')` on the model-serial column).

    With `uniquify_atom_names=True` (ligand inputs) each ATOM/HETATM name is
    rewritten to `<element><running-index>` (C1, C2, ..., O1, ...). obabel
    names every ligand atom by bare element, which OpenMM collapses as
    duplicate / alt-loc atoms within the single ligand residue (66 atoms ->
    ~4). Must stay False for the receptor: OpenMM matches protein atoms to
    force-field residue templates BY NAME, so those names must be preserved.
    Atom serial numbers are never touched, so CONECT bond records stay valid.
    """
    kept: List[str] = []
    elem_counts: dict = {}
    for line in src.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("ENDMDL"):
            break  # first model only
        if line.startswith(("ATOM  ", "HETATM")):
            if uniquify_atom_names:
                line = _rename_atom(line, elem_counts)
                # Force HETATM — small-molecule ligand, not a protein residue;
                # RDKit's PDB parser treats ATOM as residue-template lookups.
                if line.startswith("ATOM  "):
                    line = "HETATM" + line[6:]
            kept.append(line)
        elif line.startswith(("TER", "CONECT")):
            kept.append(line)
    dst.write_text("\n".join(kept) + "\nEND\n", encoding="utf-8")


def _rename_atom(line: str, elem_counts: dict) -> str:
    """Rewrite an ATOM/HETATM line's name column (cols 13-16) to a
    residue-unique `<element><index>` AND overwrite cols 77-78 with the
    correct element symbol. Some Vina pose PDBs have `*` in the element
    column, which makes RDKit refuse the line."""
    raw_elem = line[76:78].strip() if len(line) >= 78 else ""
    if not raw_elem or not raw_elem.isalpha():
        raw_name = line[12:16].strip()
        raw_elem = "".join(c for c in raw_name if c.isalpha())[:2] or "X"
    element = raw_elem.capitalize()
    elem_counts[element] = elem_counts.get(element, 0) + 1
    name = f"{element}{elem_counts[element]}"
    new_line = line[:12] + f"{name:<4}"[:4] + line[16:]
    if len(new_line) < 78:
        new_line = new_line.ljust(78)
    return new_line[:76] + f"{element:>2}" + new_line[78:]


def _strip_artifact_metals_from_modeller(
    modeller,
    *,
    docking_job_dir: Optional[Path],
    cofactor_info: Optional[List[dict]],
    keep_metals: Optional[set],
    say,
) -> List[dict]:
    """Engine-side safety net for the crystallographic-artifact metal strip.

    Runs the same classifier as `prepare_md_grade_receptor` against the loaded
    modeller's monoatomic metal ions, and removes any flagged as
    crystallographic packing artifacts (no side-chain N/S/O donor within 3 Å
    and no organic-cofactor heavy atom within 4 Å). Mutates the modeller in
    place via `modeller.delete(residue_list)`. Returns the provenance list.

    Idempotent: receptor PDBs filtered at receptor-prep time have no
    artifact ions left and this returns [].

    Required only for docking jobs prepared BEFORE the receptor-prep strip
    landed; for new jobs it is a no-op. The dual-call-site design (prep +
    engine) guarantees consistent behavior across the existing and future
    job populations without re-prepping legacy artifacts.
    """
    from backend.app.docking.md_receptor_prep import (
        _METAL_IONS, _classify_metal_role,
    )
    from openmm import unit as mm_unit

    keep = set(keep_metals or ())

    # Collect cofactor heavy-atom positions (Angstrom) from the sidecar PDBs
    # the docking job emitted (GTP/GDP/ATP/...). The classifier needs these
    # to distinguish a Mg²⁺ chelating a phosphate (functional) from a Zn²⁺
    # perched against backbone (artifact).
    cofactor_heavy_xyz: List[Tuple[float, float, float]] = []
    for cof in (cofactor_info or []):
        sidecar_rel = cof.get("sidecar_pdb")
        if not sidecar_rel:
            continue  # metal-only entry; no organic to load
        sidecar = ((docking_job_dir / sidecar_rel)
                   if docking_job_dir else Path(sidecar_rel))
        if not sidecar.exists():
            continue
        for ln in sidecar.read_text(errors="ignore").splitlines():
            if not ln.startswith(("ATOM  ", "HETATM")) or len(ln) < 54:
                continue
            elem = ln[76:78].strip() if len(ln) >= 78 else ""
            if elem.upper() == "H":
                continue
            try:
                cofactor_heavy_xyz.append(
                    (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                )
            except ValueError:
                continue

    # Index metal residues in the topology + grab positions in Angstrom.
    metal_residues: List[Tuple[Any, Tuple[float, float, float], Tuple[str, str, str]]] = []
    pos = modeller.positions
    for residue in modeller.topology.residues():
        if residue.name.upper() not in _METAL_IONS:
            continue
        atoms = list(residue.atoms())
        if len(atoms) != 1:
            continue  # only monoatomic residues qualify
        atom = atoms[0]
        p = pos[atom.index].value_in_unit(mm_unit.angstrom)
        xyz = (p.x, p.y, p.z)
        chain_id = residue.chain.id or " "
        key = (chain_id, residue.name.upper(), str(residue.id))
        metal_residues.append((residue, xyz, key))

    if not metal_residues:
        return []

    metal_atom_indices = {list(r.atoms())[0].index for r, _, _ in metal_residues}

    # Build receptor heavy-atom records (Angstrom). Exclude metals so an ion
    # isn't counted as its own shell neighbor.
    receptor_atoms: List[Dict[str, Any]] = []
    for atom in modeller.topology.atoms():
        if atom.index in metal_atom_indices:
            continue
        elem = atom.element
        if elem is not None and elem.symbol == "H":
            continue
        p = pos[atom.index].value_in_unit(mm_unit.angstrom)
        receptor_atoms.append({
            "resname": atom.residue.name.upper(),
            "atom_name": atom.name,
            "x": p.x, "y": p.y, "z": p.z,
        })

    to_delete = []
    stripped: List[dict] = []
    for residue, xyz, key in metal_residues:
        if key in keep:
            continue
        classification = _classify_metal_role(
            xyz,
            receptor_atoms=receptor_atoms,
            cofactor_heavy_atoms_xyz=cofactor_heavy_xyz,
        )
        if classification["role"] == "artifact":
            ch, rn, rs = key
            stripped.append({
                "chain": ch,
                "residue_name": rn,
                "residue_number": rs,
                "xyz_input": [round(c, 3) for c in xyz],
                "rule_version": classification["rule_version"],
                "reason": classification["reason"],
                "shell": classification["shell"],
                "n_sidechain_donor_contacts":
                    classification["n_sidechain_donor_contacts"],
                "n_cofactor_heavy_atom_contacts":
                    classification["n_cofactor_heavy_atom_contacts"],
                "source": "engine_openmm",
            })
            to_delete.append(residue)
            say(f"Stripping crystallographic-artifact metal: "
                f"{rn} {ch}{rs} — {classification['reason']}")

    if to_delete:
        modeller.delete(to_delete)
    return stripped


def _truncate_chain_residues_with_oxt_cap(
    modeller, ranges: List[Dict[str, Any]], *, mm, app, unit, say,
) -> List[Dict[str, Any]]:
    """Delete a contiguous range of residues from a chain and terminate the
    new C-terminus naturally by inserting an OXT atom on the anchor
    residue. Mutates the modeller in place via a PDB-text round-trip
    (delete → write → text-insert OXT → reload).

    Each range is `{"chain": str, "from_resseq": int, "to_resseq": int}`.

    Cap choice — natural COO- (OXT) vs amide (NME):

      The original [B8] design used an NME amide cap, the textbook
      "honest" choice for an artificial truncation (it tells the
      simulation that this isn't the real protein end). In practice
      the OpenMM Modeller.addHydrogens path could not reliably complete
      a hand-built NME residue (its createSystem prelude fails on
      "NME missing N H atoms" even when the N + CH3 heavy atoms are
      named correctly per amber14 — the addHydrogens template walker
      partial-matches and short-circuits before adding the missing Hs).

      The OXT (natural C-terminus) approach is parameterized cleanly by
      amber14's C-terminal residue templates (CASP/CGLN/CARG/...), and
      the local charge perturbation from a free carboxylate is
      negligible at the >50 Å distance from the taxane pocket where
      the α-tubulin E-hook truncation occurs. The standard published
      tubulin-MD treatment of E-hook truncations uses OXT termini for
      exactly this reason.

      Provenance records `cap == "OXT"` for forward compatibility; an
      NME variant can be added later if a problem with hand-built NME
      residues + OpenMM addHydrogens is resolved upstream.

    OXT placement geometry (within the CA-C-O plane, ~120° from O across C):
      OXT = C - (O − C)   (collinear-opposite of O across C, |C-OXT| ≈ |C-O|)
    Minimize fixes the precise trigonal-planar angle before heat.

    Returns provenance list: one entry per range with the truncation
    extent, anchor residue, and cap atom info.
    """
    if not ranges:
        return []

    import io
    import numpy as np

    provenance: List[Dict[str, Any]] = []

    pos_nm = modeller.positions

    to_delete: List[Any] = []
    # Queue of OXT insertions: (chain_id, anchor_resseq, oxt_xyz_nm). After
    # the delete pass, we round-trip through PDB text to insert OXT lines.
    oxt_specs: List[Dict[str, Any]] = []

    # IMPORTANT: do not key chains by chain.id. md_receptor_prep appends
    # metal ions via a separate Topology with chain.id == "A" (same letter
    # as the protein chain A), so a flat dict overwrites the protein chain
    # with the ion chain. Walk chains explicitly per range and select the
    # one that holds the requested residue numbers.
    all_chains = list(modeller.topology.chains())

    for rng in ranges:
        ch_id = rng["chain"]
        a = int(rng["from_resseq"])
        b = int(rng["to_resseq"])
        if a > b:
            raise ValueError(
                f"truncation range invalid: from_resseq={a} > to_resseq={b}"
            )
        candidates = [c for c in all_chains if c.id == ch_id]
        if not candidates:
            present_ids = sorted({c.id for c in all_chains})
            raise ValueError(
                f"truncation chain {ch_id!r} not present in modeller "
                f"(have: {present_ids})"
            )
        # Pick the chain whose residue numbers cover the requested range.
        target_chain = None
        for c in candidates:
            resseqs = {int(r.id) for r in c.residues()}
            if any(rs in resseqs for rs in range(a, b + 1)):
                target_chain = c
                break
        if target_chain is None:
            # Surface ALL candidate ranges for the chain id so the error is
            # actionable when there are multiple chains with this id.
            descs = []
            for c in candidates:
                rs = sorted(int(r.id) for r in c.residues())
                if rs:
                    descs.append(f"{rs[0]}-{rs[-1]}")
                else:
                    descs.append("(empty)")
            raise ValueError(
                f"no residues in chain {ch_id} match resseq {a}-{b}; "
                f"present chain-{ch_id} ranges: {descs}"
            )
        residues = list(target_chain.residues())
        # Find the contiguous range in this chain.
        targets = [r for r in residues if int(r.id) in range(a, b + 1)]
        if not targets:
            raise ValueError(
                f"no residues in chain {ch_id} match resseq {a}-{b} "
                "(internal check; should have been caught above)"
            )
        # Anchor = the residue immediately PRECEDING the first deleted one
        # (this becomes the new C-terminus's residue 438 in α-tubulin).
        anchor = None
        for r in residues:
            if int(r.id) == a - 1:
                anchor = r
                break
        if anchor is None:
            raise ValueError(
                f"truncation requires an anchor residue at resseq {a - 1} "
                f"on chain {ch_id} (one position before from_resseq); not found"
            )
        # Pull the anchor's C / O / CA positions (nm), used to seed the
        # NME N + CH3 placements.
        anchor_atoms = {a_.name: a_ for a_ in anchor.atoms()}
        missing = [n for n in ("C", "O", "CA") if n not in anchor_atoms]
        if missing:
            raise ValueError(
                f"anchor residue {ch_id}{anchor.id} lacks backbone atoms "
                f"{missing}; cannot place NME cap"
            )
        p_c = np.asarray(
            pos_nm[anchor_atoms["C"].index].value_in_unit(unit.nanometer))
        p_o = np.asarray(
            pos_nm[anchor_atoms["O"].index].value_in_unit(unit.nanometer))
        p_ca = np.asarray(
            pos_nm[anchor_atoms["CA"].index].value_in_unit(unit.nanometer))

        # OXT placement: collinear-opposite of O across C (|C-OXT| ≈ |C-O|).
        # Minimize fixes the precise trigonal-planar angle vs the CA-C
        # axis before heat.
        v_co = p_o - p_c
        co_norm = float(np.linalg.norm(v_co))
        if co_norm < 1e-6:
            raise ValueError(
                f"degenerate anchor geometry at {ch_id}{anchor.id}; "
                "cannot place OXT cap"
            )
        p_oxt = p_c - v_co  # opposite-collinear at the same distance

        to_delete.extend(targets)
        oxt_specs.append({
            "chain_id": ch_id,
            "anchor_resseq": anchor.id,
            "anchor_resname": anchor.name,
            "oxt_xyz_nm": p_oxt.tolist(),
        })
        provenance.append({
            "chain": ch_id,
            "from_resseq": a,
            "to_resseq": b,
            "n_residues_removed": len(targets),
            "cap": "OXT",
            "cap_atom_name": "OXT",
            "anchor_residue_id": anchor.id,
            "anchor_residue_name": anchor.name,
            "anchor_C_xyz_nm": p_c.tolist(),
            "oxt_xyz_nm": p_oxt.tolist(),
            "source": "engine_openmm",
        })
        say(f"Truncating chain {ch_id} residues {a}-{b} "
            f"({len(targets)} residue(s)); adding OXT cap to "
            f"{ch_id}{anchor.id}({anchor.name}) → natural COO- C-terminus")

    # 1) Delete the truncated residues from the modeller.
    if to_delete:
        modeller.delete(to_delete)

    # 2) Round-trip through PDB text to insert OXT lines. Direct Topology
    # mutation (adding an atom to an existing residue while keeping
    # positions consistent) is brittle in OpenMM; PDB text manipulation
    # is reliable. We write the post-delete modeller to a buffer, insert
    # OXT lines at the right anchor residues, then reload as a fresh
    # modeller. The OpenMM PDBFile reader will assign element + bond
    # connectivity from amber14 residue templates when SystemGenerator
    # runs later — for the readback here we only need the atom record
    # parsed correctly.
    buf = io.StringIO()
    # keepIds=True preserves chain.id and residue.id during writeFile so
    # the cofactor_info dict's (chain, resseq) lookups stay valid through
    # the roundtrip. Without it, PDBFile renumbers chains sequentially
    # (A, B, ...) and residues 1.. — breaking the Mg ion's A:501 identity
    # that downstream metadata depends on.
    app.PDBFile.writeFile(modeller.topology, modeller.positions, buf, keepIds=True)
    pdb_text = buf.getvalue()

    # Index anchor (chain, resseq) for fast O(N) scan
    oxt_targets: Dict[Tuple[str, str], Tuple[float, float, float]] = {}
    for spec in oxt_specs:
        oxt_targets[(spec["chain_id"], str(spec["anchor_resseq"]))] = tuple(spec["oxt_xyz_nm"])

    out_lines: List[str] = []
    pending_inserts: Dict[Tuple[str, str], Tuple[float, float, float]] = dict(oxt_targets)
    # Scan PDB; after the LAST atom of an anchor residue, inject the OXT
    # ATOM line before any subsequent record. Track the current residue
    # key per scan; when it changes and the previous key was an anchor,
    # emit the OXT line right before the new record.
    prev_key: Optional[Tuple[str, str]] = None
    next_serial = 1
    for line in pdb_text.splitlines(keepends=False):
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21]
            resseq = line[22:26].strip()
            key = (chain, resseq)
            # Boundary: previous residue was an anchor and we're now in a
            # different residue → insert OXT for the previous anchor.
            if (prev_key is not None and prev_key in pending_inserts
                    and key != prev_key):
                ox, oy, oz = pending_inserts.pop(prev_key)
                # PDB cols: ATOM  NNNNN  OXT RES C RRRR    X.XXX   Y.YYY   Z.ZZZ
                ox_a = ox * 10.0; oy_a = oy * 10.0; oz_a = oz * 10.0
                # Use a placeholder serial; OpenMM PDBFile reader assigns
                # its own indices. Resname from the same anchor.
                # We need the anchor's resname — fetch from any atom of
                # the previous block in out_lines (anchor resname stays
                # constant for a residue).
                anchor_resname = ""
                for prev_line in reversed(out_lines):
                    if (prev_line.startswith(("ATOM  ", "HETATM"))
                            and prev_line[21] == prev_key[0]
                            and prev_line[22:26].strip() == prev_key[1]):
                        anchor_resname = prev_line[17:20]
                        break
                # Pad anchor_resname to 3 chars
                oxt_line = (
                    "ATOM  "
                    f"{next_serial:>5d}"
                    "  OXT "
                    f"{anchor_resname:>3s} "
                    f"{prev_key[0]}"
                    f"{int(prev_key[1]):>4d}    "
                    f"{ox_a:8.3f}{oy_a:8.3f}{oz_a:8.3f}"
                    "  1.00  0.00           O  "
                )
                out_lines.append(oxt_line)
                next_serial += 1
            prev_key = key
        elif line.startswith(("TER", "END", "ENDMDL", "MASTER")):
            # Same boundary logic at chain-end / file-end records.
            if prev_key is not None and prev_key in pending_inserts:
                ox, oy, oz = pending_inserts.pop(prev_key)
                ox_a = ox * 10.0; oy_a = oy * 10.0; oz_a = oz * 10.0
                anchor_resname = ""
                for prev_line in reversed(out_lines):
                    if (prev_line.startswith(("ATOM  ", "HETATM"))
                            and prev_line[21] == prev_key[0]
                            and prev_line[22:26].strip() == prev_key[1]):
                        anchor_resname = prev_line[17:20]
                        break
                oxt_line = (
                    "ATOM  "
                    f"{next_serial:>5d}"
                    "  OXT "
                    f"{anchor_resname:>3s} "
                    f"{prev_key[0]}"
                    f"{int(prev_key[1]):>4d}    "
                    f"{ox_a:8.3f}{oy_a:8.3f}{oz_a:8.3f}"
                    "  1.00  0.00           O  "
                )
                out_lines.append(oxt_line)
                next_serial += 1
            prev_key = None
        out_lines.append(line)

    if pending_inserts:
        raise RuntimeError(
            f"OXT insert failed: {len(pending_inserts)} anchor(s) not "
            f"matched in the post-delete PDB ({list(pending_inserts.keys())})"
        )

    # Reload as a fresh modeller. OpenMM PDBFile assigns elements from the
    # atom-name table; OXT is properly classified as O.
    new_pdb_text = "\n".join(out_lines) + "\n"
    sio = io.StringIO(new_pdb_text)
    reloaded = app.PDBFile(sio)
    # Replace the modeller's topology + positions in place. Modeller
    # exposes these as attributes; rebuilding from the reloaded PDB is
    # the cleanest way to fold the OXT insertions back into the engine
    # pipeline.
    modeller.topology = reloaded.topology
    modeller.positions = reloaded.positions

    return provenance


def _pick_platform(mm, say):
    for name in ("CUDA", "OpenCL", "CPU", "Reference"):
        try:
            p = mm.Platform.getPlatformByName(name)
            say(f"Using platform: {name}")
            return p
        except Exception:
            continue
    say("Falling back to default platform.")
    return None


def _dump_pdb(simulation, mm, unit, path: Path) -> None:
    """Write the current simulation state to a PDB."""
    from openmm.app import PDBFile  # safe — caller proved openmm is importable
    state = simulation.context.getState(getPositions=True)
    with path.open("w", encoding="utf-8") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)


def _dump_solute_pdb(
    simulation, path: Path, *, solute_topology, n_solute_atoms: int,
) -> None:
    """Write the solute portion of the current simulation state as a PDB.

    Used for explicit-solvent runs so the per-frame trajectory PDBs hold
    only receptor + cofactors + ligand — no waters, no neutralization ions.
    Per the [B2] design decision (frame format = solute-only): the analyzer,
    MM-GBSA estimator, and disk footprint all benefit; the full solvated
    state is still recoverable from the checkpoint, so explicit-mode
    publication artifacts lose nothing.

    `solute_topology` is the modeller's topology AS IT EXISTED BEFORE
    addSolvent (captured to a sidecar PDB and reloaded so it's a clean,
    independent Topology object). `n_solute_atoms` is the atom count at
    that same point. Because addSolvent appends water/ion atoms to the
    end of the topology, solute atoms always occupy indices
    [0, n_solute_atoms) in the simulation state — slicing positions
    by that range yields the solute coordinates without any per-atom
    membership check.
    """
    from openmm.app import PDBFile  # safe — caller proved openmm is importable
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    solute_positions = positions[:n_solute_atoms]
    with path.open("w", encoding="utf-8") as f:
        PDBFile.writeFile(solute_topology, solute_positions, f)


def _now_utc_iso() -> str:
    import datetime as _dt
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _write_checkpoint(
    simulation, out_dir: Path, *,
    engine_kind: str,
    final_frame_index: int,
    n_frames_this_leg: int,
    temperature_k: float,
    timestep_fs: float,
    friction_per_ps: float,
    snapshot_every_ps: float,
    production_ps: float,
    solvent: str = "implicit",
    water_model: str = "tip3p",
    water_padding_nm: float = 1.0,
    ionic_strength_molar: float = 0.15,
    pressure_bar: float = 1.0,
    barostat_frequency_steps: int = 25,
    npt_equilibration_ps: float = 0.0,
    position_restraint_k_kj_per_mol_per_nm2: float = 0.0,
    n_solute_atoms: int = 0,
    box_shape: str = "cube",
    modeller_topology=None,
    modeller_positions=None,
    system=None,
    mm=None,
    app=None,
    say,
) -> Tuple[Path, dict]:
    """Serialize the final State (positions + velocities + box + time) plus a
    metadata sidecar so a later run can continuously resume from here.

    Uses ``Simulation.saveState`` (portable XML) rather than ``saveCheckpoint``
    (binary, locked to the exact OpenMM build + platform).

    For explicit-solvent mode the function ALSO writes ``system.xml`` (the
    serialized System) and ``topology.pdb`` (the solvated topology with box
    vectors). On restart, addSolvent must NOT be re-run (it generates a
    different water shell every call); the [B3] restart path will load
    these two files plus final_state.xml to rebuild the identical context.
    """
    state_path = out_dir / CHECKPOINT_STATE_FILENAME
    simulation.saveState(str(state_path))

    if solvent == "explicit":
        ensemble = "NPT"
        implicit_solvent_label: Optional[str] = None
        force_fields = list(_explicit_forcefields(water_model))
        if system is not None and mm is not None:
            sys_path = out_dir / CHECKPOINT_SYSTEM_FILENAME
            sys_path.write_text(mm.XmlSerializer.serialize(system),
                                encoding="utf-8")
        if modeller_topology is not None and modeller_positions is not None and app is not None:
            top_path = out_dir / CHECKPOINT_TOPOLOGY_FILENAME
            with top_path.open("w", encoding="utf-8") as _f:
                app.PDBFile.writeFile(modeller_topology, modeller_positions, _f)
    else:
        ensemble = "NVT"
        implicit_solvent_label = "obc2"
        force_fields = list(_FORCE_FIELDS)

    meta: Dict[str, Any] = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "engine_kind": engine_kind,
        "state_xml": CHECKPOINT_STATE_FILENAME,
        "final_frame_index": int(final_frame_index),
        "n_frames_this_leg": int(n_frames_this_leg),
        "ensemble": ensemble,
        "integrator": "LangevinMiddleIntegrator",
        "implicit_solvent": implicit_solvent_label,
        "force_fields": force_fields,
        "small_molecule_forcefield": _SMALL_MOLECULE_FF,
        "temperature_k": float(temperature_k),
        "timestep_fs": float(timestep_fs),
        "friction_per_ps": float(friction_per_ps),
        "snapshot_every_ps": float(snapshot_every_ps),
        "production_ps_this_leg": float(production_ps),
        "solvent": solvent,
        "state_saved_at": _now_utc_iso(),
    }
    if solvent == "explicit":
        meta.update({
            "water_model": water_model,
            "water_padding_nm": float(water_padding_nm),
            "ionic_strength_molar": float(ionic_strength_molar),
            "pressure_bar": float(pressure_bar),
            "barostat_frequency_steps": int(barostat_frequency_steps),
            "npt_equilibration_ps_this_leg": float(npt_equilibration_ps),
            "position_restraint_k_kj_per_mol_per_nm2":
                float(position_restraint_k_kj_per_mol_per_nm2),
            "n_solute_atoms": int(n_solute_atoms),
            "system_xml": CHECKPOINT_SYSTEM_FILENAME,
            "topology_pdb": CHECKPOINT_TOPOLOGY_FILENAME,
            "box_shape": box_shape,
        })
    meta_path = out_dir / CHECKPOINT_META_FILENAME
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    extras = (f" + {CHECKPOINT_SYSTEM_FILENAME} + {CHECKPOINT_TOPOLOGY_FILENAME}"
              if solvent == "explicit" else "")
    say(f"Checkpoint written: {state_path.name} (+ {meta_path.name}{extras}), "
        f"final_frame_index={final_frame_index}, ensemble={ensemble}.")
    return state_path, meta


def _load_checkpoint_meta(checkpoint_xml: Path) -> dict:
    """Read the metadata sidecar that sits next to a checkpoint State XML."""
    meta_path = checkpoint_xml.parent / CHECKPOINT_META_FILENAME
    if not meta_path.exists():
        raise ValueError(
            f"checkpoint metadata sidecar not found next to {checkpoint_xml} "
            f"(expected {meta_path}); cannot validate restart compatibility."
        )
    return json.loads(meta_path.read_text(encoding="utf-8"))


def _validate_restart_compatibility(
    meta: dict, *,
    temperature_k: float,
    timestep_fs: float,
    friction_per_ps: float,
    snapshot_every_ps: float,
    solvent: str = "implicit",
    water_model: str = "tip3p",
    water_padding_nm: float = 1.0,
    ionic_strength_molar: float = 0.15,
    pressure_bar: float = 1.0,
    barostat_frequency_steps: int = 25,
    box_shape: str = "cube",
    say,
) -> None:
    """Refuse to resume unless the requested run conditions match the checkpoint.

    A continuous extension is only meaningful if the integrator, thermostat,
    timestep, and force-field stack are identical — otherwise the concatenated
    trajectory would not be a single thermodynamic ensemble. Raises ValueError
    listing every mismatch.

    Explicit-mode adds checks for water_model, water_padding_nm,
    ionic_strength_molar, pressure_bar, and barostat_frequency_steps. Cross-mode
    restart (implicit↔explicit) is rejected at the first check — the System
    layouts are fundamentally different (PME + waters vs. OBC2 + no waters)
    and there's no continuous-extension semantics that would unify them."""
    mismatches: List[str] = []

    def _chk_num(key: str, want: float, tol: float) -> None:
        got = meta.get(key)
        if got is None:
            mismatches.append(f"{key}: checkpoint has no value (requested {want})")
        elif abs(float(got) - float(want)) > tol:
            mismatches.append(f"{key}: checkpoint={got} != requested={want}")

    # Solvent mode comes first — a cross-mode mismatch makes every other check
    # meaningless, but emit it as a mismatch (not raise) so the caller sees
    # the full list in one shot.
    parent_solvent = meta.get("solvent")
    if parent_solvent is None:
        # Pre-1.1.0 checkpoints (implicit-only); treat as implicit for backcompat.
        parent_solvent = "implicit"
    if parent_solvent != solvent:
        mismatches.append(
            f"solvent: checkpoint={parent_solvent!r} != requested={solvent!r} "
            f"— cross-mode restart not supported"
        )

    _chk_num("temperature_k", temperature_k, 1e-6)
    _chk_num("timestep_fs", timestep_fs, 1e-9)
    _chk_num("friction_per_ps", friction_per_ps, 1e-9)
    _chk_num("snapshot_every_ps", snapshot_every_ps, 1e-9)

    expected_force_fields = (
        list(_explicit_forcefields(water_model)) if solvent == "explicit"
        else list(_FORCE_FIELDS)
    )
    if list(meta.get("force_fields") or []) != expected_force_fields:
        mismatches.append(
            f"force_fields: checkpoint={meta.get('force_fields')} "
            f"!= requested={expected_force_fields}")
    if meta.get("small_molecule_forcefield") != _SMALL_MOLECULE_FF:
        mismatches.append(
            f"small_molecule_forcefield: checkpoint={meta.get('small_molecule_forcefield')} "
            f"!= requested={_SMALL_MOLECULE_FF}")

    if solvent == "explicit":
        if meta.get("water_model") != water_model:
            mismatches.append(
                f"water_model: checkpoint={meta.get('water_model')!r} "
                f"!= requested={water_model!r}")
        _chk_num("water_padding_nm", water_padding_nm, 1e-9)
        _chk_num("ionic_strength_molar", ionic_strength_molar, 1e-9)
        _chk_num("pressure_bar", pressure_bar, 1e-9)
        # Barostat frequency is integer — exact-match.
        got_freq = meta.get("barostat_frequency_steps")
        if got_freq is None or int(got_freq) != int(barostat_frequency_steps):
            mismatches.append(
                f"barostat_frequency_steps: checkpoint={got_freq} "
                f"!= requested={barostat_frequency_steps}")
        got_shape = meta.get("box_shape", "cube")  # backcompat default
        if got_shape != box_shape:
            mismatches.append(
                f"box_shape: checkpoint={got_shape!r} "
                f"!= requested={box_shape!r}")

    if mismatches:
        raise ValueError(
            "restart is incompatible with the checkpoint (the extended leg would "
            "not be thermodynamically continuous):\n  - " + "\n  - ".join(mismatches)
        )
    say(f"Restart compatibility validated (solvent={solvent}, temperature, "
        f"timestep, friction, snapshot interval, force-field stack"
        + (", water model, padding, ionic strength, pressure, barostat freq"
           if solvent == "explicit" else "")
        + " all match the checkpoint).")
