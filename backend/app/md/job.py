"""
MD job orchestrator.

`run_md_job(md_job_id, *, docking_job_id, ligand, pose_rank, ...)` is the
single entry point. It:

  1. Resolves the docked pose (receptor + ligand PDBs from the docking job).
  2. Picks an engine (OpenMM if importable; surrogate otherwise).
  3. Runs the engine to produce snapshot PDBs.
  4. Runs the analyzer to compute RMSD / H-bonds / contacts.
  5. Writes the canonical artifacts under jobs/<md_job_id>/md/.
  6. Updates the md_job's metadata.json with status + verdict.

Designed to be called from a FastAPI BackgroundTask. All exceptions are
caught and written as `summary.json.status = "failed"` so the frontend has
a single happy path to render.
"""

from __future__ import annotations

import json
import logging
import shutil
import time
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from backend.app.core.config import JOBS_DIR
from backend.app.core.job_artifacts import load_metadata
from backend.app.db import repository
from backend.app.md import (
    ARTIFACT_SCHEMA_VERSION,
    VERDICT_FAILED,
)
from backend.app.md.io import (
    PoseLocator,
    PoseNotFoundError,
    read_ligand_smiles,
    resolve_pose,
)
from backend.app.md.analyze import (
    parse_pdb,
    compute_rmsd_series,
    write_rmsd_csv,
    count_hbonds,
    write_hbonds_csv,
    per_residue_contact_frequency,
    write_contacts_csv,
    classify_stability,
    build_md_to_docking_resseq_map,
    relabel_contacts,
)

logger = logging.getLogger(__name__)


@dataclass
class MdRunSettings:
    production_ps: float = 500.0
    snapshot_every_ps: float = 5.0
    temperature_k: float = 300.0
    surrogate_frames: int = 41

    # ------------------------------------------------------------------
    # Solvent model (additive; default preserves the implicit-OBC2 path)
    # ------------------------------------------------------------------
    # "implicit" → unchanged behavior: GBSA-OBC2, no waters, no barostat,
    #   50 ps NVT heat, production NVT. Frame dumps include receptor +
    #   cofactors + ligand only (no waters present).
    # "explicit" → modeller.addSolvent (TIP3P) at `water_padding_nm`
    #   with `ionic_strength_molar` neutralization; PME + 1.0 nm cutoff
    #   + HBonds + rigidWater; MonteCarloBarostat. Staged equilibration:
    #     minimize → 50 ps restrained NVT heat
    #     → npt_equilibration_ps restrained NPT
    #     → release restraints → production NPT.
    #   Frame dumps are solute-only (waters + neutralization ions stripped
    #   at write time); biological metals/cofactors/ligand retained.
    solvent: Literal["implicit", "explicit"] = "implicit"

    # Explicit-only knobs (ignored when solvent == "implicit").
    water_model: str = "tip3p"
    water_padding_nm: float = 1.0
    ionic_strength_molar: float = 0.15
    pressure_bar: float = 1.0
    barostat_frequency_steps: int = 25
    npt_equilibration_ps: float = 100.0
    # Heavy-atom positional restraint (CustomExternalForce) on the solute
    # during heat + NPT equilibration. Released to k=0 at production start.
    position_restraint_k_kj_per_mol_per_nm2: float = 1000.0

    # Window the downstream gate should treat as equilibration. Default
    # None → gating sizes it automatically (50 ps heat + npt_equilibration_ps
    # + a ~20 ps restraint-release transient for explicit; 0 for implicit
    # since heat frames are never dumped). Override to pin a value when
    # comparing across runs.
    equilibration_discard_ps: Optional[float] = None

    # Periodic-box shape for explicit-solvent runs. "cube" reproduces the
    # default behavior; "octahedron" / "dodecahedron" cut wasted water for
    # elongated solutes (~25-30% atom-count reduction for tubulin-shaped
    # systems). Passed to modeller.addSolvent(boxShape=...). Ignored when
    # solvent == "implicit".
    box_shape: Literal["cube", "octahedron", "dodecahedron"] = "cube"

    # Truncate ranges of residues from the receptor BEFORE solvation. Each
    # entry: {"chain": str, "from_resseq": int, "to_resseq": int}. Useful
    # for cutting distal disordered tails (e.g., tubulin α C-terminal
    # E-hook) that bloat the explicit-solvent water shell without
    # affecting the binding site. Each truncation gets an NME amide cap
    # appended at the new C-terminus (residue from_resseq - 1), recorded
    # in `summary.receptor_truncation`. Empty list → no truncation
    # (default for backwards compatibility).
    truncate_chain_ranges: List[Dict[str, Any]] = field(default_factory=list)

    # Reproducible-replicate plumbing — see engine_openmm.run_openmm's
    # `random_seed`. When None (default) OpenMM's clock-based defaults
    # are used and the run is non-deterministic; set per replicate so
    # the trajectory is reproducible AND statistically independent of
    # sibling replicates (distinct integrator + initial-velocity draws).
    # Persisted into summary.settings.random_seed by run_md_job.
    random_seed: Optional[int] = None


def md_dir_for(md_job_id: str) -> Path:
    return JOBS_DIR / md_job_id / "md"


def summary_path(md_job_id: str) -> Path:
    return md_dir_for(md_job_id) / "summary.json"


def run_md_job(
    md_job_id: str,
    *,
    docking_job_id: str,
    ligand: str,
    pose_rank: int = 0,
    settings: Optional[MdRunSettings] = None,
    owner_id: Optional[str] = None,
    restart_from_md_job_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Run the MD validation. Caller owns the md_job_id (created upstream).

    Wrapped in job_failure_guard so an unhandled exception (including one raised
    before the inner try, e.g. in mkdir / status writes) always resolves the job
    to a terminal 'failed' status in both the filesystem and the DB rather than
    leaving a zombie at 'queued'/'running' (design note §5).

    ``restart_from_md_job_id`` continuously extends an existing OpenMM run: the
    parent must have a ``final_state.xml`` checkpoint (only runs produced after
    the checkpoint feature landed do). The parent's frames are copied into this
    job's frames dir and the engine appends new frames, so the combined
    trajectory clears the MM-GBSA frame-count gate. Restart is OpenMM-only —
    the surrogate engine has no velocity state to resume."""
    with repository.job_failure_guard(md_job_id, owner_id=owner_id):
        return _run_md_job_impl(
            md_job_id,
            docking_job_id=docking_job_id,
            ligand=ligand,
            pose_rank=pose_rank,
            settings=settings,
            restart_from_md_job_id=restart_from_md_job_id,
        )


def _run_md_job_impl(
    md_job_id: str,
    *,
    docking_job_id: str,
    ligand: str,
    pose_rank: int = 0,
    settings: Optional[MdRunSettings] = None,
    restart_from_md_job_id: Optional[str] = None,
) -> Dict[str, Any]:
    settings = settings or MdRunSettings()
    md_dir = md_dir_for(md_job_id)
    md_dir.mkdir(parents=True, exist_ok=True)

    repository.set_job_status_inflight(
        md_job_id,
        "running",
        set_started=True,
        summary_patch={"stage": "md", "started_at": _iso_now()},
    )

    started = time.perf_counter()
    pose: Optional[PoseLocator] = None
    engine_used: Optional[str] = None
    engine_attempts: list[str] = []

    try:
        pose = resolve_pose(docking_job_id, ligand, pose_rank)
        smiles = read_ligand_smiles(docking_job_id, ligand)

        # Read cofactor parameterization provenance from the parent docking
        # job's receptor_prep. Phase 6D: PG-first (jobs.summary['receptor_prep'])
        # with FS-fallback to metadata.json — defense-in-depth per Phase 6
        # GATE 1 decision 5. Handles 3 states:
        #   (1) PG row + key present (post-6C jobs, backfilled d96d719a).
        #   (2) PG row present but no receptor_prep key (12 jobs before
        #       backfill_summary.py ran, or zombies that lack the key on disk too).
        #   (3) PG row missing (legacy / pre-Phase-3 transient state).
        docking_job_dir = JOBS_DIR / docking_job_id
        cofactor_info: list = []
        try:
            parent = repository.get_job_row(docking_job_id)
            receptor_prep = (parent.summary or {}).get("receptor_prep") if parent else None
            if not receptor_prep:
                parent_meta = load_metadata(docking_job_id)
                receptor_prep = parent_meta.get("receptor_prep") if parent_meta else None
            cofactor_info = (receptor_prep or {}).get("parameterized_cofactors") or []
        except Exception as e:
            logger.warning("could not read receptor_prep from %s metadata: %s", docking_job_id, e)

        # ---- restart setup (optional continuous extension) ----
        # Resolve the parent checkpoint and stage the parent's frames into this
        # job's frames dir so the analyzer + frame-count gate see the FULL
        # logical trajectory (parent 000..offset + new offset+1..). The parent
        # job dir is never mutated.
        restart_from_checkpoint: Optional[Path] = None
        restart_frame_offset: int = 0
        if restart_from_md_job_id:
            parent_md_dir = md_dir_for(restart_from_md_job_id)
            ckpt = parent_md_dir / "final_state.xml"
            ckpt_meta = parent_md_dir / "final_state_meta.json"
            if not ckpt.exists() or not ckpt_meta.exists():
                raise FileNotFoundError(
                    f"parent MD job {restart_from_md_job_id} has no restart "
                    f"checkpoint ({ckpt} / {ckpt_meta}); only runs produced after "
                    "the checkpoint feature can be continuously extended."
                )
            restart_frame_offset = int(
                json.loads(ckpt_meta.read_text()).get("final_frame_index", 0)
            )
            new_frames_dir = md_dir / "frames"
            new_frames_dir.mkdir(parents=True, exist_ok=True)
            parent_frames = sorted((parent_md_dir / "frames").glob("frame_*.pdb"))
            for fp in parent_frames:
                shutil.copy2(fp, new_frames_dir / fp.name)
            restart_from_checkpoint = ckpt
            logger.info(
                "restart: staged %d parent frames from %s (offset=%d)",
                len(parent_frames), restart_from_md_job_id, restart_frame_offset,
            )

        # ---- pick engine ----
        result = None
        # 1) OpenMM
        try:
            from backend.app.md.engine_openmm import run_openmm, OpenMMNotAvailable
            try:
                engine_attempts.append("openmm")
                result = run_openmm(
                    receptor_pdb=pose.receptor_pdb,
                    ligand_pdb=pose.ligand_pdb,
                    out_dir=md_dir,
                    ligand_smiles=smiles,
                    cofactor_info=cofactor_info,
                    docking_job_dir=docking_job_dir,
                    production_ps=settings.production_ps,
                    snapshot_every_ps=settings.snapshot_every_ps,
                    temperature_k=settings.temperature_k,
                    restart_from_checkpoint=restart_from_checkpoint,
                    parent_md_job_id=restart_from_md_job_id,
                    frame_offset_hint=restart_frame_offset if restart_from_md_job_id else None,
                    # Explicit-solvent threading (implicit defaults unchanged).
                    solvent=settings.solvent,
                    water_model=settings.water_model,
                    water_padding_nm=settings.water_padding_nm,
                    ionic_strength_molar=settings.ionic_strength_molar,
                    pressure_bar=settings.pressure_bar,
                    barostat_frequency_steps=settings.barostat_frequency_steps,
                    npt_equilibration_ps=settings.npt_equilibration_ps,
                    position_restraint_k_kj_per_mol_per_nm2=
                        settings.position_restraint_k_kj_per_mol_per_nm2,
                    box_shape=settings.box_shape,
                    truncate_chain_ranges=settings.truncate_chain_ranges,
                    random_seed=settings.random_seed,
                )
                engine_used = result.engine_kind
            except OpenMMNotAvailable as e:
                engine_attempts.append(f"openmm:unavailable ({e})")
                result = None
        except Exception as e:
            logger.exception("OpenMM engine attempt failed (falling back to surrogate)")
            engine_attempts.append(f"openmm:error ({e})")
            result = None

        # A restart has no surrogate fallback — the surrogate engine holds no
        # velocity state to resume, so silently continuing there would produce a
        # discontinuous, mislabeled trajectory. Fail loudly instead.
        if result is None and restart_from_md_job_id:
            raise RuntimeError(
                "restart requires the OpenMM engine, which was unavailable; "
                f"attempts={engine_attempts}"
            )

        # 2) RDKit surrogate
        if result is None:
            try:
                from backend.app.md.engine_surrogate import run_surrogate
                engine_attempts.append("surrogate")
                result = run_surrogate(
                    receptor_pdb=pose.receptor_pdb,
                    ligand_pdb=pose.ligand_pdb,
                    out_dir=md_dir,
                    n_frames=settings.surrogate_frames,
                    snapshot_every_ps=settings.snapshot_every_ps,
                )
                engine_used = result.engine_kind
            except Exception as e:
                engine_attempts.append(f"surrogate:error ({e})")
                raise

        # ---- analyze snapshots ----
        # For a restart, the engine returns ONLY the new leg's frames, but the
        # parent's frames were staged into md_dir/frames during restart setup —
        # so analyze the full contiguous trajectory off disk (frame index i sits
        # at i * snapshot_every_ps, the same convention the engine uses). Fresh
        # runs analyze exactly what the engine returned.
        if restart_from_md_job_id and engine_used and engine_used.startswith("openmm"):
            frame_paths = sorted((md_dir / "frames").glob("frame_*.pdb"))
            frame_times = [i * settings.snapshot_every_ps for i in range(len(frame_paths))]
        else:
            frame_paths = list(result.snapshots)
            frame_times = list(result.times_ps)

        snapshots = [parse_pdb(p) for p in frame_paths]
        rmsd_series = compute_rmsd_series(snapshots, frame_times)
        write_rmsd_csv(md_dir / "rmsd.csv", rmsd_series)

        hbond_counts = [(t, count_hbonds(s))
                        for t, s in zip(frame_times, snapshots)]
        if hbond_counts:
            counts_only = [c for _, c in hbond_counts]
            persistence = {
                "mean": sum(counts_only) / len(counts_only),
                "min": min(counts_only),
                "max": max(counts_only),
                "frac_with_any": sum(1 for c in counts_only if c > 0) / len(counts_only),
            }
        else:
            persistence = {"mean": 0.0, "min": 0, "max": 0, "frac_with_any": 0.0}
        write_hbonds_csv(md_dir / "hbonds.csv", hbond_counts, persistence)

        contacts = per_residue_contact_frequency(snapshots)
        # Q6b PART 2: relabel MD's contiguous chain numbering back to the
        # docking receptor's author numbering. The MD prep (PDBFixer /
        # OpenMM normalization) drops resseq gap markers and renumbers
        # chains 1..N; the contacts/consensus pipelines + the tubulin
        # literature label by author numbering (PHE272, LEU275, GLY370 on
        # β-tubulin's taxane site). Without this relabel, top_contacts
        # land on LEU272/THR273/ARG275 — the same residues with shifted
        # numbers, confusing the operator. relabel_contacts is a no-op
        # if the map can't be built deterministically (per-chain length
        # mismatch or residue-name disagreement) — see analyze.py.
        renumber_map = build_md_to_docking_resseq_map(
            md_pdb=frame_paths[0],
            docking_pdb=pose.receptor_pdb,
        )
        contacts = relabel_contacts(contacts, renumber_map)
        write_contacts_csv(md_dir / "contacts.csv", contacts)

        verdict = classify_stability(rmsd_series, hbond_counts, contacts)

        wall_s = time.perf_counter() - started

        # Checkpoint / restart provenance. Every OpenMM run now serializes a
        # final State; expose it (+ a resumable flag) so downstream tooling
        # knows an under-sampled run can be EXTENDED rather than re-run from
        # scratch. restart_info is None for fresh runs, a provenance dict for
        # continuations.
        checkpoint_meta = getattr(result, "checkpoint_meta", None)
        restart_info = getattr(result, "restart_info", None)
        receptor_strip = getattr(result, "receptor_strip", None)
        checkpoint_block = None
        if checkpoint_meta:
            checkpoint_block = {**checkpoint_meta, "path": "md/final_state.xml"}

        summary = {
            "schema_version": ARTIFACT_SCHEMA_VERSION,
            "md_job_id": md_job_id,
            "docking_job_id": docking_job_id,
            "ligand": ligand,
            "pose_rank": pose_rank,
            "smiles": smiles,
            "engine": {
                "kind": engine_used,
                "attempts": engine_attempts,
            },
            "settings": asdict(settings),
            "n_frames": len(frame_paths),
            # Restart support: resumable=True means a final_state.xml checkpoint
            # exists and this run can be continuously extended (engine_openmm
            # checkpoint/restart). restart documents this run's own lineage.
            "resumable": checkpoint_block is not None,
            "checkpoint": checkpoint_block,
            "restart": restart_info,
            # Engine-side artifact-metal strip provenance. None when the
            # receptor was clean (no monoatomic metals, or already filtered
            # at receptor-prep time). Non-empty list documents which ions
            # the engine removed, with rule_version + reason per ion.
            "receptor_strip": receptor_strip,
            # Tail-truncation provenance. None when no truncation was
            # configured. Each entry documents the chain, range, the NME
            # cap insertion site, and the residue count.
            "receptor_truncation": getattr(result, "receptor_truncation", None),
            "wall_seconds": round(wall_s, 2),
            "status": "completed",
            "verdict": verdict.verdict,
            "rationale": verdict.rationale,
            "metrics": {
                "rmsd_backbone_final_a":         verdict.rmsd_backbone_final_a,
                # Q6b primary metrics (receptor-frame pose displacement —
                # what soft gate E reads, what classify_stability keys on).
                "rmsd_ligand_pose_final_a":      verdict.rmsd_ligand_pose_final_a,
                "rmsd_ligand_pose_max_a":        verdict.rmsd_ligand_pose_max_a,
                # Q6b diagnostic (legacy ligand-on-ligand Kabsch, renamed):
                "rmsd_ligand_internal_final_a":  verdict.rmsd_ligand_internal_final_a,
                "rmsd_ligand_internal_max_a":    verdict.rmsd_ligand_internal_max_a,
                "hbond_persistence_frac":        verdict.hbond_persistence_frac,
            },
            "top_contacts": [
                {"chain": c, "resseq": r, "resname": n, "frac": frac}
                for (c, r, n, frac) in verdict.top_contacts
            ],
            # Q6b PART 2 provenance: surface whether the MD→docking
            # author-numbering relabel ran. False means the map could
            # not be built (per-chain length mismatch or residue-name
            # disagreement) and labels are still MD-numbered — operator
            # should compare manually to docking/dockready_receptor.pdb.
            "receptor_renumbered": renumber_map is not None,
            "artifacts": {
                "summary":   "md/summary.json",
                "rmsd_csv":  "md/rmsd.csv",
                "hbonds_csv":"md/hbonds.csv",
                "contacts_csv":"md/contacts.csv",
                "log":       "md/log.txt",
                "frames_dir":"md/frames/",
                # Present only when an OpenMM checkpoint was written.
                **({"checkpoint_state": "md/final_state.xml",
                    "checkpoint_meta": "md/final_state_meta.json"}
                   if checkpoint_block else {}),
            },
            "free_energy": {
                "status": "planned",
                "reason": "MM/GBSA + MM/PBSA require an explicit MD/parameterization stack "
                          "(AmberTools or OpenMM-with-explicit-solvent) that isn't bundled "
                          "with this repo. Re-enable once the engine is in place.",
            },
        }
        summary_path(md_job_id).write_text(json.dumps(summary, indent=2),
                                           encoding="utf-8")
        repository.set_job_status_inflight(
            md_job_id,
            "completed",
            set_completed=True,
            summary_patch={
                "stage": "md",
                "md_verdict": verdict.verdict,
                "md_engine": engine_used,
                "wall_seconds": summary["wall_seconds"],
            },
        )
        repository.update_md_run(
            md_job_id,
            engine=engine_used,
            verdict=verdict.verdict,
            wall_seconds=summary["wall_seconds"],
        )
        return summary

    except PoseNotFoundError as e:
        return _write_failure(md_job_id, docking_job_id, ligand, pose_rank,
                              engine_attempts, str(e),
                              category="missing_input")
    except Exception as e:
        logger.exception("MD job failed")
        return _write_failure(md_job_id, docking_job_id, ligand, pose_rank,
                              engine_attempts,
                              f"{e.__class__.__name__}: {e}\n{traceback.format_exc()}",
                              category="engine_error")


def _write_failure(md_job_id, docking_job_id, ligand, pose_rank,
                   attempts, message, *, category: str):
    payload = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "md_job_id": md_job_id,
        "docking_job_id": docking_job_id,
        "ligand": ligand,
        "pose_rank": pose_rank,
        "engine": {"kind": None, "attempts": attempts},
        "status": "failed",
        "verdict": VERDICT_FAILED,
        "error_category": category,
        "error": message,
    }
    p = summary_path(md_job_id)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    repository.set_job_status_inflight(
        md_job_id,
        "failed",
        message=message,
        summary_patch={"stage": "md", "error": message[:500]},
    )
    repository.update_md_run(md_job_id, verdict=VERDICT_FAILED)
    return payload


def _iso_now():
    import datetime as _dt
    return _dt.datetime.now().isoformat(timespec="seconds")
