"""Sweep 3 explicit-TIP3P MD replicates × 3 taxol poses (Pose 0/1/2).

Per the 2026-06-16 replicate-experiment PRD. Goal: upgrade the single-run
intra-trajectory SEM (which is intrinsically too tight) into a proper
between-replicate σ on the MM-GBSA ΔG_bind for taxol on this pocket.

  9 fresh MD jobs, distinct from the 3 reference qualifying runs (never
  touched here):
      pose 0 → 5bc61f59 (reference; not modified)
      pose 1 → 34840aa1 (reference; not modified)
      pose 2 → 80e53d8a (reference; not modified)

  All replicates derive from docking job
  4a37bb0c-7655-411e-a238-1ff5b7bab910 (the canonical taxol-on-1JFF
  multi-ligand dock) and reuse its persisted GAFF/AM1-BCC template
  cache (the engine never re-derives charges in this sweep).

Independence (PRD requirement, option (a)):

  Each replicate runs a fresh minimisation + 50 ps NVT heat + 100 ps
  restrained NPT + restraint release + 1 ns production. Distinct
  `MdRunSettings.random_seed` per replicate decorrelates BOTH the
  Maxwell-Boltzmann initial velocity draw AND the Langevin noise
  through the whole trajectory — the two random sources `engine_openmm`
  plumbs the seed into. Same starting structure for all three
  replicates of a pose; different post-heat velocities + different
  noise sequences → different equilibration paths → statistically
  independent trajectories.

Seed scheme: SEED_BASE + pose_rank * 1000 + rep_idx. Deterministic +
recoverable from the manifest, distinct across all 9 runs by
construction.

Resumability: skips replicates whose `summary.json` already exists
(reuse the trajectory; the FE + aggregation steps re-key off the
manifest).

Outputs:

  jobs/<md_id>/md/...                              — per replicate
  docs/replicates_taxol_manifest.json              — (pose, rep)→md_id+seed
  docs/replicates_taxol_progress.md                — live progress table
"""

from __future__ import annotations

# CRITICAL (1/2): pre-load the conda libexpat BEFORE openmm / any backend
# import. `import openmm` pulls in a conflicting libexpat at runtime; once that
# copy wins the process-global expat symbols, the first `ElementTree.parse`
# (pyexpat) segfaults on the ABI mismatch — and OpenMM itself triggers that
# parse inside `PDBFile._loadNameReplacementTables`, so every MD dies at
# receptor load (after "Loading receptor ... and ligand ...", before "Using
# platform: CUDA", no Python traceback — a bare native SIGSEGV). ctypes-loading
# $CONDA_PREFIX/lib/libexpat.so.1 with RTLD_GLOBAL first binds the correct expat
# symbols globally, so pyexpat resolves against them no matter what openmm loads
# afterwards. This is the libexpat analogue of the libstdc++ load-order fix
# below (commit 15b3f16) and removes the need for a manual LD_PRELOAD.
import ctypes as _ctypes
import os as _os

_libexpat = _os.path.join(_os.environ.get("CONDA_PREFIX", ""),
                          "lib", "libexpat.so.1")
if _os.path.exists(_libexpat):
    _ctypes.CDLL(_libexpat, mode=_ctypes.RTLD_GLOBAL)

# CRITICAL (2/2): pre-load openmm BEFORE any backend import. The backend pulls
# in system libstdc++.so.6 (via psycopg/sqlalchemy), which lacks the
# CXXABI_1.3.15 symbol that openmm 8.x's libOpenMM.so requires. Once the
# system libstdc++ is loaded, openmm's import fails with ImportError, the
# engine raises OpenMMNotAvailable, and run_md_job silently falls back to
# the RDKit surrogate engine (41-frame MMFF wiggle) — a sneakily-wrong
# default for a replicate sweep meant to measure real explicit-TIP3P MD
# variability. Importing openmm first wins the libstdc++ ABI selection
# and lets CUDA / explicit-MD work end-to-end.
import openmm as _preload_openmm  # noqa: F401  pylint: disable=unused-import

import json
import os
import sys
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from backend.app.core.config import JOBS_DIR
from backend.app.core.job_artifacts import init_job
from backend.app.db import repository
from backend.app.md.job import MdRunSettings, md_dir_for, run_md_job

DOCKING_JOB_ID = "4a37bb0c-7655-411e-a238-1ff5b7bab910"
LIGAND_NAME = "lig_000000"
SEED_BASE = 20260616  # ISO date — deterministic + traceable
REPLICATES_PER_POSE = 3
POSES = [0, 1, 2]

REFERENCE_RUNS = {
    0: "5bc61f59-834f-4e71-a492-d32ddfdc7326",
    1: "34840aa1-cbf5-4c6b-a665-ea4f52110f5d",
    2: "80e53d8a-1926-4525-b2e1-55cb1e30eedd",
}

DOCS = Path("/home/xchem/projects/docking_app2/docs")
MANIFEST_PATH = DOCS / "replicates_taxol_manifest.json"
PROGRESS_MD = DOCS / "replicates_taxol_progress.md"

# Mirror the [B9] explicit-solvent protocol the reference runs used
# (5bc61f59 / 34840aa1 / 80e53d8a) so the replicate system is IDENTICAL to
# them. Two settings must match or the solvation box regresses:
#   - box_shape="octahedron": "cube" inflates the elongated tubulin system to
#     ~680k atoms (~3.5x); octahedron alone is ~522k.
#   - truncate_chain_ranges A439-451: cuts the disordered alpha-tubulin
#     C-terminal E-hook (VEGEGEEEGEEY), which otherwise extends far into solvent
#     and bloats the box. Removing it (OXT-capping A438) collapses the system to
#     the reference ~191k atoms (13,526 solute / 177,989 solvent / 191,515 total).
# With both matched, the ONLY deliberate per-replicate difference is random_seed.
SETTINGS_TEMPLATE = dict(
    production_ps=1000.0,
    snapshot_every_ps=5.0,
    temperature_k=300.0,
    solvent="explicit",
    water_model="tip3p",
    water_padding_nm=1.0,
    ionic_strength_molar=0.15,
    pressure_bar=1.0,
    barostat_frequency_steps=25,
    npt_equilibration_ps=100.0,
    position_restraint_k_kj_per_mol_per_nm2=1000.0,
    box_shape="octahedron",
    truncate_chain_ranges=[{"chain": "A", "from_resseq": 439, "to_resseq": 451}],
)


def make_settings(*, random_seed: int) -> MdRunSettings:
    return MdRunSettings(random_seed=random_seed, **SETTINGS_TEMPLATE)


def seed_for(pose_rank: int, rep_idx: int) -> int:
    return SEED_BASE + pose_rank * 1000 + rep_idx


def load_manifest() -> Dict[str, Dict[str, str]]:
    if MANIFEST_PATH.is_file():
        return json.loads(MANIFEST_PATH.read_text())
    return {}


def save_manifest(m: Dict[str, Dict[str, str]]) -> None:
    MANIFEST_PATH.write_text(json.dumps(m, indent=2, sort_keys=True))


def progress_table(manifest: Dict[str, Dict[str, str]]) -> str:
    lines = [
        "# Replicate sweep progress — taxol 3 reps × 3 poses",
        "",
        "| Pose | Rep | md_job_id | seed | status | n_frames | wall (s) |",
        "|---|---|---|---|---|---|---|",
    ]
    for pose in POSES:
        for rep in range(1, REPLICATES_PER_POSE + 1):
            key = f"pose{pose}_rep{rep}"
            entry = manifest.get(key, {})
            md_id = entry.get("md_job_id", "—")
            seed = entry.get("seed", "—")
            status = entry.get("status", "queued")
            n = entry.get("n_frames", "—")
            wall = entry.get("wall_seconds", "—")
            lines.append(
                f"| {pose} | {rep} | `{md_id[:8] if md_id != '—' else '—'}` "
                f"| {seed} | {status} | {n} | {wall} |"
            )
    return "\n".join(lines) + "\n"


def write_progress(manifest: Dict[str, Dict[str, str]]) -> None:
    PROGRESS_MD.write_text(progress_table(manifest))


def replicate_already_done(md_id: str) -> Optional[Dict]:
    """Returns the summary.json dict if a previous run completed."""
    p = JOBS_DIR / md_id / "md" / "summary.json"
    if not p.is_file():
        return None
    try:
        s = json.loads(p.read_text())
    except Exception:
        return None
    if s.get("status") == "completed":
        return s
    return None


def run_one_replicate(pose: int, rep: int, manifest: Dict[str, Dict[str, str]]) -> Dict:
    key = f"pose{pose}_rep{rep}"
    seed = seed_for(pose, rep)

    # Resume support — manifest holds the md_id we previously allocated.
    entry = manifest.get(key, {})
    md_id = entry.get("md_job_id")
    if md_id:
        existing = replicate_already_done(md_id)
        if existing is not None:
            print(f"  resuming {key}: existing {md_id[:8]} completed "
                  f"({existing.get('n_frames')} frames)")
            entry.update({
                "status": "completed",
                "n_frames": existing.get("n_frames"),
                "wall_seconds": round(existing.get("wall_seconds") or 0, 1),
                "verdict": existing.get("verdict"),
            })
            return entry

    # Fresh allocation.
    if not md_id:
        md_id = str(uuid.uuid4())

    print(f"\n=== {key}: md_job_id={md_id} seed={seed} ===")
    owner_id = repository.get_bootstrap_owner_id()

    # Mirror the API route's create-then-run flow.
    init_job(md_id)
    md_summary = {
        "stage": "md",
        "parent_docking_job_id": DOCKING_JOB_ID,
        "ligand": LIGAND_NAME,
        "pose_rank": pose,
        "replicate": {"pose": pose, "rep_idx": rep, "seed": seed,
                      "sweep": "replicates_taxol_2026-06-16"},
    }
    repository.set_job_status(md_id, "queued", summary_patch=md_summary)
    repository.create_job_row(
        job_id=md_id, owner_id=owner_id, kind="md", fs_root=md_id,
        status="queued", parent_job_id=DOCKING_JOB_ID,
        summary=md_summary,
    )
    repository.create_md_run(
        md_job_id=md_id,
        docking_job_id=DOCKING_JOB_ID,
        ligand_name=LIGAND_NAME,
        pose_rank=pose,
        production_ps=SETTINGS_TEMPLATE["production_ps"],
        snapshot_every_ps=SETTINGS_TEMPLATE["snapshot_every_ps"],
        temperature_k=SETTINGS_TEMPLATE["temperature_k"],
    )
    md_dir_for(md_id).mkdir(parents=True, exist_ok=True)

    settings = make_settings(random_seed=seed)

    # Update manifest now so a crash mid-run leaves the md_id recoverable.
    entry = {
        "md_job_id": md_id,
        "seed": seed,
        "pose": pose,
        "rep_idx": rep,
        "docking_job_id": DOCKING_JOB_ID,
        "ligand_name": LIGAND_NAME,
        "status": "running",
        "started_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    manifest[key] = entry
    save_manifest(manifest)
    write_progress(manifest)

    t0 = time.perf_counter()
    try:
        run_md_job(
            md_id,
            owner_id=str(owner_id) if owner_id else None,
            docking_job_id=DOCKING_JOB_ID,
            ligand=LIGAND_NAME,
            pose_rank=pose,
            settings=settings,
        )
    except Exception as e:
        wall = time.perf_counter() - t0
        print(f"  ✗ {key} FAILED after {wall:.0f}s — {type(e).__name__}: {e}")
        entry.update({"status": "failed", "error": str(e),
                      "wall_seconds": round(wall, 1)})
        save_manifest(manifest); write_progress(manifest)
        return entry
    wall = time.perf_counter() - t0

    s = json.loads((JOBS_DIR / md_id / "md" / "summary.json").read_text())
    entry.update({
        "status": s.get("status", "?"),
        "n_frames": s.get("n_frames"),
        "wall_seconds": round(wall, 1),
        "verdict": s.get("verdict"),
        "finished_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })
    save_manifest(manifest); write_progress(manifest)
    print(f"  ✓ {key} done in {wall:.0f}s ({s.get('n_frames')} frames; "
          f"verdict={s.get('verdict')})")
    return entry


def main():
    print(f"=== Taxol replicate sweep ===")
    print(f"  docking_job_id : {DOCKING_JOB_ID}")
    print(f"  ligand_name    : {LIGAND_NAME}")
    print(f"  poses          : {POSES}")
    print(f"  reps per pose  : {REPLICATES_PER_POSE}")
    print(f"  seed scheme    : {SEED_BASE} + pose*1000 + rep_idx")
    print(f"  reference runs : (NEVER touched) " + ", ".join(
        f"pose{p}={ref[:8]}" for p, ref in REFERENCE_RUNS.items()))
    print(f"  manifest       : {MANIFEST_PATH}")
    print()

    manifest = load_manifest()
    write_progress(manifest)

    for pose in POSES:
        for rep in range(1, REPLICATES_PER_POSE + 1):
            run_one_replicate(pose, rep, manifest)

    # Final summary.
    print("\n=== Sweep complete — manifest summary ===")
    for pose in POSES:
        for rep in range(1, REPLICATES_PER_POSE + 1):
            key = f"pose{pose}_rep{rep}"
            e = manifest.get(key, {})
            print(f"  {key:14s}  md_id={e.get('md_job_id','—')[:8]}  "
                  f"seed={e.get('seed','—')}  status={e.get('status','?')}  "
                  f"frames={e.get('n_frames','?')}  wall={e.get('wall_seconds','?')}s")
    print(f"\nManifest: {MANIFEST_PATH}")
    print(f"Progress : {PROGRESS_MD}")


if __name__ == "__main__":
    main()
