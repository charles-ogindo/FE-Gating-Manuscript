"""
Locate the receptor + ligand pose PDB files that an MD run will consume.

The MD stage reads from the docking job's viewer-exported PDBs (already in
PDB format, already pose-ranked) so no PDBQT->PDB conversion is needed here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from backend.app.core.config import JOBS_DIR


@dataclass
class PoseLocator:
    docking_job_id: str
    ligand: str
    pose_rank: int
    receptor_pdb: Path
    ligand_pdb: Path


class PoseNotFoundError(RuntimeError):
    pass


def resolve_pose(
    docking_job_id: str,
    ligand: str,
    pose_rank: int = 0,
) -> PoseLocator:
    """
    Return a fully-resolved PoseLocator pointing at PDB files on disk.

    Raises PoseNotFoundError with an actionable message if anything is missing.
    """
    job_dir = JOBS_DIR / docking_job_id
    if not job_dir.is_dir():
        raise PoseNotFoundError(f"docking job {docking_job_id!r} not found")

    viewer = job_dir / "docking" / "viewer"
    if not viewer.is_dir():
        raise PoseNotFoundError(
            f"viewer artifacts missing under {viewer.relative_to(JOBS_DIR)} — "
            "did docking complete?"
        )

    receptor = viewer / "receptor.pdb"
    if not receptor.is_file():
        raise PoseNotFoundError(
            "receptor.pdb missing under docking/viewer/ — re-export viewer artifacts"
        )

    lig_dir = viewer / ligand
    if not lig_dir.is_dir():
        raise PoseNotFoundError(
            f"ligand {ligand!r} has no viewer directory at "
            f"{lig_dir.relative_to(JOBS_DIR)}"
        )

    rank_pdb = lig_dir / f"{pose_rank:02d}.pdb"
    if not rank_pdb.is_file():
        available = sorted(p.stem for p in lig_dir.glob("*.pdb"))
        raise PoseNotFoundError(
            f"pose rank {pose_rank} missing for ligand {ligand!r}; "
            f"available ranks={available}"
        )

    return PoseLocator(
        docking_job_id=docking_job_id,
        ligand=ligand,
        pose_rank=pose_rank,
        receptor_pdb=receptor,
        ligand_pdb=rank_pdb,
    )


def list_available_ligands(docking_job_id: str) -> List[str]:
    """List ligand subfolders under docking/viewer/ for a completed docking job."""
    viewer = JOBS_DIR / docking_job_id / "docking" / "viewer"
    if not viewer.is_dir():
        return []
    return sorted(
        p.name for p in viewer.iterdir()
        if p.is_dir() and not p.name.startswith(".")
    )


def read_ligand_smiles(docking_job_id: str, ligand: str) -> Optional[str]:
    """Best-effort: pull SMILES from docking results so the MD summary can echo it."""
    job_dir = JOBS_DIR / docking_job_id
    candidates = [
        job_dir / "docking" / "scores.json",
        job_dir / "docking" / "dock_scores.json",
    ]
    for c in candidates:
        if not c.is_file():
            continue
        try:
            data = json.loads(c.read_text(encoding="utf-8"))
            rows = data.get("results", data) if isinstance(data, dict) else data
            for row in rows or []:
                if row.get("ligand") == ligand and row.get("smiles"):
                    return row["smiles"]
        except Exception:
            continue
    return None
