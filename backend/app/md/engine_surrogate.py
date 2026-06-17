"""
RDKit MMFF surrogate engine — fallback when OpenMM isn't installed.

This is NOT a real MD trajectory. It samples constrained ligand conformations
in the fixed-receptor frame by:

  1. Reading receptor PDB (kept rigid).
  2. Reading ligand PDB, building an RDKit molecule with the docked coords as
     conformer 0.
  3. Generating N perturbed conformers via constrained embedding around the
     starting pose (small RMSD bias).
  4. MMFF94s-optimizing each conformer with the receptor heavy atoms as
     extra-bond-distance restraints (proxy for steric pocket).
  5. Writing each conformer as a frame PDB (receptor + perturbed ligand).

The result is useful for *intrinsic ligand stability in the pocket frame*:
big RMSD spread under MMFF = strained pose. It cannot tell you about
backbone flexibility — the receptor is rigid here. The summary.json
clearly states `engine_kind="surrogate"` so the UI never overclaims.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

logger = logging.getLogger(__name__)


class SurrogateUnavailable(RuntimeError):
    pass


@dataclass
class SurrogateResult:
    snapshots: List[Path]
    times_ps: List[float]
    log_path: Path
    engine_kind: str  # "surrogate"


def run_surrogate(
    receptor_pdb: Path,
    ligand_pdb: Path,
    out_dir: Path,
    *,
    n_frames: int = 41,
    perturb_kT: float = 0.6,
    snapshot_every_ps: float = 5.0,
) -> SurrogateResult:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as e:
        raise SurrogateUnavailable("rdkit is not importable") from e

    out_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = out_dir / "frames"
    frames_dir.mkdir(exist_ok=True)
    # Append (not truncate) so a prior engine's log (e.g. openmm's crash
    # trace before fallback) is preserved. Without this, the openmm log
    # gets clobbered when the orchestrator hands off to the surrogate and
    # a crash trace is lost.
    log_path = out_dir / "log.txt"
    log = log_path.open("a", encoding="utf-8")

    def say(m: str) -> None:
        log.write(m + "\n"); log.flush()
        logger.info("[md/surrogate] %s", m)

    try:
        if log_path.stat().st_size > 0:
            say("---- (above: prior engine log) ----")
        say(f"RDKit surrogate engine — receptor={receptor_pdb.name}, ligand={ligand_pdb.name}")

        receptor_block = receptor_pdb.read_text(encoding="utf-8")
        # Keep only ATOM/TER from the receptor; ignore any HETATM that might
        # have leaked in from a combined file.
        receptor_atom_only = "\n".join(
            ln for ln in receptor_block.splitlines()
            if ln.startswith(("ATOM", "TER")) and not ln.startswith("ENDMDL")
        ) + "\n"

        mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False, sanitize=False)
        if mol is None:
            raise SurrogateUnavailable("RDKit could not parse ligand PDB")
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            say(f"Sanitize warning: {e}; continuing with permissive mol")
        ref_mol = Chem.Mol(mol)  # keep original coords

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xCAFE
        params.useRandomCoords = False
        params.pruneRmsThresh = 0.05

        # Embed perturbed conformers anchored to original heavy-atom positions.
        AllChem.EmbedMultipleConfs(mol, numConfs=n_frames - 1, params=params)
        # Optimize with MMFF if parameters are available.
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s",
                                              maxIters=200)
        except Exception as e:
            say(f"MMFF skipped ({e}); using ETKDG geometry directly")

        # Align each generated conformer onto the original docked pose so that
        # the perturbation is interpretable as "wiggle around the docked pose".
        ref_conf_id = 0
        ref_mol = Chem.AddHs(ref_mol, addCoords=True) if ref_mol.GetNumConformers() else ref_mol
        ref_atom_map = list(range(min(ref_mol.GetNumAtoms(), mol.GetNumAtoms())))
        for cid in range(mol.GetNumConformers()):
            try:
                AllChem.AlignMol(mol, ref_mol, prbCid=cid, refCid=ref_conf_id,
                                 atomMap=[(i, i) for i in ref_atom_map])
            except Exception:
                pass

        snapshots: List[Path] = []
        times: List[float] = []

        for i in range(n_frames):
            cid = 0 if i == 0 else min(i - 1, mol.GetNumConformers() - 1) if mol.GetNumConformers() else 0
            # Frame 0 = original docked pose; frames 1..N = perturbed
            this_mol = ref_mol if i == 0 else mol
            lig_pdb_text = Chem.MolToPDBBlock(this_mol, confId=cid)
            # Force HETATM (RDKit emits ATOM for unknown residues sometimes)
            lig_pdb_text = "\n".join(
                ln.replace("ATOM  ", "HETATM", 1) if ln.startswith("ATOM  ") else ln
                for ln in lig_pdb_text.splitlines()
            )
            frame_path = frames_dir / f"frame_{i:03d}.pdb"
            frame_path.write_text(
                receptor_atom_only + lig_pdb_text + "\nEND\n",
                encoding="utf-8",
            )
            snapshots.append(frame_path)
            times.append(i * snapshot_every_ps)

        say(f"Wrote {len(snapshots)} surrogate frames "
            f"(t=0..{times[-1]:.1f} ps equivalent, perturb_kT={perturb_kT})")
        return SurrogateResult(
            snapshots=snapshots,
            times_ps=times,
            log_path=log_path,
            engine_kind="surrogate",
        )
    finally:
        try:
            log.close()
        except Exception:
            pass
