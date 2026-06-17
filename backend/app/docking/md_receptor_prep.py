"""
Unified MD-grade receptor preparation — cofactors PARAMETERIZED, not stripped.

The receptor is made amber14-parameterizable AND cofactors are retained:
  - PROTEIN: PDBFixer models gaps, completes atoms, adds H at pH 7.4
  - METAL IONS (MG, ZN, CA, MN, FE, ...): stay in the receptor PDB,
    parameterized by amber14/tip3p.xml's ion residue templates at MD time
  - ORGANIC COFACTORS (GTP, GDP, ATP, NAD, ...): extracted to per-cofactor
    sidecar PDBs; SMILES looked up from a curated library; OpenFF Molecules
    built (RDKit `AssignBondOrdersFromTemplate` + `AddHs`) for SystemGenerator
    parameterization (GAFF-2.11) at MD time
  - VALIDATED end-to-end by building a SystemGenerator-backed system
    (amber14 + tip3p ions + obc2 + GAFF-2.11 for cofactors) and running
    create_system; on failure the output is NOT written.

If an organic cofactor isn't in the SMILES library MdReceptorPrepError is
raised — cofactors are never stripped silently.

Provenance written to metadata.json `receptor_prep`:
  - modeled_regions, chain_breaks_in_input
  - parameterized_cofactors: list of
      {chain, residue_name, residue_number, smiles, sidecar_pdb,
       atom_count, parameterization_method}
  - validated_with: forcefields + small_molecule_forcefield used
  - validated: true
"""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "MSE", "SEC",
}
_WATER = {"HOH", "WAT", "H2O", "TIP3", "TIP", "SOL"}

# Single-atom metal ions that amber14/tip3p.xml templates parameterize directly.
_METAL_IONS = {
    "MG", "ZN", "CA", "MN", "FE", "FE2", "FE3", "CO", "NI", "CU", "CD",
    "HG", "NA", "K", "CL", "LI", "RB", "CS", "SR", "BA",
}

# Curated SMILES for common biological cofactors. Heavy-atom topology must
# match what the PDB delivers (the standard PDB CCD definitions). Neutral
# forms — protonation/charges resolved by the force-field stack at MD time.
_COFACTOR_SMILES: Dict[str, str] = {
    "GTP": "Nc1nc2c(ncn2[C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1",
    "GDP": "Nc1nc2c(ncn2[C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1",
    "ATP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "ADP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "AMP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O",
}


class MdReceptorPrepError(RuntimeError):
    """Raised when MD-grade prep cannot produce a createSystem-passing receptor."""


# ---------------------------------------------------------------------
# Crystallographic-artifact metal-ion classifier
#
# Some structural files include metal ions that exist only because of the
# crystallization protocol — e.g. the Zn²⁺ in Nogales/Downing zinc-induced
# tubulin sheets (PDB 1JFF lineage), which holds protofilaments together in
# the 2D crystal but is absent in physiological microtubules. These ions
# are NOT part of the biological structure; carrying them into MD distorts
# the receptor's surface energetics and (for Zn²⁺ in particular) the
# fixed-charge 12-6 LJ parameterization is a poor model that's only "good
# enough" when the metal is doing nothing important.
#
# The classifier flags such ions by looking at their coordination shell:
#   - Functional metals (catalytic Zn, structural Ca, GTP/ATP-Mg, etc.) have
#     side-chain N/S/O donor contacts within ~3 Å (HIS Nε, CYS Sγ, ASP/GLU
#     OD/OE, etc.) OR organic-cofactor heavy-atom contacts within ~4 Å
#     (Mg²⁺ coordinating phosphate oxygens).
#   - Crystallographic packing ions have ONLY backbone (N/Cα/C/O) atoms in
#     their shell at >3 Å — characteristic of a loose electrostatic perch
#     against the protein surface rather than a true binding pocket.
#
# Classification is conservative: a metal must miss BOTH side-chain donors
# AND organic-cofactor neighborhood to be flagged. False negatives (real
# weakly-coordinated metals retained) are preferred to false positives
# (functional metals stripped).
#
# Operator override: receptor prep accepts `keep_metals` (a set of
# (chain, resname, resseq) tuples) that bypass the classifier — for the
# rare case where the heuristic mis-classifies a metal you know is real.
# ---------------------------------------------------------------------

# Side-chain N/S/O donor atom names per residue. These are the standard
# metal-coordination donors that distinguish a true binding pocket from a
# backbone-only packing contact.
_SIDECHAIN_METAL_DONORS: Dict[str, set] = {
    "HIS": {"ND1", "NE2"},
    "CYS": {"SG"},
    "MET": {"SD"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "ASN": {"OD1", "ND2"},
    "GLN": {"OE1", "NE2"},
    "ARG": {"NE", "NH1", "NH2"},
    "LYS": {"NZ"},
    "SER": {"OG"},
    "THR": {"OG1"},
    "TYR": {"OH"},
    # Phosphorylated / modified residues sometimes carry their own donors;
    # extend here as future inputs surface them.
}

# Classification cutoffs (Angstrom). 3.0 Å for side-chain donors matches the
# upper edge of physiological metal-coordination distances (catalytic Zn-N(His)
# is typically 2.0-2.3 Å; Mg-O(carboxylate) is 2.0-2.5 Å). 4.0 Å for organic
# cofactors catches Mg²⁺ at GTP/ATP phosphate sites where O1G-Mg is typically
# 2.0-2.5 Å but second-shell phosphate carbons sit at 3.5-4.0 Å.
_ARTIFACT_SIDECHAIN_DONOR_CUTOFF_A = 3.0
_ARTIFACT_COFACTOR_HEAVY_ATOM_CUTOFF_A = 4.0
# 4.5 Å neighborhood is reported for audit even when classification is "artifact"
# so the provenance shows WHICH backbone atoms it was perched against.
_ARTIFACT_NEIGHBORHOOD_REPORT_A = 4.5

# Bump when the classification rule changes so existing provenance entries
# remain interpretable in the catalog. Embedded in every strip record.
ARTIFACT_STRIP_RULE_VERSION = "1.0.0"


def _classify_metal_role(
    ion_xyz: Tuple[float, float, float],
    *,
    receptor_atoms: List[Dict[str, Any]],
    cofactor_heavy_atoms_xyz: List[Tuple[float, float, float]],
) -> Dict[str, Any]:
    """Classify a monoatomic metal ion as 'functional' or 'crystallographic artifact'.

    Pure function — no I/O, no openmm import. Returns a dict with the role
    label, the audit shell, and a human-readable reason. Stripping decisions
    live in the caller (see prepare_md_grade_receptor + engine_openmm).

    Args:
        ion_xyz: (x, y, z) of the metal ion in Angstroms.
        receptor_atoms: heavy-atom records for the receptor, each
            ``{"resname": str, "atom_name": str, "x": float, "y": float, "z": float}``.
            MUST exclude the ion being classified itself (and any H atoms).
        cofactor_heavy_atoms_xyz: heavy-atom positions for ALL organic cofactor
            atoms (GTP, GDP, ATP, NAD, …), gathered from the cofactor sidecars.

    Returns:
        ``{"role": "functional" | "artifact", "rule_version": str,
           "reason": str, "n_sidechain_donor_contacts": int,
           "n_cofactor_heavy_atom_contacts": int, "shell": list}``
    """
    sx, sy, sz = ion_xyz

    sidechain_donor_contacts: List[Tuple[float, str, str]] = []
    shell_atoms: List[Tuple[float, Dict[str, Any]]] = []
    sd_cut_sq = _ARTIFACT_SIDECHAIN_DONOR_CUTOFF_A ** 2
    shell_cut_sq = _ARTIFACT_NEIGHBORHOOD_REPORT_A ** 2

    for atom in receptor_atoms:
        dx, dy, dz = atom["x"] - sx, atom["y"] - sy, atom["z"] - sz
        d2 = dx * dx + dy * dy + dz * dz
        if d2 <= shell_cut_sq:
            shell_atoms.append((math.sqrt(d2), atom))
        donors = _SIDECHAIN_METAL_DONORS.get(atom["resname"], set())
        if atom["atom_name"] in donors and d2 <= sd_cut_sq:
            sidechain_donor_contacts.append(
                (math.sqrt(d2), atom["resname"], atom["atom_name"])
            )

    cof_cut_sq = _ARTIFACT_COFACTOR_HEAVY_ATOM_CUTOFF_A ** 2
    cof_contacts = 0
    for cx, cy, cz in cofactor_heavy_atoms_xyz:
        dx, dy, dz = cx - sx, cy - sy, cz - sz
        if dx * dx + dy * dy + dz * dz <= cof_cut_sq:
            cof_contacts += 1

    is_artifact = (not sidechain_donor_contacts) and (cof_contacts == 0)

    shell_summary = [
        {
            "distance_a": round(d, 2),
            "resname": a["resname"],
            "atom_name": a["atom_name"],
        }
        for d, a in sorted(shell_atoms, key=lambda da: da[0])[:8]
    ]

    if is_artifact:
        nearest = ", ".join(
            f"{a['resname']}.{a['atom_name']}@{d:.2f}Å"
            for d, a in sorted(shell_atoms, key=lambda da: da[0])[:5]
        )
        reason = (
            f"no side-chain N/S/O donor within "
            f"{_ARTIFACT_SIDECHAIN_DONOR_CUTOFF_A:.1f} Å and no organic-cofactor "
            f"heavy atom within {_ARTIFACT_COFACTOR_HEAVY_ATOM_CUTOFF_A:.1f} Å; "
            f"{len(shell_atoms)}-atom backbone/no-donor shell within "
            f"{_ARTIFACT_NEIGHBORHOOD_REPORT_A:.1f} Å"
            + (f" (nearest: {nearest})" if nearest else "")
        )
    else:
        parts: List[str] = []
        if sidechain_donor_contacts:
            ex = ", ".join(
                f"{r}.{a}@{d:.2f}Å" for d, r, a in sorted(sidechain_donor_contacts)[:3]
            )
            parts.append(
                f"{len(sidechain_donor_contacts)} side-chain donor(s) "
                f"≤{_ARTIFACT_SIDECHAIN_DONOR_CUTOFF_A:.1f} Å ({ex})"
            )
        if cof_contacts:
            parts.append(
                f"{cof_contacts} organic-cofactor heavy atom(s) "
                f"≤{_ARTIFACT_COFACTOR_HEAVY_ATOM_CUTOFF_A:.1f} Å"
            )
        reason = " + ".join(parts) + " — functional coordination"

    return {
        "role": "artifact" if is_artifact else "functional",
        "rule_version": ARTIFACT_STRIP_RULE_VERSION,
        "reason": reason,
        "n_sidechain_donor_contacts": len(sidechain_donor_contacts),
        "n_cofactor_heavy_atom_contacts": cof_contacts,
        "shell": shell_summary,
    }


def _parse_receptor_heavy_atoms_for_classifier(
    pdb_body: str, *, exclude_metal_keys: Optional[set] = None
) -> List[Dict[str, Any]]:
    """Pull heavy-atom records out of a PDB body for the metal classifier.

    Excludes hydrogens and (optionally) a set of metal-ion (chain, resname,
    resseq) keys so an ion isn't counted as its own shell neighbor when
    multiple metals are present.
    """
    out: List[Dict[str, Any]] = []
    exclude = exclude_metal_keys or set()
    for ln in pdb_body.splitlines():
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        if len(ln) < 54:
            continue
        try:
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
        except ValueError:
            continue
        elem = ln[76:78].strip() if len(ln) >= 78 else ""
        if elem.upper() == "H":
            continue
        resname = ln[17:20].strip().upper()
        atom_name = ln[12:16].strip()
        chain = ln[21] if len(ln) > 21 else " "
        resseq = ln[22:26].strip()
        if (chain, resname, resseq) in exclude:
            continue
        out.append({
            "resname": resname,
            "atom_name": atom_name,
            "x": x, "y": y, "z": z,
        })
    return out


def _read_seqres(pdb: Optional[Path]) -> List[str]:
    if not pdb or not pdb.exists():
        return []
    return [ln for ln in pdb.read_text(errors="ignore").splitlines()
            if ln.startswith("SEQRES")]


def _detect_internal_gaps(body: str) -> List[Dict[str, Any]]:
    per_chain: Dict[str, List[int]] = {}
    seen = set()
    for ln in body.splitlines():
        if not ln.startswith("ATOM  "):
            continue
        ch = ln[21]
        try:
            rs = int(ln[22:26])
        except ValueError:
            continue
        if (ch, rs) in seen:
            continue
        seen.add((ch, rs))
        per_chain.setdefault(ch, []).append(rs)
    gaps: List[Dict[str, Any]] = []
    for ch, nums in per_chain.items():
        for i in range(1, len(nums)):
            if nums[i] != nums[i - 1] + 1:
                gaps.append({
                    "chain": ch,
                    "after_residue": nums[i - 1],
                    "before_residue": nums[i],
                    "length": nums[i] - nums[i - 1] - 1,
                })
    return gaps


def _extract_residue_pdb_block(body: str, chain: str, resname: str, resnum: str) -> str:
    """Extract ATOM/HETATM lines for one residue (chain, resname, resnum)."""
    lines = []
    for ln in body.splitlines():
        if ln.startswith(("ATOM  ", "HETATM")):
            if (ln[17:20].strip() == resname
                    and (ln[21].strip() or " ") == chain
                    and ln[22:26].strip() == resnum):
                lines.append(ln)
    return "\n".join(lines) + "\nEND\n"


def _build_cofactor_molecule(pdb_block: str, smiles: str, label: str):
    """OpenFF Molecule from a heavy-atom PDB block + SMILES template.

    Uses RDKit AssignBondOrdersFromTemplate to copy SMILES-derived bond
    orders onto the PDB-parsed heavy-atom structure, then adds H with
    generated coordinates. The resulting Molecule has full topology
    (incl. H) and a conformer with the PDB heavy-atom coords + H coords.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from openff.toolkit import Molecule

    rdmol_pdb = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
    if rdmol_pdb is None:
        raise MdReceptorPrepError(f"{label}: RDKit could not parse PDB block")
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise MdReceptorPrepError(f"{label}: invalid SMILES {smiles!r}")
    try:
        rdmol = AllChem.AssignBondOrdersFromTemplate(template, rdmol_pdb)
    except Exception as e:
        raise MdReceptorPrepError(
            f"{label}: AssignBondOrdersFromTemplate failed "
            f"(heavy-atom topology mismatch between PDB and SMILES?): {e}"
        ) from e
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    try:
        Chem.SanitizeMol(rdmol)
    except Exception as e:
        raise MdReceptorPrepError(f"{label}: RDKit sanitize failed: {e}") from e
    try:
        return Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
    except Exception as e:
        raise MdReceptorPrepError(f"{label}: OpenFF Molecule.from_rdkit failed: {e}") from e


# Force-field stack used by SystemGenerator for both validation and MD.
_FORCEFIELDS = ["amber14-all.xml", "amber14/tip3p.xml", "implicit/obc2.xml"]
_SMALL_MOLECULE_FF = "gaff-2.11"


# ---------------------------------------------------------------------
# Item 10 (2026-06-01): atom-name normalization for downstream PDB→PDBQT tools.
#
# PDBFixer follows OpenMM/Amber atom-naming conventions that differ from
# PDB-standard / RDKit-residue-template conventions in subtle places. When a
# PDBFixer-produced PDB is fed to MGLTools' `prepare_receptor4.py` or Meeko's
# `mk_prepare_receptor.py`, the template-aware logic in each tool can fail
# (MGLTools: spurious "duplicate HG" collision; Meeko: RDKit AtomValenceException
# during _aux_altloc_mol_build). The fix is upstream of both: rewrite the file
# in place to use the conventions downstream tools expect.
#
# The transform table below encodes empirically-verified renames. The
# "minimal-now-extensible-later" philosophy: each entry encodes verified fact,
# not speculation. Add new entries as future receptor inputs surface failures.
#
# Each entry: (transform_name, predicate, rename_fn)
#   transform_name: short identifier (used in audit / logging)
#   predicate:      callable(atom_names_set) -> bool; True if transform applies
#   rename_fn:      callable(atom_names_set) -> dict[old_name, new_name]
# ---------------------------------------------------------------------
_ATOM_NAME_TRANSFORMS: List[tuple] = [
    (
        "n_terminal_nh3_plus",
        # PDBFixer's NH3+ N-terminus convention: the residue has both a lone
        # 'H' (the canonical backbone amide name) AND an 'H2' or 'H3', which
        # together form the trio of N-terminal hydrogens. PDB-standard +
        # RDKit residue templates expect H1/H2/H3. Rename H → H1.
        # Surfaced 2026-05-31 in item 10's tubulin-taxol failure; confirmed
        # via standalone Meeko CLI tests as the empirical root cause.
        lambda atom_names: ("H" in atom_names)
                           and (("H2" in atom_names) or ("H3" in atom_names)),
        lambda atom_names: {"H": "H1"},
    ),
    # --- Additional transforms added as failures surface them on future inputs. ---
]


def _format_atom_name_field(name: str) -> str:
    """Format an atom name into the PDB cols-13-16 4-char field.

    Standard PDB convention: 4-char names fill cols 13-16 with no leading
    space; shorter names (1-3 char) have a leading space + left-justified
    content in cols 14-16.
        'H'    -> ' H  '
        'H1'   -> ' H1 '
        'HG2'  -> ' HG2'
        'HG12' -> 'HG12'
    """
    if len(name) >= 4:
        return name[:4]
    return f" {name:<3}"


def _normalize_atom_names_for_downstream(pdb_path: Path) -> Dict[str, List[tuple]]:
    """Normalize PDBFixer's atom-naming idioms to conventions downstream
    PDB→PDBQT tools (MGLTools, Meeko) expect.

    Rewrites pdb_path in place. Returns a provenance dict
        {"<chain>:<resnum>": [(old_name, new_name), ...], ...}
    suitable for embedding under receptor_prep["normalized_atom_names"] in
    the prep summary for audit.

    Returns an empty dict when no transform applied (the receptor either had
    no N-terminus to rename, or used standard naming already).
    """
    text = pdb_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    # Pass 1: collect atom names per residue.
    residue_atoms: Dict[tuple, set] = {}
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21]
            resnum = line[22:26].strip()
            atom_name = line[12:16].strip()
            residue_atoms.setdefault((chain, resnum), set()).add(atom_name)

    # Pass 2: compute rename plan per residue by walking the transform table.
    rename_plan: Dict[tuple, Dict[str, str]] = {}
    audit: Dict[str, List[tuple]] = {}
    for key, names in residue_atoms.items():
        per_residue: Dict[str, str] = {}
        for _name, predicate, rename_fn in _ATOM_NAME_TRANSFORMS:
            if predicate(names):
                per_residue.update(rename_fn(names))
        if per_residue:
            rename_plan[key] = per_residue
            audit[f"{key[0]}:{key[1]}"] = sorted(
                [(old, new) for old, new in per_residue.items()]
            )

    if not rename_plan:
        return {}

    # Pass 3: rewrite the file applying renames.
    out_lines = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21]
            resnum = line[22:26].strip()
            atom_name = line[12:16].strip()
            key = (chain, resnum)
            if key in rename_plan and atom_name in rename_plan[key]:
                new_name = rename_plan[key][atom_name]
                new_field = _format_atom_name_field(new_name)
                line = line[:12] + new_field + line[16:]
        out_lines.append(line)

    new_text = "\n".join(out_lines)
    if text.endswith("\n"):
        new_text += "\n"
    pdb_path.write_text(new_text, encoding="utf-8")

    return audit


def prepare_md_grade_receptor(
    split_receptor_pdb: Path,
    out_receptor_pdb: Path,
    out_cofactors_dir: Path,
    *,
    reference_complex_pdb: Optional[Path] = None,
    ph: float = 7.4,
    keep_metals: Optional[set] = None,
) -> Dict[str, Any]:
    """Build an MD-grade receptor + parameterized cofactor sidecars.

    Returns provenance dict; raises MdReceptorPrepError (without writing
    output) on force-field validation failure or unknown cofactor.

    keep_metals: optional set of (chain, resname, resseq_str) tuples that
        bypass the crystallographic-artifact classifier. Use when the
        heuristic mis-classifies a metal you have reason to know is real.
    """
    from pdbfixer import PDBFixer
    from openmm import unit
    from openmm.app import CutoffNonPeriodic, HBonds, Modeller, PDBFile
    from openmmforcefields.generators import SystemGenerator

    if not split_receptor_pdb.exists() or split_receptor_pdb.stat().st_size == 0:
        raise MdReceptorPrepError(f"missing input receptor: {split_receptor_pdb}")

    out_receptor_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_cofactors_dir.mkdir(parents=True, exist_ok=True)

    # Re-attach SEQRES so PDBFixer can model internal gaps.
    seqres = _read_seqres(reference_complex_pdb)
    body = split_receptor_pdb.read_text(errors="ignore")
    work_pdb = out_receptor_pdb.parent / "_mdprep_input.pdb"
    work_pdb.write_text(
        ("\n".join(seqres) + "\n" if seqres else "") + body, encoding="utf-8"
    )

    chain_breaks = _detect_internal_gaps(body)

    # PDBFixer runs on PROTEIN ONLY (heterogens removed first) so its
    # findMissingResidues alignment is clean and addMissingAtoms properly
    # completes chain termini (OXT, terminal H). Cofactors are re-added
    # afterwards from the raw body — ions back into the receptor, organics
    # as sidecars with OpenFF Molecules.
    fixer = PDBFixer(filename=str(work_pdb))
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    modeled: List[Dict[str, Any]] = []
    for (ci, ins), resnames in sorted(fixer.missingResidues.items()):
        ch = chains[ci]
        rl = list(ch.residues())
        after = rl[ins - 1].id if ins > 0 else None
        before = rl[ins].id if ins < len(rl) else None
        modeled.append({
            "chain": ch.id,
            "after_residue": after, "before_residue": before,
            "length": len(resnames),
            "residue_names": list(resnames),
            "method": "pdbfixer.findMissingResidues+addMissingAtoms",
        })
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    # Protein-only modeller with gaps filled, termini capped, H added.
    modeller = Modeller(fixer.topology, fixer.positions)

    # --- Classify heterogens from the RAW body (single source of truth for
    #     cofactor identity): metals re-added to modeller; organics → sidecars ---
    parameterized: List[Dict[str, Any]] = []
    cofactor_mols: List[Any] = []
    ion_lines_by_key: Dict[tuple, List[str]] = {}
    organic_lines_by_key: Dict[tuple, List[str]] = {}
    seen_keys: set = set()
    for ln in body.splitlines():
        if not ln.startswith("HETATM"):
            continue
        rn = ln[17:20].strip().upper()
        ch = (ln[21].strip() or " ")
        rs = ln[22:26].strip()
        key = (ch, rn, rs)
        if rn in _WATER:
            continue
        if rn in _METAL_IONS:
            ion_lines_by_key.setdefault(key, []).append(ln)
        else:
            organic_lines_by_key.setdefault(key, []).append(ln)

    docking_root = out_receptor_pdb.parent  # e.g. jobs/<id>/docking

    # --- Crystallographic-artifact metal strip pass ---
    # Some structural inputs carry metals that exist only because of the
    # crystallization protocol (e.g. the Zn²⁺ in zinc-induced tubulin sheets).
    # Carrying them into MD distorts surface energetics, and the fixed-charge
    # 12-6 LJ parameterization for transition-metal ions is poorly suited to
    # any non-functional contact. Classify each candidate ion against its
    # coordination shell; strip artifacts before the receptor PDB is written.
    # `keep_metals` is the operator override for misclassifications.
    stripped_artifacts: List[Dict[str, Any]] = []
    if ion_lines_by_key:
        keep = set(keep_metals or ())
        cofactor_heavy_xyz: List[Tuple[float, float, float]] = []
        for lines in organic_lines_by_key.values():
            for ln in lines:
                if len(ln) < 54:
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
        # Build receptor heavy-atom records once, excluding the ions themselves
        # so an ion can't be its own coordination shell neighbor when multiple
        # metals share a chain.
        all_ion_keys = set(ion_lines_by_key.keys())
        receptor_heavy = _parse_receptor_heavy_atoms_for_classifier(
            body, exclude_metal_keys=all_ion_keys,
        )
        kept_ion_lines: Dict[tuple, List[str]] = {}
        for key, lines in ion_lines_by_key.items():
            ch, rn, rs = key
            if key in keep:
                kept_ion_lines[key] = lines
                continue
            line = lines[0]
            try:
                ion_xyz = (float(line[30:38]), float(line[38:46]),
                           float(line[46:54]))
            except ValueError:
                # Unparseable ion line — keep it (defensive: don't silently drop).
                kept_ion_lines[key] = lines
                continue
            classification = _classify_metal_role(
                ion_xyz,
                receptor_atoms=receptor_heavy,
                cofactor_heavy_atoms_xyz=cofactor_heavy_xyz,
            )
            if classification["role"] == "artifact":
                stripped_artifacts.append({
                    "chain": ch,
                    "residue_name": rn,
                    "residue_number": rs,
                    "xyz_input": [round(c, 3) for c in ion_xyz],
                    "rule_version": classification["rule_version"],
                    "reason": classification["reason"],
                    "shell": classification["shell"],
                    "n_sidechain_donor_contacts":
                        classification["n_sidechain_donor_contacts"],
                    "n_cofactor_heavy_atom_contacts":
                        classification["n_cofactor_heavy_atom_contacts"],
                })
                logger.info(
                    "Stripping crystallographic-artifact metal: %s %s%s — %s",
                    rn, ch, rs, classification["reason"],
                )
            else:
                kept_ion_lines[key] = lines
        ion_lines_by_key = kept_ion_lines

    # Re-add metal ions into the receptor modeller (amber14/tip3p templates
    # them). Build the topology programmatically — relying on PDB parsing of
    # the original HETATM line is unsafe because element-from-name inference
    # on a 2-char atom name like "ZN" gives 'Z'+'N' (nitrogen), not Zn.
    if ion_lines_by_key:
        from openmm import Vec3
        from openmm.app import Element, Topology
        from openmm.unit import angstrom
        ion_topo = Topology()
        ion_positions: List[Any] = []
        # one chain per chain-id to preserve chain identity
        ion_chains: Dict[str, Any] = {}
        for (ch, rn, rs), lines in ion_lines_by_key.items():
            line = lines[0]  # ions are single-atom residues
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            chain_obj = ion_chains.get(ch) or ion_topo.addChain(id=ch)
            ion_chains[ch] = chain_obj
            residue = ion_topo.addResidue(rn, chain_obj, id=rs)
            element = Element.getBySymbol(rn.capitalize())  # MG->Mg, ZN->Zn
            ion_topo.addAtom(rn, element, residue)
            ion_positions.append(Vec3(x, y, z) * angstrom)
            parameterized.append({
                "chain": ch, "residue_name": rn, "residue_number": rs,
                "smiles": None, "atom_count": 1, "sidecar_pdb": None,
                "parameterization_method": "amber14/tip3p ion template",
            })
        modeller.add(ion_topo, ion_positions)

    # Organic cofactors → sidecars + OpenFF Molecules.
    cofactor_sidecars: List[Path] = []
    for (ch, rn, rs), lines in organic_lines_by_key.items():
        if rn not in _COFACTOR_SMILES:
            raise MdReceptorPrepError(
                f"unknown organic cofactor {rn} (chain {ch} resseq {rs}): "
                f"add its SMILES to _COFACTOR_SMILES — cofactors must not be "
                f"silently stripped"
            )
        block = "\n".join(lines) + "\nEND\n"
        label = f"{rn}_{ch}{rs}"
        mol = _build_cofactor_molecule(block, _COFACTOR_SMILES[rn], label)
        sidecar = out_cofactors_dir / f"{label}.pdb"
        mol.to_file(str(sidecar), file_format="PDB")
        # OpenFF's PDB writer (via RDKit) emits heavy atoms as ATOM (correct
        # chain/resname/resseq) and H atoms as HETATM in a SEPARATE residue
        # named "UNL" with blank chain and resseq 1 — OpenMM's PDBFile reader
        # then sees two residues, only the heavy one matches the Molecule's
        # template, and createSystem fails on "missing 11 H atoms". Force
        # every atom into a single (HETATM, chain, resname, resseq) residue.
        _fixed = []
        for _ln in sidecar.read_text().splitlines():
            if _ln.startswith(("ATOM  ", "HETATM")):
                _ln = "HETATM" + _ln[6:]
                _ln = _ln[:17] + f"{rn:>3}" + _ln[20:]      # resname cols 18-20
                _ln = _ln[:21] + ch + _ln[22:]              # chain col 22
                _ln = _ln[:22] + f"{int(rs):>4}" + _ln[26:] # resseq cols 23-26
            _fixed.append(_ln)
        sidecar.write_text("\n".join(_fixed) + "\n")
        try:
            rel_sidecar = str(sidecar.relative_to(docking_root.parent))
        except ValueError:
            rel_sidecar = str(sidecar)
        parameterized.append({
            "chain": ch, "residue_name": rn, "residue_number": rs,
            "smiles": _COFACTOR_SMILES[rn],
            "atom_count": mol.n_atoms, "sidecar_pdb": rel_sidecar,
            "parameterization_method": f"{_SMALL_MOLECULE_FF} via SystemGenerator",
        })
        cofactor_mols.append(mol)
        cofactor_sidecars.append(sidecar)

    # --- Validate: build a system from (receptor + ions) + cofactor Molecules ---
    sg = SystemGenerator(
        forcefields=_FORCEFIELDS,
        small_molecule_forcefield=_SMALL_MOLECULE_FF,
        molecules=cofactor_mols,
        forcefield_kwargs={"constraints": HBonds},
        nonperiodic_forcefield_kwargs={
            "nonbondedMethod": CutoffNonPeriodic,
            "nonbondedCutoff": 1.0 * unit.nanometer,
        },
    )
    val_modeller = Modeller(modeller.topology, modeller.positions)
    # Add each cofactor by loading its sidecar PDB (which has H + CONECT
    # written by Molecule.to_file). This is the same path MD will use, so
    # validation here matches the runtime topology exactly.
    from openmm.app import PDBFile as _PDBFile
    for sidecar in cofactor_sidecars:
        pdbf = _PDBFile(str(sidecar))
        val_modeller.add(pdbf.topology, pdbf.positions)
    try:
        sg.create_system(val_modeller.topology)
    except Exception as e:
        raise MdReceptorPrepError(
            f"prepared receptor + cofactors failed SystemGenerator.create_system: {e}"
        ) from e

    # Validation passed — write the receptor PDB (protein + ions only).
    with out_receptor_pdb.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)

    # Item 10 fix (2026-06-01): normalize PDBFixer's atom-name idioms in
    # place. Downstream PDB→PDBQT tools (MGLTools, Meeko) consult RDKit
    # residue templates that don't tolerate PDBFixer's NH3+ N-terminus
    # convention ('H' for the lone backbone amide H alongside H2/H3); they
    # need 'H1'. See _normalize_atom_names_for_downstream + transform table.
    normalized_atom_names = _normalize_atom_names_for_downstream(out_receptor_pdb)

    try:
        work_pdb.unlink()
    except OSError:
        pass

    n_atoms = modeller.topology.getNumAtoms()
    logger.info(
        "MD-grade receptor: %d atoms | %d modeled regions | %d input breaks | "
        "%d parameterized cofactors (%d organic + %d ions) | "
        "%d residues with atom-name renames | %d stripped artifact metals",
        n_atoms, len(modeled), len(chain_breaks), len(parameterized),
        len(cofactor_mols), len(parameterized) - len(cofactor_mols),
        len(normalized_atom_names), len(stripped_artifacts),
    )
    return {
        "modeled_regions": modeled,
        "chain_breaks_in_input": chain_breaks,
        "parameterized_cofactors": parameterized,
        "stripped_artifacts": stripped_artifacts,
        "normalized_atom_names": normalized_atom_names,
        "validated_with": {
            "forcefields": _FORCEFIELDS,
            "small_molecule_forcefield": _SMALL_MOLECULE_FF,
        },
        "validated": True,
        "n_atoms": n_atoms,
        "ph": ph,
    }
