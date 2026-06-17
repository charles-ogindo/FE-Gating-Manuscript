"""
Confirms backend.app.md.analyze.parse_pdb excludes waters + neutralization
ions from both the receptor and ligand atom lists, which is the load-bearing
guarantee that pose-RMSD (and the downstream classify_stability verdict)
isn't contaminated by explicit-solvent water/ion atoms.

The explicit-solvent engine [B2] writes solute-only frames so waters
shouldn't reach the analyzer in the happy path, but a) any future variant
that dumps full-solvent frames must still produce correct RMSD numbers,
and b) the pre-existing _NON_LIGAND_HETATMS exclusion is what the [B5]
gating + provenance assertion ("Confirm pose-RMSD excludes water/ions")
explicitly relies on. This test pins that behavior so a future innocent-
looking rename in analyze.py can't silently flip it.
"""

from __future__ import annotations

import textwrap

from backend.app.md.analyze import parse_pdb


def _write(p, content):
    p.write_text(textwrap.dedent(content))


def test_parse_pdb_excludes_waters_from_both_receptor_and_ligand(tmp_path):
    """HOH/WAT/TIP3 atoms in a frame must be classified as neither
    receptor nor ligand — they should be dropped from BOTH lists."""
    pdb = tmp_path / "frame_explicit.pdb"
    _write(pdb, """\
        ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
        HETATM    3  C1  LIG X   1       5.000   5.000   5.000  1.00  0.00           C
        HETATM    4  O   HOH B   1       7.000   7.000   7.000  1.00  0.00           O
        HETATM    5  H1  HOH B   1       7.500   7.500   7.500  1.00  0.00           H
        HETATM    6  H2  HOH B   1       7.500   6.500   7.500  1.00  0.00           H
        HETATM    7  O   WAT B   2      -7.000  -7.000  -7.000  1.00  0.00           O
        HETATM    8  O   TIP3 B  3       3.000   3.000   3.000  1.00  0.00           O
        END
        """)
    snap = parse_pdb(pdb)
    # Receptor = ATOM lines only (NOT waters, NOT ligand, NOT cofactors).
    receptor_resnames = {a["resname"] for a in snap.receptor_atoms}
    assert "ALA" in receptor_resnames
    assert "HOH" not in receptor_resnames
    assert "WAT" not in receptor_resnames
    assert "TIP3" not in receptor_resnames
    # Ligand = HETATM minus _NON_LIGAND_HETATMS (waters, cofactors, metals).
    ligand_resnames = {a["resname"] for a in snap.ligand_atoms}
    assert ligand_resnames == {"LIG"}
    assert "HOH" not in ligand_resnames
    assert "WAT" not in ligand_resnames


def test_parse_pdb_excludes_neutralization_ions_from_both_lists(tmp_path):
    """NA/CL/K — the typical addSolvent neutralization ions — must also be
    excluded. They're in _METAL_IONS, which is a subset of _NON_LIGAND_HETATMS."""
    pdb = tmp_path / "frame_ions.pdb"
    _write(pdb, """\
        ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
        HETATM    2  C1  LIG X   1       2.000   2.000   2.000  1.00  0.00           C
        HETATM    3 NA   NA  B   1       5.000   0.000   0.000  1.00  0.00          Na
        HETATM    4 CL   CL  B   2      -5.000   0.000   0.000  1.00  0.00          Cl
        HETATM    5  K   K   B   3       0.000   5.000   0.000  1.00  0.00           K
        END
        """)
    snap = parse_pdb(pdb)
    receptor_resnames = {a["resname"] for a in snap.receptor_atoms}
    ligand_resnames = {a["resname"] for a in snap.ligand_atoms}
    # Neutralization ions absent from both sides.
    for ion in ("NA", "CL", "K"):
        assert ion not in receptor_resnames
        assert ion not in ligand_resnames
    # Sanity: the legitimate solute atoms ARE present.
    assert "ALA" in receptor_resnames
    assert "LIG" in ligand_resnames


def test_parse_pdb_keeps_cofactors_off_ligand_list(tmp_path):
    """Cofactors (GTP/GDP/...) must not be misclassified as the docked ligand.
    _COFACTOR_SMILES.keys() ⊂ _NON_LIGAND_HETATMS."""
    pdb = tmp_path / "frame_cofactors.pdb"
    _write(pdb, """\
        ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
        HETATM    2  C1  LIG X   1       2.000   2.000   2.000  1.00  0.00           C
        HETATM    3  PA  GTP A 500      -5.000   0.000   0.000  1.00  0.00           P
        HETATM    4  PA  GDP B 600       5.000   0.000   0.000  1.00  0.00           P
        END
        """)
    snap = parse_pdb(pdb)
    ligand_resnames = {a["resname"] for a in snap.ligand_atoms}
    assert ligand_resnames == {"LIG"}
    # Cofactors are HETATMs but not the docked ligand.
    assert "GTP" not in ligand_resnames
    assert "GDP" not in ligand_resnames
