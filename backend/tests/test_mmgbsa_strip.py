"""
Tests for the MM-GBSA per-frame water/ion strip in
backend.app.free_energy.mmgbsa_runner._collect_frame_positions.

Hermetic: tiny synthetic PDB frames + a small synthetic Topology built via
the openmm.app primitives (no MD, no SystemGenerator, no sqm). Validates
that the strip:

  - filters HOH / WAT / TIP3 + neutralization ions (NA/CL/K) out of an
    explicit-solvent frame so the returned position array matches the
    GBSA complex topology atom-for-atom.
  - is a no-op when frames already contain only solute (the [B2]
    solute-only-frame writer's happy path) — counts match, positions
    pass through.
  - raises loudly when a frame is missing solute atoms the GBSA
    topology expects (catches receptor/atom-name drift between MD
    prep and GBSA build).
  - falls back to legacy "no mask" behavior when expected_topology=None.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pytest

from backend.app.free_energy.mmgbsa_runner import _collect_frame_positions


# ---------------------------------------------------------------------------
# Synthetic helpers
# ---------------------------------------------------------------------------

def _make_solute_topology():
    """Build a tiny ALA-residue topology (5 heavy atoms: N, CA, C, O, CB).

    Just enough atoms to exercise the mask without pulling in OpenMM's
    residue-template library. The classifier doesn't care about element
    correctness — it walks (chain, resname, resseq, atom_name) tuples.
    """
    from openmm.app import Topology, Element
    top = Topology()
    chain = top.addChain("A")
    res = top.addResidue("ALA", chain, id="1")
    top.addAtom("N",  Element.getBySymbol("N"), res)
    top.addAtom("CA", Element.getBySymbol("C"), res)
    top.addAtom("C",  Element.getBySymbol("C"), res)
    top.addAtom("O",  Element.getBySymbol("O"), res)
    top.addAtom("CB", Element.getBySymbol("C"), res)
    return top


def _write_solute_only_frame(frame_path: Path) -> None:
    """[B2]-format frame: 5 ALA heavy atoms, no waters, no ions."""
    pdb = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1       2.300   1.000   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1       2.000   2.100   0.000  1.00  0.00           O
        ATOM      5  CB  ALA A   1       1.900  -1.000   0.500  1.00  0.00           C
        END
        """)
    frame_path.write_text(pdb)


def _write_solvated_frame(frame_path: Path) -> None:
    """A frame as a future variant might write it: the 5 ALA atoms + 2 HOH
    waters + 1 NA neutralization ion. The mask must keep the 5 solute
    atoms and drop the 3 solvent atoms."""
    pdb = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1       2.300   1.000   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1       2.000   2.100   0.000  1.00  0.00           O
        ATOM      5  CB  ALA A   1       1.900  -1.000   0.500  1.00  0.00           C
        HETATM    6  O   HOH B   1       5.000   5.000   5.000  1.00  0.00           O
        HETATM    7  H1  HOH B   1       5.500   5.500   5.500  1.00  0.00           H
        HETATM    8  H2  HOH B   1       5.500   4.500   5.500  1.00  0.00           H
        HETATM    9  O   HOH B   2      -5.000  -5.000  -5.000  1.00  0.00           O
        HETATM   10  H1  HOH B   2      -5.500  -5.500  -5.500  1.00  0.00           H
        HETATM   11  H2  HOH B   2      -5.500  -4.500  -5.500  1.00  0.00           H
        HETATM   12 NA   NA  B   3       7.000   0.000   0.000  1.00  0.00          Na
        END
        """)
    frame_path.write_text(pdb)


def _write_incomplete_frame(frame_path: Path) -> None:
    """A buggy frame missing the CB atom — should raise."""
    pdb = textwrap.dedent("""\
        ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
        ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
        ATOM      3  C   ALA A   1       2.300   1.000   0.000  1.00  0.00           C
        ATOM      4  O   ALA A   1       2.000   2.100   0.000  1.00  0.00           O
        END
        """)
    frame_path.write_text(pdb)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSoluteOnlyFrameHappyPath:
    """The [B2] frame writer dumps solute-only PDBs — mask is a no-op."""

    def test_solute_only_frame_passes_through_unchanged(self, tmp_path):
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_solute_only_frame(frames / "frame_000.pdb")
        _write_solute_only_frame(frames / "frame_001.pdb")

        top = _make_solute_topology()
        positions = _collect_frame_positions(
            md, expected_topology=top, solvent_mode="explicit",
        )
        assert len(positions) == 2
        assert positions[0].shape == (5, 3)
        # Order: N, CA, C, O, CB (from the PDB / topology).
        np.testing.assert_allclose(positions[0][0], [0.0, 0.0, 0.0], atol=1e-3)
        np.testing.assert_allclose(positions[0][1], [0.15, 0.0, 0.0], atol=1e-3)
        np.testing.assert_allclose(positions[0][4], [0.19, -0.1, 0.05], atol=1e-3)


class TestSolvatedFrameStrip:
    """A solvated frame must have its waters + ions filtered out."""

    def test_waters_and_ions_dropped_solute_retained(self, tmp_path):
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_solvated_frame(frames / "frame_000.pdb")

        top = _make_solute_topology()
        positions = _collect_frame_positions(
            md, expected_topology=top, solvent_mode="explicit",
        )
        assert len(positions) == 1
        # 5 solute atoms remain; HOH x2 + NA x1 = 7 atoms dropped.
        assert positions[0].shape == (5, 3)
        np.testing.assert_allclose(positions[0][0], [0.0, 0.0, 0.0], atol=1e-3)
        np.testing.assert_allclose(positions[0][4], [0.19, -0.1, 0.05], atol=1e-3)

    def test_strip_preserves_topology_atom_order(self, tmp_path):
        """The output positions array must be indexed by the expected
        topology's atom indices, not by source-PDB line order."""
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_solvated_frame(frames / "frame_000.pdb")

        top = _make_solute_topology()
        positions = _collect_frame_positions(
            md, expected_topology=top, solvent_mode="explicit",
        )
        # CB is the last topology atom (index 4); confirm its position
        # is in slot [4] regardless of where in the PDB it landed.
        expected_cb = [0.19, -0.1, 0.05]  # nm
        np.testing.assert_allclose(positions[0][4], expected_cb, atol=1e-3)


class TestIncompleteFrameRaises:
    """Missing solute atoms should raise — never silently truncate."""

    def test_missing_atom_raises_runtime_error(self, tmp_path):
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_incomplete_frame(frames / "frame_000.pdb")

        top = _make_solute_topology()  # expects 5 atoms; PDB has 4
        with pytest.raises(RuntimeError) as exc_info:
            _collect_frame_positions(
                md, expected_topology=top, solvent_mode="explicit",
            )
        msg = str(exc_info.value)
        assert "frame_000.pdb" in msg
        assert "expected 5" in msg
        assert "matched 4" in msg
        assert "solvent_mode=explicit" in msg


class TestLegacyNoMaskPath:
    """expected_topology=None preserves the pre-[B4] behavior."""

    def test_no_mask_returns_all_positions(self, tmp_path):
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_solvated_frame(frames / "frame_000.pdb")

        positions = _collect_frame_positions(md)  # no expected_topology
        assert len(positions) == 1
        # All 12 atoms are returned (5 solute + 7 solvent).
        assert positions[0].shape == (12, 3)

    def test_no_mask_implicit_mode_is_also_no_op(self, tmp_path):
        md = tmp_path / "md"
        frames = md / "frames"
        frames.mkdir(parents=True)
        _write_solute_only_frame(frames / "frame_000.pdb")
        positions = _collect_frame_positions(md, solvent_mode="implicit")
        assert positions[0].shape == (5, 3)


class TestEmptyTrajectory:
    """No frame files → empty list, no crash."""

    def test_no_frames_returns_empty_list(self, tmp_path):
        md = tmp_path / "md"
        (md / "frames").mkdir(parents=True)
        top = _make_solute_topology()
        positions = _collect_frame_positions(
            md, expected_topology=top, solvent_mode="implicit",
        )
        assert positions == []
