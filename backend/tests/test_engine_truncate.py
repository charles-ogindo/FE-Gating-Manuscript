"""
Hermetic tests for engine_openmm._truncate_chain_residues_with_oxt_cap.

Synthetic 5-residue mini-protein topology (no SystemGenerator, no sqm)
exercising:
  - the truncated residues are deleted from the modeller
  - the anchor residue gains an OXT atom (natural COO- C-terminus, the
    [B9] design decision after the NME hand-built path proved brittle
    with OpenMM's Modeller.addHydrogens template walker — see the helper's
    docstring)
  - provenance records chain, range, n removed, anchor id, cap atom, etc.
  - geometric placement of OXT is roughly opposite to O across C at the
    same distance (within ~0.05 Å of the |C-O| reference)
  - invalid inputs (bad range, missing anchor, unknown chain) raise
  - multi-chain disambiguation: a request for chain A picks the protein
    chain even when md_receptor_prep has appended a second chain A
    holding only metal ions
"""

from __future__ import annotations

import math

import pytest


# ---------------------------------------------------------------------------
# Synthetic 5-residue chain (ALA-ALA-ALA-ALA-ALA) with realistic backbone
# coordinates so the OXT placement geometry is exercised on something other
# than a degenerate triangle.
# ---------------------------------------------------------------------------

def _build_modeller_5ala():
    import openmm
    from openmm.app import Topology, Element, Modeller
    from openmm import Vec3, unit

    top = Topology()
    chain = top.addChain("A")
    positions = []
    # Place residues along the x-axis at ~3.8 Å spacing (Cα-Cα).
    # Per residue, layout: N, CA, C, O, CB (heavy atoms only — H added later
    # by the engine via addHydrogens against amber14).
    for i in range(1, 6):  # residues 1..5
        res = top.addResidue("ALA", chain, id=str(i))
        base = i * 0.38  # 0.38 nm (~3.8 Å)
        n_pos  = Vec3(base + 0.00, 0.00, 0.00)
        ca_pos = Vec3(base + 0.14, 0.05, 0.00)
        c_pos  = Vec3(base + 0.28, 0.00, 0.00)
        o_pos  = Vec3(base + 0.28, -0.12, 0.05)
        cb_pos = Vec3(base + 0.14, 0.05, 0.15)
        n_atom  = top.addAtom("N",  Element.getBySymbol("N"), res)
        ca_atom = top.addAtom("CA", Element.getBySymbol("C"), res)
        c_atom  = top.addAtom("C",  Element.getBySymbol("C"), res)
        o_atom  = top.addAtom("O",  Element.getBySymbol("O"), res)
        cb_atom = top.addAtom("CB", Element.getBySymbol("C"), res)
        top.addBond(n_atom, ca_atom)
        top.addBond(ca_atom, c_atom)
        top.addBond(c_atom, o_atom)
        top.addBond(ca_atom, cb_atom)
        positions.extend([n_pos, ca_pos, c_pos, o_pos, cb_pos])
    # Apply units to the LIST, not to each Vec3 — Quantity-of-list expects
    # raw Vec3 elements; wrapping each element with units double-wraps and
    # blows the OpenMM Quantity __getitem__ assertion.
    pos_quantity = positions * unit.nanometer
    return Modeller(top, pos_quantity)


def _call_truncate(modeller, ranges):
    """Wrap the engine helper; pass mm/app/unit so the function can build
    Vec3/Quantity objects + the cap Topology without re-importing openmm."""
    import openmm as mm
    from openmm import app, unit
    from backend.app.md.engine_openmm import _truncate_chain_residues_with_oxt_cap

    msgs = []
    def say(m): msgs.append(m)

    prov = _truncate_chain_residues_with_oxt_cap(
        modeller, ranges, mm=mm, app=app, unit=unit, say=say,
    )
    return prov, msgs


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestTruncateHappyPath:
    """Delete residues 4-5, cap anchor residue 3 with OXT."""

    def test_deletes_target_residues_and_adds_oxt(self):
        modeller = _build_modeller_5ala()
        before_count = modeller.topology.getNumResidues()
        assert before_count == 5

        prov, msgs = _call_truncate(
            modeller,
            [{"chain": "A", "from_resseq": 4, "to_resseq": 5}],
        )

        # 5 ALA - 2 truncated = 3 ALA residues remain. OXT is added to
        # residue 3 — NOT as a new residue, but as a new atom on the anchor.
        after_residues = list(modeller.topology.residues())
        residue_names = [r.name for r in after_residues]
        assert residue_names.count("ALA") == 3
        # The OXT atom should appear inside the anchor residue (residue 3).
        anchor = None
        for r in after_residues:
            if r.id == "3" and r.name == "ALA":
                anchor = r
                break
        assert anchor is not None, "anchor residue 3 missing after truncation"
        atom_names = sorted(a.name for a in anchor.atoms())
        assert "OXT" in atom_names, f"OXT not found on anchor; atoms={atom_names}"

    def test_provenance_records_truncation_details(self):
        modeller = _build_modeller_5ala()
        prov, _ = _call_truncate(
            modeller,
            [{"chain": "A", "from_resseq": 4, "to_resseq": 5}],
        )
        assert len(prov) == 1
        p = prov[0]
        assert p["chain"] == "A"
        assert p["from_resseq"] == 4
        assert p["to_resseq"] == 5
        assert p["n_residues_removed"] == 2
        assert p["cap"] == "OXT"
        assert p["cap_atom_name"] == "OXT"
        assert p["anchor_residue_id"] == "3"
        assert p["anchor_residue_name"] == "ALA"
        assert p["source"] == "engine_openmm"
        assert "anchor_C_xyz_nm" in p
        assert "oxt_xyz_nm" in p
        assert len(p["anchor_C_xyz_nm"]) == 3
        assert len(p["oxt_xyz_nm"]) == 3

    def test_oxt_geometry_collinear_opposite_of_O_across_C(self):
        modeller = _build_modeller_5ala()
        prov, _ = _call_truncate(
            modeller,
            [{"chain": "A", "from_resseq": 4, "to_resseq": 5}],
        )
        # Pull C and OXT positions in nm from the reloaded modeller.
        from openmm import unit
        positions = modeller.positions
        anchor_c_xyz = None
        anchor_o_xyz = None
        oxt_xyz = None
        for r in modeller.topology.residues():
            if r.id == "3" and r.name == "ALA":
                for a in r.atoms():
                    p = positions[a.index].value_in_unit(unit.nanometer)
                    if a.name == "C":
                        anchor_c_xyz = (p.x, p.y, p.z)
                    elif a.name == "O":
                        anchor_o_xyz = (p.x, p.y, p.z)
                    elif a.name == "OXT":
                        oxt_xyz = (p.x, p.y, p.z)
        assert anchor_c_xyz is not None
        assert anchor_o_xyz is not None
        assert oxt_xyz is not None

        cx, cy, cz = anchor_c_xyz
        ox, oy, oz = anchor_o_xyz
        xx, xy, xz = oxt_xyz

        # |C - OXT| should equal |C - O| (collinear-opposite placement)
        d_co = math.sqrt((cx - ox) ** 2 + (cy - oy) ** 2 + (cz - oz) ** 2)
        d_cx = math.sqrt((cx - xx) ** 2 + (cy - xy) ** 2 + (cz - xz) ** 2)
        # PDB writer rounds to 3 decimal places in Å (0.0001 nm precision);
        # tolerate ±0.001 nm for the round-trip + the geometry test.
        assert abs(d_co - d_cx) < 0.001, (
            f"|C-O|={d_co:.4f} nm vs |C-OXT|={d_cx:.4f} nm "
            "(should be equal for opposite-collinear placement)"
        )
        # OXT direction should be opposite to O direction relative to C.
        # Check via dot product on the C-O vs C-OXT vectors: should be ≈ -1
        # times |C-O|² (anti-parallel of equal length).
        co_vec = (ox - cx, oy - cy, oz - cz)
        cx_vec = (xx - cx, xy - cy, xz - cz)
        dot = sum(co_vec[i] * cx_vec[i] for i in range(3))
        # dot = -|C-O|² for opposite-collinear
        expected_dot = -(d_co * d_co)
        assert abs(dot - expected_dot) < 0.001, (
            f"C-O · C-OXT = {dot:.4f}, expected ≈ {expected_dot:.4f}"
        )


class TestTruncateWithMultipleSameIdChains:
    """Two chain objects sharing id='A' — the protein chain and an ion-only
    chain that md_receptor_prep appends. Truncation must select the chain
    that holds the requested residue numbers, NOT the first chain with
    matching id."""

    def _build_modeller_with_ion_chain(self):
        """5-ALA chain A + a second chain A with one MG residue at 501."""
        from openmm.app import Topology, Element, Modeller
        from openmm import Vec3, unit
        top = Topology()
        # Protein chain A: 5 ALAs at residues 1..5
        chain1 = top.addChain("A")
        positions = []
        for i in range(1, 6):
            res = top.addResidue("ALA", chain1, id=str(i))
            base = i * 0.38
            n_atom  = top.addAtom("N",  Element.getBySymbol("N"), res)
            ca_atom = top.addAtom("CA", Element.getBySymbol("C"), res)
            c_atom  = top.addAtom("C",  Element.getBySymbol("C"), res)
            o_atom  = top.addAtom("O",  Element.getBySymbol("O"), res)
            cb_atom = top.addAtom("CB", Element.getBySymbol("C"), res)
            top.addBond(n_atom, ca_atom)
            top.addBond(ca_atom, c_atom)
            top.addBond(c_atom, o_atom)
            top.addBond(ca_atom, cb_atom)
            positions.extend([
                Vec3(base + 0.00, 0.00, 0.00),
                Vec3(base + 0.14, 0.05, 0.00),
                Vec3(base + 0.28, 0.00, 0.00),
                Vec3(base + 0.28, -0.12, 0.05),
                Vec3(base + 0.14, 0.05, 0.15),
            ])
        # Second chain A: an MG ion at residue 501 — mirrors how
        # md_receptor_prep appends the ion modeller.
        chain2 = top.addChain("A")
        ion_res = top.addResidue("MG", chain2, id="501")
        top.addAtom("MG", Element.getBySymbol("Mg"), ion_res)
        positions.append(Vec3(10.0, 10.0, 10.0))
        pos_quantity = positions * unit.nanometer
        return Modeller(top, pos_quantity)

    def test_truncation_picks_protein_chain_over_ion_chain(self):
        """Range 4-5 on chain A should target the protein chain, NOT the
        ion chain. Pre-fix, this raised 'no residues match'."""
        modeller = self._build_modeller_with_ion_chain()
        before_protein_residues = sum(
            1 for c in modeller.topology.chains() if c.id == "A"
            for r in c.residues() if r.name == "ALA"
        )
        assert before_protein_residues == 5

        prov, _ = _call_truncate(
            modeller,
            [{"chain": "A", "from_resseq": 4, "to_resseq": 5}],
        )
        assert len(prov) == 1
        # 3 ALAs remain (residues 1-3), MG remains (ion chain).
        residue_summary = []
        for c in modeller.topology.chains():
            for r in c.residues():
                residue_summary.append((c.id, r.name, r.id))
        ala_remaining = [r for r in residue_summary if r[1] == "ALA"]
        mg_kept = [r for r in residue_summary if r[1] == "MG"]
        assert len(ala_remaining) == 3
        assert len(mg_kept) == 1
        assert mg_kept[0] == ("A", "MG", "501")
        # OXT should be on the protein chain's anchor (residue 3)
        anchor = None
        for c in modeller.topology.chains():
            for r in c.residues():
                if r.id == "3" and r.name == "ALA":
                    anchor = r
        assert anchor is not None
        atom_names = sorted(a.name for a in anchor.atoms())
        assert "OXT" in atom_names

    def test_error_message_lists_present_chain_ranges_on_miss(self):
        """If the request misses every chain-A range, the error must list
        each candidate chain's residue range so the operator can debug."""
        modeller = self._build_modeller_with_ion_chain()
        with pytest.raises(ValueError) as exc:
            _call_truncate(
                modeller,
                [{"chain": "A", "from_resseq": 100, "to_resseq": 200}],
            )
        msg = str(exc.value)
        assert "no residues in chain A match resseq 100-200" in msg
        assert "1-5" in msg or "501-501" in msg


class TestTruncateErrorPaths:
    """Bad inputs must fail loudly, never silently."""

    def test_inverted_range_raises(self):
        modeller = _build_modeller_5ala()
        with pytest.raises(ValueError) as exc:
            _call_truncate(
                modeller,
                [{"chain": "A", "from_resseq": 5, "to_resseq": 3}],
            )
        assert "from_resseq=5 > to_resseq=3" in str(exc.value)

    def test_unknown_chain_raises(self):
        modeller = _build_modeller_5ala()
        with pytest.raises(ValueError) as exc:
            _call_truncate(
                modeller,
                [{"chain": "Z", "from_resseq": 4, "to_resseq": 5}],
            )
        assert "chain 'Z' not present" in str(exc.value)

    def test_missing_anchor_raises(self):
        """If from_resseq == first_residue_in_chain, there is no anchor."""
        modeller = _build_modeller_5ala()
        with pytest.raises(ValueError) as exc:
            _call_truncate(
                modeller,
                [{"chain": "A", "from_resseq": 1, "to_resseq": 5}],
            )
        msg = str(exc.value)
        assert "anchor" in msg

    def test_empty_range_no_op(self):
        """Empty ranges list returns [] without mutating the modeller."""
        modeller = _build_modeller_5ala()
        before = modeller.topology.getNumResidues()
        prov, _ = _call_truncate(modeller, [])
        assert prov == []
        assert modeller.topology.getNumResidues() == before


class TestTruncateLogsAreInformative:
    """The engine's `say` log must surface the truncation extent + cap atom."""

    def test_log_mentions_chain_range_and_oxt_cap(self):
        modeller = _build_modeller_5ala()
        _, msgs = _call_truncate(
            modeller,
            [{"chain": "A", "from_resseq": 4, "to_resseq": 5}],
        )
        # At least one message about the truncation
        combined = "\n".join(msgs)
        assert "Truncating chain A residues 4-5" in combined
        assert "2 residue(s)" in combined
        assert "OXT" in combined
        assert "COO-" in combined or "C-terminus" in combined
