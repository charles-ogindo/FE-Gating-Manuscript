"""Tests for the single-trajectory MM-GBSA estimator + FE artifact runner.

Two scopes:
  1. `mmgbsa.estimate_mmgbsa` correctness — synthetic systems, fast.
     Includes the BIT-FOR-BIT refactor regression (fast-path reuse-Context
     vs slow-path build-Context-per-call must produce identical numerics).
  2. `mmgbsa_runner.compute_md_fe` artifact-shape — gate-blocked-but-
     computed path with a synthetic fixture (no real MD frames needed).

The real-data taxol smoke is NOT a unit test (it's `sqm`-bound on a
cold cache and runs for ~hours). The bit-for-bit refactor regression
below exercises the same OpenMM force classes the taxol path uses
(HarmonicBond + NonbondedForce + GBSAOBCForce + force-group dispatch)
on a small synthetic system in seconds, which is strictly stronger
than the aggregate-within-tolerance check the slow-path smoke would
provide.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest


# ---------------------------------------------------------------------
# Helpers for tiny OpenMM fixtures
# ---------------------------------------------------------------------

def _make_topology(n_atoms, resnames):
    import openmm.app as app
    t = app.Topology()
    chain = t.addChain()
    for i in range(n_atoms):
        res = t.addResidue(resnames[i], chain)
        t.addAtom("C", app.element.carbon, res)
    return t


def _make_system_with_obc2(n_atoms, with_bonds=True):
    """Real force mix matching the taxol path: HarmonicBondForce +
    NonbondedForce + GBSAOBCForce."""
    import openmm as mm
    from openmm import unit
    s = mm.System()
    for _ in range(n_atoms):
        s.addParticle(12.0)
    if with_bonds and n_atoms >= 2:
        bond = mm.HarmonicBondForce()
        for i in range(n_atoms - 1):
            bond.addBond(
                i, i + 1,
                0.15 * unit.nanometer,
                50000.0 * unit.kilojoule_per_mole / unit.nanometer**2,
            )
        s.addForce(bond)
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(1.0 * unit.nanometer)
    for i in range(n_atoms):
        nb.addParticle(
            0.1 if i % 2 == 0 else -0.1,
            0.35 * unit.nanometer,
            0.5 * unit.kilojoule_per_mole,
        )
    s.addForce(nb)
    gb = mm.GBSAOBCForce()
    gb.setNonbondedMethod(mm.GBSAOBCForce.CutoffNonPeriodic)
    gb.setCutoffDistance(1.0 * unit.nanometer)
    for i in range(n_atoms):
        gb.addParticle(0.1 if i % 2 == 0 else -0.1, 0.17, 0.85)
    s.addForce(gb)
    return s


# =====================================================================
# Test 1: refactor regression — fast vs slow path bit-for-bit
# =====================================================================

class TestRefactorRegression:
    """The mmgbsa.py refactor swapped one-Context-per-call for one-Context
    per-Subsystem (303 → 3 builds on the taxol smoke). This test asserts
    the swap is numerically transparent — identical (system, positions,
    integrator config) must yield identical (total + per-bucket) energies
    regardless of Context lifecycle. OpenMM CPU is deterministic.

    Same force mix the real taxol estimator uses, so a pass here is a
    valid proxy for a pass on the full trajectory."""

    def _slow_path_eval(self, system, topology, pos_nm_array):
        """Emulate the pre-refactor slow path: one Context per call."""
        import openmm as mm
        from openmm import app, unit
        from backend.app.free_energy.mmgbsa import _FORCE_GROUP_BUCKET

        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1.0 / unit.picosecond, 2 * unit.femtosecond,
        )
        platform = mm.Platform.getPlatformByName("CPU")
        sim = app.Simulation(topology, system, integrator, platform)
        sim.context.setPositions(pos_nm_array * unit.nanometer)
        e_total = (
            sim.context.getState(getEnergy=True)
            .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        )
        buckets = {"bonded": 0.0, "nonbonded": 0.0, "solvation": 0.0, "other": 0.0}
        for f in system.getForces():
            g = f.getForceGroup()
            b = _FORCE_GROUP_BUCKET.get(f.__class__.__name__, "other")
            buckets[b] += (
                sim.context.getState(getEnergy=True, groups={g})
                .getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            )
        return e_total, buckets

    def test_fast_path_matches_slow_path_bit_for_bit(self):
        from openmm import unit
        from backend.app.free_energy.mmgbsa import (
            _assign_force_groups, _build_context, _evaluate_context,
        )

        N_REC, N_LIG = 6, 2
        N_CPX = N_REC + N_LIG
        rec_sys = _make_system_with_obc2(N_REC)
        lig_sys = _make_system_with_obc2(N_LIG)
        cpx_sys = _make_system_with_obc2(N_CPX)
        rec_top = _make_topology(N_REC, ["RCP"] * N_REC)
        lig_top = _make_topology(N_LIG, ["LIG"] * N_LIG)
        cpx_top = _make_topology(N_CPX, ["RCP"] * N_REC + ["LIG"] * N_LIG)

        # Force-group assignment BEFORE Context construction (same as
        # estimate_mmgbsa does it).
        for s in (cpx_sys, rec_sys, lig_sys):
            _assign_force_groups(s)

        cpx_ctx = _build_context(cpx_sys, cpx_top)
        rec_ctx = _build_context(rec_sys, rec_top)
        lig_ctx = _build_context(lig_sys, lig_top)

        receptor_idx = list(range(N_REC))
        ligand_idx = list(range(N_REC, N_CPX))

        # Deterministic 10 frames at fixed seed.
        rng = np.random.default_rng(42)
        frames = [rng.uniform(-0.5, 0.5, size=(N_CPX, 3)) for _ in range(10)]

        max_diff = 0.0
        for pos in frames:
            cpx_pos = pos
            rec_pos = pos[receptor_idx]
            lig_pos = pos[ligand_idx]

            e_cpx_s, b_cpx_s = self._slow_path_eval(cpx_sys, cpx_top, cpx_pos)
            e_rec_s, b_rec_s = self._slow_path_eval(rec_sys, rec_top, rec_pos)
            e_lig_s, b_lig_s = self._slow_path_eval(lig_sys, lig_top, lig_pos)

            e_cpx_f, b_cpx_f = _evaluate_context(cpx_ctx, cpx_pos * unit.nanometer)
            e_rec_f, b_rec_f = _evaluate_context(rec_ctx, rec_pos * unit.nanometer)
            e_lig_f, b_lig_f = _evaluate_context(lig_ctx, lig_pos * unit.nanometer)

            # Total energies match.
            assert e_cpx_s == e_cpx_f, f"complex total drift: {e_cpx_s} != {e_cpx_f}"
            assert e_rec_s == e_rec_f, f"receptor total drift: {e_rec_s} != {e_rec_f}"
            assert e_lig_s == e_lig_f, f"ligand total drift: {e_lig_s} != {e_lig_f}"

            # Per-bucket energies match.
            for bucket in ("bonded", "nonbonded", "solvation"):
                assert b_cpx_s[bucket] == b_cpx_f[bucket], \
                    f"complex {bucket} drift: {b_cpx_s[bucket]} != {b_cpx_f[bucket]}"

            dg_slow = e_cpx_s - e_rec_s - e_lig_s
            dg_fast = e_cpx_f - e_rec_f - e_lig_f
            max_diff = max(max_diff, abs(dg_slow - dg_fast))

        assert max_diff == 0.0, f"ΔG drift > 0 over 10 frames: max |Δ|={max_diff}"


# =====================================================================
# Test 2: synthetic two-frame estimator sanity
# =====================================================================

class TestEstimatorSanity:
    """Build a deterministic 3-atom 'complex' (Ar2 receptor + Ar1 ligand)
    with a single NonbondedForce (LJ, no charges, no GB). Verify the
    estimator returns ΔG=0 at infinite separation and -ε at the LJ
    contact. Catches contract regressions on the estimator surface
    (return shape, aggregation, kJ→kcal conversion)."""

    def _make_lj_system(self, n_atoms):
        import openmm as mm
        from openmm import unit
        s = mm.System()
        for _ in range(n_atoms):
            s.addParticle(40.0)  # Ar
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
        nb.setCutoffDistance(2.0 * unit.nanometer)
        for _ in range(n_atoms):
            nb.addParticle(
                0.0, 0.34 * unit.nanometer, 0.997 * unit.kilojoule_per_mole,
            )
        s.addForce(nb)
        return s

    def test_far_apart_and_lj_contact(self):
        from backend.app.free_energy.mmgbsa import estimate_mmgbsa

        cpx = self._make_lj_system(3)
        rec = self._make_lj_system(2)
        lig = self._make_lj_system(1)
        ct = _make_topology(3, ["RCP", "RCP", "LIG"])
        rt = _make_topology(2, ["RCP", "RCP"])
        lt = _make_topology(1, ["LIG"])

        sig = 0.34
        rmin = sig * (2 ** (1 / 6))
        far = np.array([[0, 0, 0], [sig * 1.5, 0, 0], [10, 10, 10]])
        close = np.array([[0, 0, 0], [sig * 1.5, 0, 0], [sig * 1.5 + rmin, 0, 0]])

        result = estimate_mmgbsa(
            complex_system=cpx, complex_topology=ct,
            receptor_system=rec, receptor_topology=rt,
            ligand_system=lig, ligand_topology=lt,
            frame_positions=[far, close],
            receptor_idx=[0, 1], ligand_idx=[2],
            gate_can_run=True, gate_reason=None,
            method_meta={"test": "lj"},
        )
        assert result.n_frames_used == 2
        assert result.n_frames_skipped == 0
        # Far frame: ΔG ~ 0 (no interaction beyond cutoff).
        assert abs(result.per_frame[0].delta_g_total) < 1e-3
        # Close frame: ΔG ≈ -ε in kcal/mol.
        expected = -0.997 / 4.184
        assert abs(result.per_frame[1].delta_g_total - expected) < 1e-2
        # Result is preliminary=False, sampling_adequate=True given gate True.
        assert result.preliminary is False
        assert result.sampling_adequate is True
        assert result.gate_reason is None


# =====================================================================
# Test 3: gate-blocked-but-computed — the integrity-of-gating contract
# =====================================================================

class TestGateBlockedStillComputes:
    """If the FE gate blocks (e.g., taxol's INSUFFICIENT_FRAMES), the
    estimator MUST still compute on the available frames AND flag the
    artifact preliminary=True / sampling_adequate=False / gate_reason=<str>.
    The threshold (MIN_TOTAL_FRAMES) must NEVER be lowered to make a
    real-data fixture "pass" — that threshold is the integrity of the
    gating discipline."""

    def test_blocked_but_computed_marks_preliminary(self):
        from backend.app.free_energy.mmgbsa import estimate_mmgbsa

        sys_cpx = _make_system_with_obc2(3)
        sys_rec = _make_system_with_obc2(2)
        sys_lig = _make_system_with_obc2(1)
        top_cpx = _make_topology(3, ["RCP", "RCP", "LIG"])
        top_rec = _make_topology(2, ["RCP", "RCP"])
        top_lig = _make_topology(1, ["LIG"])

        rng = np.random.default_rng(7)
        frames = [rng.uniform(-0.5, 0.5, size=(3, 3)) for _ in range(5)]

        gate_reason = (
            "MD produced 5 frames; need at least 200 total. Re-run MD with "
            "longer production (or smaller snapshot interval)."
        )
        result = estimate_mmgbsa(
            complex_system=sys_cpx, complex_topology=top_cpx,
            receptor_system=sys_rec, receptor_topology=top_rec,
            ligand_system=sys_lig, ligand_topology=top_lig,
            frame_positions=frames,
            receptor_idx=[0, 1], ligand_idx=[2],
            gate_can_run=False, gate_reason=gate_reason,
            method_meta={"test": "gate-blocked"},
        )
        # Still computed.
        assert result.n_frames_used == 5
        # Honestly flagged.
        assert result.sampling_adequate is False
        assert result.preliminary is True
        assert result.gate_reason == gate_reason
        # Numbers are present (not NaN) even though preliminary.
        import math
        assert math.isfinite(result.delta_g_mean)


# =====================================================================
# Test 4: runner artifact shape — gate-blocked taxol-shaped fixture
# =====================================================================

class TestArtifactShape:
    """compute_md_fe's output is the `free_energy` block that gets merged
    into MD summary.json. Smoke the shape end-to-end against a fully
    synthetic fixture (no real MD frames or sqm). Uses monkeypatches to
    swap the heavy System construction for tiny synthetic Systems so the
    test runs in seconds."""

    def test_artifact_block_shape_and_preliminary_flagging(
        self, tmp_path, monkeypatch,
    ):
        # Stand up a tiny MD job dir: summary.json + 5 frames + receptor
        # docking dir with required metadata.
        md_id = "11111111-2222-3333-4444-555555555555"
        docking_id = "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
        ligand_name = "lig_000000"
        jobs = tmp_path / "jobs"

        # Docking job — minimum to satisfy compute_md_fe's reads.
        dock_dir = jobs / docking_id
        (dock_dir / "docking").mkdir(parents=True)
        (dock_dir / "metadata.json").write_text(json.dumps({
            "receptor_prep": {"parameterized_cofactors": []},
        }))
        # Note: dock_dir / "docking/viewer/receptor.pdb" is NOT created;
        # we monkeypatch prepare_mmgbsa_systems so build path is skipped.

        # MD job dir + 5 frames + summary.json with the bits compute_md_fe reads.
        md_dir = jobs / md_id / "md"
        (md_dir / "frames").mkdir(parents=True)
        (md_dir / "summary.json").write_text(json.dumps({
            "schema_version": "1.1.0", "md_job_id": md_id,
            "docking_job_id": docking_id, "ligand": ligand_name, "pose_rank": 0,
            "status": "completed", "verdict": "stable",
            "engine": {"kind": "openmm_full"},
            "settings": {"snapshot_every_ps": 5.0},
            "n_frames": 5,
            "metrics": {
                "rmsd_ligand_pose_max_a": 1.5,
                "rmsd_backbone_final_a": 1.5,
            },
            "free_energy": {"status": "planned"},  # the pre-FE marker
        }))
        # Frame PDBs — minimum that parses; 3 atoms each.
        for i in range(5):
            (md_dir / "frames" / f"frame_{i:03d}.pdb").write_text(
                "HEADER  test\n"
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
                "ATOM      2  CA  ALA A   2       0.300   0.000   0.000  1.00  0.00           C\n"
                "ATOM      3  C   LIG B   1       0.150   0.150   0.300  1.00  0.00           C\n"
                "END\n"
            )

        # Repoint JOBS_DIR + the gating module's JOBS_DIR so validate()
        # finds our synthetic summary. Mirrors the test_md_report.py pattern.
        from backend.app.core import config as cfg
        from backend.app.free_energy import gating as gating_mod
        from backend.app.free_energy import mmgbsa_runner as runner_mod
        monkeypatch.setattr(cfg, "JOBS_DIR", jobs)
        monkeypatch.setattr(gating_mod, "JOBS_DIR", jobs)
        monkeypatch.setattr(runner_mod, "JOBS_DIR", jobs)

        # Replace the heavy System-construction path with tiny synthetic
        # Systems matching the 3-atom frame topology. Skips sqm + openff.
        cpx_sys = _make_system_with_obc2(3)
        rec_sys = _make_system_with_obc2(2)
        lig_sys = _make_system_with_obc2(1)
        cpx_top = _make_topology(3, ["ALA", "ALA", "LIG"])
        rec_top = _make_topology(2, ["ALA", "ALA"])
        lig_top = _make_topology(1, ["LIG"])

        def fake_prepare(md_id_arg, **kw):
            return (cpx_sys, rec_sys, lig_sys,
                    cpx_top, rec_top, lig_top,
                    [0, 1], [2])
        monkeypatch.setattr(runner_mod, "prepare_mmgbsa_systems", fake_prepare)

        from backend.app.free_energy.mmgbsa_runner import compute_md_fe
        fe = compute_md_fe(md_id, jobs_dir=jobs)

        # --- Shape checks (the self-describing-artifact contract) ---
        assert fe["status"] == "completed"
        assert fe["schema_version"] == "1.0.0"
        # method block fully populated
        assert fe["method"]["name"] == "single-trajectory MM-GBSA"
        assert fe["method"]["implicit_solvent"] == "OBC2"
        assert "gaff-2.11" in fe["method"]["small_molecule_forcefield"]
        assert "omitted" in fe["method"]["configurational_entropy"]
        # criteria.gate carries the gate output verbatim
        assert "gate" in fe["criteria"]
        # Gate must be BLOCKED on this fixture (5 frames < 200) — that's
        # the load-bearing integrity check.
        assert fe["criteria"]["gate"]["can_run"] is False
        assert fe["criteria"]["gate"]["reasons"][0]["key"] == "INSUFFICIENT_FRAMES"
        # provenance has the keys methods sections care about
        for key in ("computed_at_utc", "n_frames_used", "n_frames_skipped",
                    "equilibration_discard_frames", "wall_seconds",
                    "estimator_module"):
            assert key in fe["provenance"], f"missing provenance.{key}"
        # result is preliminary because the gate blocked
        assert fe["result"]["preliminary"] is True
        assert fe["result"]["sampling_adequate"] is False
        assert fe["result"]["gate_reason"] is not None
        assert "200" in fe["result"]["gate_reason"]
        # Numerics are present and finite (computed-despite-block).
        import math
        assert math.isfinite(fe["result"]["delta_g_mean_kcal_per_mol"])
        assert "bonded" in fe["result"]["components_mean_kcal_per_mol"]
        assert "nonbonded" in fe["result"]["components_mean_kcal_per_mol"]
        assert "solvation" in fe["result"]["components_mean_kcal_per_mol"]
        # Per-frame series populated, one entry per frame
        assert len(fe["per_frame"]) == 5
