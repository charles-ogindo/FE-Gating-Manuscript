"""Q6b — pose RMSD via receptor-frame superposition.

Synthetic regression check that internal RMSD (ligand-on-ligand Kabsch)
and pose RMSD (receptor-frame ligand displacement) capture different
quantities:

  - A rigidly-translated ligand with NO conformational change should have
    internal RMSD ≈ 0 (Kabsch removes the translation) AND pose RMSD ≈ the
    translation magnitude (Kabsch on the fixed receptor recovers identity,
    so the ligand-translation survives).
  - A receptor that translates with the ligand should have pose RMSD ≈ 0
    (the alignment cancels out the joint motion) AND internal RMSD ≈ 0.

This is the load-bearing invariant for the Q6b verdict shift: pre-Q6b
"stable" labels keyed on internal RMSD missed pose drift; post-Q6b uses
the pose metric so pocket displacement actually downgrades the verdict.
"""

from __future__ import annotations

import numpy as np

from backend.app.md.analyze import Snapshot, compute_rmsd_series
from backend.app.utils.rmsd import kabsch_fit


# ---------------------------------------------------------------------
# Snapshot fixtures
# ---------------------------------------------------------------------

def _ca_atom(resseq: int, xyz: tuple[float, float, float]) -> dict:
    return {
        "chain": "A",
        "resseq": resseq,
        "resname": "ALA",
        "name": "CA",
        "element": "C",
        "xyz": xyz,
    }


def _lig_atom(name: str, xyz: tuple[float, float, float]) -> dict:
    return {
        "chain": "X",
        "resseq": 999,
        "resname": "LIG",
        "name": name,
        "element": "C",
        "xyz": xyz,
    }


def _receptor_5ca(offset: tuple[float, float, float] = (0, 0, 0)) -> list[dict]:
    """Five backbone Cα atoms scattered in a non-degenerate (non-collinear)
    pattern; enough to fix the rotation in Kabsch_fit."""
    template = [
        (0.0, 0.0, 0.0),
        (3.8, 0.0, 0.0),
        (5.7, 3.3, 0.0),
        (3.8, 6.6, 0.0),
        (0.0, 6.6, 0.0),
    ]
    ox, oy, oz = offset
    return [
        _ca_atom(i + 1, (x + ox, y + oy, z + oz))
        for i, (x, y, z) in enumerate(template)
    ]


def _ligand_4heavy(offset: tuple[float, float, float] = (0, 0, 0)) -> list[dict]:
    """Four heavy atoms ~5 Å away from the receptor 'pocket' so any
    translation we apply is visible in pose RMSD."""
    template = [
        (10.0, 3.0, 0.0),
        (11.4, 3.0, 0.0),
        (10.7, 4.2, 0.0),
        (10.7, 1.8, 0.0),
    ]
    ox, oy, oz = offset
    return [
        _lig_atom(f"C{i+1}", (x + ox, y + oy, z + oz))
        for i, (x, y, z) in enumerate(template)
    ]


# =====================================================================
# Test 1 — rigid ligand translation: internal ≈ 0, pose ≈ |Δ|
# =====================================================================

class TestRigidLigandTranslation:
    def test_ligand_translated_3a_internal_zero_pose_three(self):
        # Frame 0: ref. Frame 1: same receptor, ligand translated +3 Å along x.
        frame0 = Snapshot(
            receptor_atoms=_receptor_5ca(),
            ligand_atoms=_ligand_4heavy(),
        )
        frame1 = Snapshot(
            receptor_atoms=_receptor_5ca(),
            ligand_atoms=_ligand_4heavy(offset=(3.0, 0.0, 0.0)),
        )
        series = compute_rmsd_series([frame0, frame1], [0.0, 1.0])

        # Frame 0 self-comparison: all metrics ~0.
        t0, bb0, lig_int0, lig_pose0 = series[0]
        assert bb0 < 1e-9
        assert lig_int0 < 1e-9
        assert lig_pose0 < 1e-9

        # Frame 1: receptor unchanged → bb RMSD ~0; ligand rigid-translated.
        t1, bb1, lig_int1, lig_pose1 = series[1]
        assert bb1 < 1e-9, f"backbone moved: {bb1}"
        # Internal RMSD removes the translation → ~0.
        assert lig_int1 < 1e-6, (
            f"internal RMSD should be ~0 on a rigid translation, got {lig_int1}"
        )
        # Pose RMSD preserves the translation magnitude → ~3.0 Å.
        assert abs(lig_pose1 - 3.0) < 1e-6, (
            f"pose RMSD should be the translation magnitude ~3.0 Å, got {lig_pose1}"
        )

    def test_joint_translation_both_metrics_zero(self):
        # Frame 1: receptor AND ligand translated together by +5 Å along z.
        # The receptor superposition cancels the joint motion, so pose RMSD ~0.
        # Internal RMSD already ignores translation, so ~0 too.
        frame0 = Snapshot(
            receptor_atoms=_receptor_5ca(),
            ligand_atoms=_ligand_4heavy(),
        )
        frame1 = Snapshot(
            receptor_atoms=_receptor_5ca(offset=(0.0, 0.0, 5.0)),
            ligand_atoms=_ligand_4heavy(offset=(0.0, 0.0, 5.0)),
        )
        series = compute_rmsd_series([frame0, frame1], [0.0, 1.0])
        t1, bb1, lig_int1, lig_pose1 = series[1]
        # Backbone Kabsch absorbs the translation → ~0.
        assert bb1 < 1e-9
        assert lig_int1 < 1e-6
        assert lig_pose1 < 1e-6, (
            f"joint translation should cancel: pose RMSD got {lig_pose1}"
        )


# =====================================================================
# Test 2 — kabsch_fit transform sanity (the helper Q6b adds)
# =====================================================================

class TestKabschFit:
    def test_identity_returns_identity_transform(self):
        coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        R, t, rmsd = kabsch_fit(coords, coords)
        assert np.allclose(R, np.eye(3), atol=1e-9)
        assert np.allclose(t, np.zeros(3), atol=1e-9)
        assert rmsd < 1e-9

    def test_pure_translation_recovered(self):
        ref = np.array([[0, 0, 0], [3.8, 0, 0], [5.7, 3.3, 0], [0, 6.6, 0]], dtype=float)
        # mob is ref translated by (2, -1, 0.5).
        delta = np.array([2.0, -1.0, 0.5])
        mob = ref + delta
        R, t, rmsd = kabsch_fit(ref, mob)
        # R should be identity (no rotation); applying mob @ R + t should
        # recover ref exactly.
        assert np.allclose(R, np.eye(3), atol=1e-9)
        recovered = mob @ R + t
        assert np.allclose(recovered, ref, atol=1e-9)
        assert rmsd < 1e-9

    def test_pure_rotation_recovered(self):
        # 90° rotation around z.
        theta = np.pi / 2
        c, s = np.cos(theta), np.sin(theta)
        R_true = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)
        ref = np.array([[0, 0, 0], [3.8, 0, 0], [5.7, 3.3, 0], [0, 6.6, 1.2]], dtype=float)
        # Note: applying R as (ref @ R_true.T) produces ref rotated by R_true.
        mob = ref @ R_true.T
        R, t, rmsd = kabsch_fit(ref, mob)
        # Applying mob @ R + t should recover ref.
        recovered = mob @ R + t
        assert np.allclose(recovered, ref, atol=1e-9)
        assert rmsd < 1e-9
