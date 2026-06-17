import numpy as np
import pytest

from backend.app.utils.rmsd import kabsch_rmsd


def test_empty_input_raises():
    with pytest.raises(ValueError, match="empty"):
        kabsch_rmsd(np.empty((0, 3)), np.empty((0, 3)))


def test_mismatched_shape_raises():
    with pytest.raises(ValueError, match="equal shape"):
        kabsch_rmsd(np.zeros((4, 3)), np.zeros((5, 3)))


def test_wrong_dimensionality_raises():
    with pytest.raises(ValueError, match="N,3"):
        kabsch_rmsd(np.zeros(9), np.zeros(9))


def test_identical_coords_zero_rmsd():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert kabsch_rmsd(coords, coords) == pytest.approx(0.0, abs=1e-9)
