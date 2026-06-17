# FILE: backend/app/utils/rmsd.py

import numpy as np


def kabsch_rmsd(ref_coords, mob_coords) -> float:
    """
    Compute RMSD between two sets of coordinates using the Kabsch algorithm.

    ref_coords: (N, 3)
    mob_coords: (N, 3)
    """
    ref = np.array(ref_coords, dtype=float)
    mob = np.array(mob_coords, dtype=float)
    if ref.shape != mob.shape or ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError(
            f"Expected (N,3) arrays of equal shape, got {ref.shape} and {mob.shape}"
        )
    if ref.shape[0] == 0:
        raise ValueError("Cannot compute RMSD on empty coordinate set")

    # Center
    ref_cent = ref.mean(axis=0)
    mob_cent = mob.mean(axis=0)
    ref_c = ref - ref_cent
    mob_c = mob - mob_cent

    # Covariance
    C = mob_c.T @ ref_c
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1.0, 1.0, d]) @ Wt

    # Rotate
    mob_rot = np.dot(mob_c, U)
    diff = ref_c - mob_rot
    rmsd = np.sqrt((diff * diff).sum() / ref.shape[0])
    return float(rmsd)


def kabsch_fit(ref_coords, mob_coords):
    """
    Return the Kabsch transform that aligns mob → ref, plus the post-fit RMSD.

    The caller decides what to do with the transform. Q6b uses this on the
    receptor backbone Cα to recover the rigid-body (R, t) that brings each MD
    frame into the reference receptor frame; that SAME (R, t) is then applied
    to the ligand heavy atoms to compute pose RMSD (pocket displacement)
    without performing a second alignment. The scalar `kabsch_rmsd` above
    cannot do this — it discards the transform.

    The transform satisfies:
        mob_aligned = (mob - mob_centroid) @ R + ref_centroid
                    = mob @ R + t      where   t = ref_centroid - mob_centroid @ R

    Returns
    -------
    R : (3, 3) ndarray
        Proper rotation matrix (det(R) = +1) such that mob @ R aligns the
        mob-centered points onto the ref-centered points.
    t : (3,) ndarray
        Translation vector such that mob @ R + t aligns mob onto ref.
    rmsd : float
        Post-fit RMSD between ref and the aligned mob, in the same units as
        the input coordinates.

    Raises
    ------
    ValueError
        On empty input or shape mismatch (mirrors kabsch_rmsd's contract).
    """
    ref = np.array(ref_coords, dtype=float)
    mob = np.array(mob_coords, dtype=float)
    if ref.shape != mob.shape or ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError(
            f"Expected (N,3) arrays of equal shape, got {ref.shape} and {mob.shape}"
        )
    if ref.shape[0] == 0:
        raise ValueError("Cannot compute RMSD on empty coordinate set")

    ref_cent = ref.mean(axis=0)
    mob_cent = mob.mean(axis=0)
    ref_c = ref - ref_cent
    mob_c = mob - mob_cent

    C = mob_c.T @ ref_c
    V, S, Wt = np.linalg.svd(C)
    # Reflection guard — keep det(R) = +1.
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    t = ref_cent - mob_cent @ R

    mob_aligned = mob_c @ R
    diff = ref_c - mob_aligned
    rmsd = float(np.sqrt((diff * diff).sum() / ref.shape[0]))
    return R, t, rmsd