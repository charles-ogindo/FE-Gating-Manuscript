"""
Molecular-dynamics validation stage.

This package is a downstream refinement step for docking: given a docked pose
(receptor + ligand PDB pair from jobs/<docking_id>/docking/viewer/...), it runs
a short stability simulation and writes the standard MD artifacts under
jobs/<md_id>/md/.

Engines are lazy-imported by `job.run_md_job` so this package can always be
imported by the API layer even when no MD engine is installed.
"""

# Q6b (2026-06-04): 1.0.0 → 1.1.0 adds rmsd_ligand_pose_{final,max}_a as the
# primary ligand metrics (receptor-frame pose displacement); renames the
# legacy ligand metrics to rmsd_ligand_internal_{final,max}_a (preserved
# verbatim as ligand-on-ligand Kabsch RMSD diagnostics). The verdict is
# now keyed on pose RMSD. Free-energy soft gate E follows.
ARTIFACT_SCHEMA_VERSION = "1.1.0"

# Public verdict labels used by both engines and the frontend.
VERDICT_STABLE   = "stable"
VERDICT_DRIFTING = "drifting"
VERDICT_UNSTABLE = "unstable"
VERDICT_FAILED   = "failed"
