"""
Free-energy (MM/GBSA, MM/PBSA) refinement stage.

This stage is strictly downstream of MD and is gated by the MD verdict and
sampling adequacy. The gating function in `gating.py` is the *single source
of truth* — the GET /free-energy/gating endpoint, the POST /free-energy/run
guard, and the UI all call into it so frontend and backend never disagree
on whether a run is permitted.

The actual per-frame MM/GBSA energy evaluation is not bundled with this
repo (it requires AmberTools / OpenMM with explicit-solvent parameterization
of the ligand). When a run is permitted, this stage persists the gating
result, the resolved sampling plan, and the protocol/frames_used artifacts
so a later engine can pick them up; it then marks the calculation itself as
`status="planned"` rather than fabricating numbers.
"""

ARTIFACT_SCHEMA_VERSION = "1.0.0"

# --- Default thresholds (configurable per call) ---
DEFAULT_MIN_TOTAL_FRAMES = 200
DEFAULT_MIN_WINDOW_FRAMES = 50
DEFAULT_MAX_WINDOW_FRAMES_TARGET = 200
DEFAULT_LAST_FRACTION = 0.30
LAST_FRACTION_MIN = 0.10
LAST_FRACTION_MAX = 0.60

# Soft-gate (warning) thresholds.
LIGAND_MAX_RMSD_WARN_A   = 5.0
LIGAND_MAX_RMSD_BLOCK_A  = 10.0
BACKBONE_FINAL_RMSD_WARN_A = 2.5
WINDOW_DURATION_WARN_PS  = 1000.0   # < 1 ns post-eq → warn

# Reason / warning keys (machine-readable). Keep stable — the frontend
# tooltip-renders these and an external test harness asserts on them.
REASON_KEYS = {
    "MD_SUMMARY_MISSING":     "MD summary file not found.",
    "MD_NOT_COMPLETE":        "MD job did not complete.",
    "MD_NOT_STABLE":          "MD stability verdict is not stable.",
    # 2026-06-05 curation pass: a surrogate-engine trajectory is degenerate
    # (no real dynamics; the surrogate pre-aligns each conformer to the
    # docked pose and holds the receptor rigid). Free-energy estimates
    # derived from such a trajectory are not meaningful, so we hard-block
    # at the gating layer rather than waiting for a downstream "looks OK
    # but isn't" surprise. The check fires AFTER MD_NOT_STABLE so a real
    # OpenMM job that happened to fail stability still reports the right
    # blocker (verdict, not engine kind).
    "MD_NOT_OPENMM_FULL":     "MD trajectory engine is not OpenMM-full; surrogate runs are degenerate and not eligible for free-energy.",
    "MD_FRAMES_MISSING":      "MD frames directory not found.",
    "MD_FRAMES_EMPTY":        "MD frames directory is empty.",
    "MD_FRAME_UNREADABLE":    "First MD frame could not be opened/parsed.",
    "INSUFFICIENT_FRAMES":    "Not enough sampling frames after applying the window/stride.",
    "INVALID_PARAMETERS":     "Caller supplied invalid sampling parameters.",
}

WARNING_KEYS = {
    "LIGAND_DRIFT":           "Ligand RMSD trends suggest the pose may unbind during the window.",
    "BACKBONE_NOT_EQUILIBRATED": "Backbone RMSD has not plateaued; the chosen window may include drift.",
    "SHORT_WINDOW_TIME":      "Sampling window covers < 1 ns of simulation time; estimate is approximate.",
    "STRIDE_MANUAL":          "User-supplied stride; auto-stride was overridden.",
    "EXPERIMENTAL_ENGINE":    "No bundled MM/GBSA engine — the calculation step itself is planned.",
}
