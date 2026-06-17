# FILE: backend/app/core/config.py
#
# REVIEWER-REPO TRIMMED SHIM.
# The full application's config.py (scenario rules, safety-pipeline weights,
# database URL, CORS, etc.) is product configuration and is intentionally NOT
# included in this evidence repository. The copied scientific modules
# (free_energy/*, md/*, docking/md_receptor_prep, utils/rmsd) only import
# JOBS_DIR from here, so this shim provides exactly that and nothing more.
#
# JOBS_DIR is where per-job artifacts (MD frames, summary.json, etc.) live.
# Override with the JOBS_DIR env var if you keep job outputs elsewhere.

import os
from pathlib import Path

# Project root (the directory that contains backend/, scripts/, examples/).
BASE_DIR = Path(__file__).resolve().parents[3]

# Runtime data directory for per-job artifacts.
JOBS_DIR = Path(os.environ.get("JOBS_DIR", str(BASE_DIR / "jobs")))
JOBS_DIR.mkdir(parents=True, exist_ok=True)
