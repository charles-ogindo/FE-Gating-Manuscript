"""
Free-energy protocol + frames_used artifact writers.

When POST /free-energy/run is called and gating passes, the route writes:

  jobs/<md_job_id>/free_energy/protocol.json    — method + dielectrics +
                                                  ionic strength + sampling
                                                  knobs + thresholds.
  jobs/<md_job_id>/free_energy/frames_used.json — exact frame indices,
                                                  relative paths, and times.
  jobs/<md_job_id>/free_energy/summary.json     — top-level run state. Marked
                                                  status="planned" because no
                                                  MM/GBSA engine is bundled;
                                                  a downstream engine can flip
                                                  it to "completed" later.
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from backend.app.core.config import JOBS_DIR
from . import ARTIFACT_SCHEMA_VERSION
from .gating import GatingResult
from .sampling import SamplingPlan


@dataclass
class FreeEnergyProtocol:
    method: str = "MM-GBSA"                      # "MM-GBSA" | "MM-PBSA"
    implicit_solvent: str = "OBC2"
    solute_dielectric: float = 1.0
    solvent_dielectric: float = 80.0
    ionic_strength_mM: float = 150.0
    salt_radii_model: str = "mbondi2"
    surface_tension_kcal_mol_A2: float = 0.0072  # standard MM-GBSA SA term

    def to_dict(self) -> dict:
        return asdict(self)


def fe_dir(md_job_id: str, jobs_dir: Optional[Path] = None) -> Path:
    base = jobs_dir if jobs_dir is not None else JOBS_DIR
    return Path(base) / md_job_id / "free_energy"


def write_run_artifacts(
    md_job_id: str,
    gating: GatingResult,
    protocol: FreeEnergyProtocol,
    *,
    jobs_dir: Optional[Path] = None,
) -> Dict[str, str]:
    """
    Persist the three canonical artifacts. Returns the relative paths under
    jobs/<md_job_id>/ so the API response can hand them straight to the UI.
    """
    if not gating.can_run:
        raise RuntimeError(
            "write_run_artifacts called on a failed gating result; "
            "caller must check can_run first."
        )
    plan = gating.resolved_sampling_plan
    if plan is None:
        raise RuntimeError("gating result is missing the sampling plan")

    out = fe_dir(md_job_id, jobs_dir=jobs_dir)
    out.mkdir(parents=True, exist_ok=True)

    now = _dt.datetime.now().isoformat(timespec="seconds")

    protocol_doc = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "md_job_id": md_job_id,
        "created_at": now,
        "method": protocol.method,
        "implicit_solvent": protocol.implicit_solvent,
        "solute_dielectric": protocol.solute_dielectric,
        "solvent_dielectric": protocol.solvent_dielectric,
        "ionic_strength_mM": protocol.ionic_strength_mM,
        "salt_radii_model": protocol.salt_radii_model,
        "surface_tension_kcal_mol_A2": protocol.surface_tension_kcal_mol_A2,
        "sampling": plan.to_dict(),
        "thresholds": gating.thresholds,
        "warnings_at_gating": [w.to_dict() for w in gating.warnings],
    }
    (out / "protocol.json").write_text(
        json.dumps(protocol_doc, indent=2), encoding="utf-8")

    frames_doc = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "md_job_id": md_job_id,
        "created_at": now,
        "n_frames": plan.n_frames_sampled,
        "snapshot_every_ps": plan.snapshot_every_ps,
        "frame_indices": plan.frame_indices,
        "times_ps": plan.times_ps,
        "frame_paths": [
            f"md/frames/frame_{i:03d}.pdb" for i in plan.frame_indices
        ],
        "window_start_frame": plan.window_start_frame,
        "window_end_frame": plan.window_end_frame,
        "stride": plan.stride,
        "stride_was_auto": plan.stride_was_auto,
    }
    (out / "frames_used.json").write_text(
        json.dumps(frames_doc, indent=2), encoding="utf-8")

    summary_doc = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "md_job_id": md_job_id,
        "created_at": now,
        "status": "planned",
        "verdict": "not_implemented",
        "method": protocol.method,
        "reason": (
            "Gating passed and the sampling plan + protocol have been "
            "persisted, but no MM/GBSA / MM/PBSA engine is bundled with this "
            "repo. A downstream engine should read protocol.json + "
            "frames_used.json and update this summary in place."
        ),
        "artifacts": {
            "protocol":      "free_energy/protocol.json",
            "frames_used":   "free_energy/frames_used.json",
            "summary":       "free_energy/summary.json",
        },
        "gating": {
            "can_run": gating.can_run,
            "warnings": [w.to_dict() for w in gating.warnings],
            "thresholds": gating.thresholds,
        },
    }
    (out / "summary.json").write_text(
        json.dumps(summary_doc, indent=2), encoding="utf-8")

    return {
        "protocol":      f"free_energy/protocol.json",
        "frames_used":   f"free_energy/frames_used.json",
        "summary":       f"free_energy/summary.json",
    }
