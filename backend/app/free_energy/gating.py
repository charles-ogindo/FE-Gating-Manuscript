"""
Free-energy gating — single source of truth.

`validate(md_job_id, **knobs)` returns a `GatingResult` consumed by:

  - GET  /free-energy/gating   — UI calls this on every knob change
  - POST /free-energy/run      — refuses to start if can_run is False
  - tests/test_free_energy_gating.py — exercises the rules directly

Hard gates (block):
  A) MD job exists and summary.json says status == "completed".
  B) summary.verdict == "stable".
  B') summary.engine.kind == "openmm_full" — surrogate-engine trajectories
      are degenerate (the surrogate pre-aligns each conformer to the
      docked pose + holds the receptor rigid, so pose RMSD ≡ internal
      RMSD and there is no real dynamics to integrate); a surrogate
      verdict==stable is not eligible for free-energy. Added 2026-06-05.
  C) md/frames/ directory exists, is non-empty, and the first frame opens.
  D) Total frames >= MIN_TOTAL_FRAMES and the resolved sampling window has
     at least MIN_WINDOW_FRAMES *sampled* frames.

Soft gates (warn but allow):
  E) Ligand stays bound: warn if `metrics.rmsd_ligand_pose_max_a` exceeds
     LIGAND_MAX_RMSD_WARN_A; only *block* if it exceeds the much-stricter
     LIGAND_MAX_RMSD_BLOCK_A (10 Å — well past unbinding). Reads the Q6b
     POSE metric (receptor-frame ligand displacement) — the metric that
     actually answers "did the ligand stay in the pocket?". Pre-Q6b
     summary.json files emit `rmsd_ligand_max_a` (ligand-on-ligand); this
     gate falls back to that key when the pose key is absent, so pre-Q6b
     re-analysis is not strictly required for the gate to function — but
     the legacy metric is blind to pocket displacement, so a re-analysis
     pass is strongly recommended for accurate gating.
  F) Backbone equilibration: warn if `metrics.rmsd_backbone_final_a`
     exceeds BACKBONE_FINAL_RMSD_WARN_A.
  G) Outliers / unreadable frames: scan a small sample of frame PDBs, warn
     about how many couldn't be parsed.
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from backend.app.core.config import JOBS_DIR
from . import (
    BACKBONE_FINAL_RMSD_WARN_A,
    DEFAULT_LAST_FRACTION,
    DEFAULT_MIN_TOTAL_FRAMES,
    DEFAULT_MIN_WINDOW_FRAMES,
    LIGAND_MAX_RMSD_BLOCK_A,
    LIGAND_MAX_RMSD_WARN_A,
    REASON_KEYS,
    WARNING_KEYS,
    WINDOW_DURATION_WARN_PS,
)
from .sampling import SamplingPlan, compute_sampling_plan

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class GatingReason:
    key: str           # stable machine key (see REASON_KEYS / WARNING_KEYS)
    message: str       # human-readable, includes the offending values
    severity: str      # "blocker" | "warning"
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class GatingResult:
    can_run: bool
    reasons: List[GatingReason]
    warnings: List[GatingReason]
    resolved_sampling_plan: Optional[SamplingPlan]
    md_job_id: str
    md_status: Optional[str] = None
    md_verdict: Optional[str] = None
    thresholds: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "can_run": self.can_run,
            "reasons": [r.to_dict() for r in self.reasons],
            "warnings": [w.to_dict() for w in self.warnings],
            "resolved_sampling_plan":
                self.resolved_sampling_plan.to_dict()
                if self.resolved_sampling_plan else None,
            "md_job_id": self.md_job_id,
            "md_status": self.md_status,
            "md_verdict": self.md_verdict,
            "thresholds": self.thresholds,
        }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def validate(
    md_job_id: str,
    *,
    last_fraction: Optional[float] = None,
    stride: Optional[int] = None,
    min_total_frames: int = DEFAULT_MIN_TOTAL_FRAMES,
    min_window_frames: int = DEFAULT_MIN_WINDOW_FRAMES,
    jobs_dir: Optional[Path] = None,
) -> GatingResult:
    """
    Gate a free-energy run for the given MD job. Returns a structured result
    regardless of pass/fail — callers do not need to handle exceptions for
    the normal failure modes (the result describes them in `reasons`).
    """
    base = jobs_dir if jobs_dir is not None else JOBS_DIR
    thresholds = {
        "min_total_frames": int(min_total_frames),
        "min_window_frames": int(min_window_frames),
        "ligand_max_rmsd_warn_a": LIGAND_MAX_RMSD_WARN_A,
        "ligand_max_rmsd_block_a": LIGAND_MAX_RMSD_BLOCK_A,
        "backbone_final_rmsd_warn_a": BACKBONE_FINAL_RMSD_WARN_A,
    }

    reasons: List[GatingReason] = []
    warnings: List[GatingReason] = []

    # ------------------------------------------------------------------
    # Hard gate (A): summary exists and status == "completed"
    # ------------------------------------------------------------------
    md_root = Path(base) / md_job_id / "md"
    summary_path = md_root / "summary.json"
    if not summary_path.is_file():
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_SUMMARY_MISSING",
            f"No MD summary at {summary_path}",
            details={"summary_path": str(summary_path)},
        )

    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as e:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_SUMMARY_MISSING",
            f"MD summary could not be parsed: {e}",
            details={"summary_path": str(summary_path)},
        )

    md_status = summary.get("status")
    md_verdict = summary.get("verdict")

    if md_status != "completed":
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_NOT_COMPLETE",
            f"MD job status is {md_status!r}; need 'completed'.",
            details={"status": md_status},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Hard gate (B): verdict == "stable"
    # ------------------------------------------------------------------
    if md_verdict != "stable":
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_NOT_STABLE",
            f"MD stability verdict is {md_verdict!r}; need 'stable'. "
            f"Inspect MD results before requesting a free-energy run.",
            details={"verdict": md_verdict},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Hard gate (B'): engine.kind == "openmm_full" (2026-06-05).
    # Blocks surrogate trajectories from passing free-energy gating even
    # when their verdict happens to be "stable" — the surrogate's
    # pre-aligned, rigid-receptor construction makes pose RMSD ≡ internal
    # RMSD by construction, and any free-energy estimate from such a
    # trajectory is meaningless. Pre-Q9 summaries that pre-date the
    # engine.kind field (none in this codebase but possible from older
    # imports) treat a missing engine.kind as a block: we cannot tell
    # whether they were openmm_full or not.
    # ------------------------------------------------------------------
    engine_block = summary.get("engine") or {}
    md_engine_kind = engine_block.get("kind")
    if md_engine_kind != "openmm_full":
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_NOT_OPENMM_FULL",
            (
                f"MD engine is {md_engine_kind!r}; need 'openmm_full'. "
                "Surrogate trajectories are degenerate (rigid receptor, "
                "pre-aligned conformers) and not eligible for free-energy "
                "calculation."
            ),
            details={"engine_kind": md_engine_kind},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Hard gate (C): frames directory exists + first frame opens
    # ------------------------------------------------------------------
    frames_dir = md_root / "frames"
    if not frames_dir.is_dir():
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_FRAMES_MISSING",
            f"Frames directory missing at {frames_dir}",
            details={"frames_dir": str(frames_dir)},
            md_status=md_status, md_verdict=md_verdict,
        )

    frame_paths = sorted(frames_dir.glob("frame_*.pdb"))
    n_total = len(frame_paths)
    if n_total == 0:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_FRAMES_EMPTY",
            f"Frames directory at {frames_dir} contains no frame_*.pdb files.",
            details={"frames_dir": str(frames_dir)},
            md_status=md_status, md_verdict=md_verdict,
        )

    # First-frame openable check (also catches truncated PDB).
    try:
        head = frame_paths[0].read_text(encoding="utf-8", errors="ignore")
        if not any(line.startswith(("ATOM", "HETATM")) for line in head.splitlines()):
            raise RuntimeError("no ATOM/HETATM records")
    except Exception as e:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "MD_FRAME_UNREADABLE",
            f"First MD frame {frame_paths[0].name} could not be read: {e}",
            details={"frame_path": str(frame_paths[0])},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Hard gate (D-1): minimum *total* frames before any windowing
    # ------------------------------------------------------------------
    if n_total < min_total_frames:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "INSUFFICIENT_FRAMES",
            f"MD produced {n_total} frames; need at least {min_total_frames} total. "
            "Re-run MD with longer production (or smaller snapshot interval).",
            details={"n_total": n_total, "min_total_frames": min_total_frames},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Resolve the sampling plan once; downstream checks use it.
    # ------------------------------------------------------------------
    dt_ps = (summary.get("settings") or {}).get("snapshot_every_ps")
    try:
        plan = compute_sampling_plan(
            n_total_frames=n_total,
            last_fraction=last_fraction,
            stride=stride,
            snapshot_every_ps=dt_ps,
        )
    except ValueError as e:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "INVALID_PARAMETERS",
            f"Invalid sampling parameters: {e}",
            details={"last_fraction": last_fraction, "stride": stride},
            md_status=md_status, md_verdict=md_verdict,
        )

    # ------------------------------------------------------------------
    # Hard gate (D-2): post-stride sampled count >= min_window_frames
    # ------------------------------------------------------------------
    if plan.n_frames_sampled < min_window_frames:
        return _fail(
            md_job_id, reasons, warnings, thresholds,
            "INSUFFICIENT_FRAMES",
            (
                f"Resolved sampling plan yields {plan.n_frames_sampled} sampled "
                f"frame(s); need at least {min_window_frames}. "
                f"Window: last {plan.window_last_fraction*100:.0f}% "
                f"({plan.n_frames_in_window} raw frames), stride={plan.stride}. "
                "Try a larger window or a smaller manual stride."
            ),
            details={
                "n_total": n_total,
                "n_in_window": plan.n_frames_in_window,
                "stride": plan.stride,
                "n_sampled": plan.n_frames_sampled,
                "min_window_frames": min_window_frames,
            },
            md_status=md_status, md_verdict=md_verdict,
            resolved_sampling_plan=plan,
        )

    # ------------------------------------------------------------------
    # Soft gate (E): ligand drift (warn or strict-block at huge values).
    # Q6b: prefers the pose metric (receptor-frame displacement). Falls
    # back to the pre-Q6b ligand-on-ligand key so pre-Q6b summaries that
    # haven't been re-analyzed don't error — see module docstring.
    # ------------------------------------------------------------------
    metrics = summary.get("metrics") or {}
    lig_max = metrics.get("rmsd_ligand_pose_max_a")
    metric_key = "rmsd_ligand_pose_max_a"
    if lig_max is None:
        lig_max = metrics.get("rmsd_ligand_max_a")  # pre-Q6b fallback
        metric_key = "rmsd_ligand_max_a"
    if isinstance(lig_max, (int, float)) and math.isfinite(lig_max):
        if lig_max >= LIGAND_MAX_RMSD_BLOCK_A:
            return _fail(
                md_job_id, reasons, warnings, thresholds,
                "MD_NOT_STABLE",
                (
                    f"Ligand max RMSD in trajectory is {lig_max:.2f} Å "
                    f"(≥ {LIGAND_MAX_RMSD_BLOCK_A:.1f} Å) — ligand has likely "
                    "unbound. Free-energy estimate would be meaningless."
                ),
                details={metric_key: lig_max},
                md_status=md_status, md_verdict=md_verdict,
                resolved_sampling_plan=plan,
            )
        if lig_max >= LIGAND_MAX_RMSD_WARN_A:
            warnings.append(GatingReason(
                key="LIGAND_DRIFT",
                message=(
                    f"Ligand max RMSD {lig_max:.2f} Å in trajectory — pose may "
                    "be sliding; treat the free-energy estimate as approximate."
                ),
                severity="warning",
                details={metric_key: lig_max},
            ))

    # ------------------------------------------------------------------
    # Soft gate (F): backbone equilibration heuristic
    # ------------------------------------------------------------------
    bb_final = metrics.get("rmsd_backbone_final_a")
    if (isinstance(bb_final, (int, float))
            and math.isfinite(bb_final)
            and bb_final >= BACKBONE_FINAL_RMSD_WARN_A):
        warnings.append(GatingReason(
            key="BACKBONE_NOT_EQUILIBRATED",
            message=(
                f"Backbone Cα RMSD at trajectory end is {bb_final:.2f} Å "
                f"(≥ {BACKBONE_FINAL_RMSD_WARN_A:.1f} Å) — the chosen window "
                "may still include drift; consider tightening last_fraction."
            ),
            severity="warning",
            details={"rmsd_backbone_final_a": bb_final},
        ))

    # ------------------------------------------------------------------
    # Soft gate (G): outliers in the sampled frames
    # ------------------------------------------------------------------
    n_bad = _count_unreadable_frames(frame_paths, plan.frame_indices)
    if n_bad > 0:
        warnings.append(GatingReason(
            key="OUTLIER_FRAMES",
            message=(
                f"{n_bad} of {plan.n_frames_sampled} sampled frame(s) could not "
                "be parsed and will be filtered if MM/GBSA is run."
            ),
            severity="warning",
            details={"n_unreadable": n_bad, "n_sampled": plan.n_frames_sampled},
        ))

    # ------------------------------------------------------------------
    # Time-based warning: short post-eq window
    # ------------------------------------------------------------------
    if dt_ps is not None and isinstance(dt_ps, (int, float)) and dt_ps > 0:
        window_duration_ps = plan.n_frames_in_window * float(dt_ps)
        if window_duration_ps < WINDOW_DURATION_WARN_PS:
            warnings.append(GatingReason(
                key="SHORT_WINDOW_TIME",
                message=(
                    f"Sampling window covers {window_duration_ps:.0f} ps "
                    f"(< {WINDOW_DURATION_WARN_PS:.0f} ps); free-energy "
                    "estimate is approximate."
                ),
                severity="warning",
                details={"window_duration_ps": window_duration_ps},
            ))

    # Manual-stride note (informational warning so the protocol records it).
    if not plan.stride_was_auto:
        warnings.append(GatingReason(
            key="STRIDE_MANUAL",
            message=(
                f"Stride {plan.stride} supplied by the caller; auto-stride was "
                "overridden."
            ),
            severity="warning",
            details={"stride": plan.stride},
        ))

    # Standing flag that no engine is bundled — frontend renders this as the
    # "Experimental / approximate" chip so the user can't miss it.
    warnings.append(GatingReason(
        key="EXPERIMENTAL_ENGINE",
        message=WARNING_KEYS["EXPERIMENTAL_ENGINE"],
        severity="warning",
    ))

    return GatingResult(
        can_run=True,
        reasons=reasons,
        warnings=warnings,
        resolved_sampling_plan=plan,
        md_job_id=md_job_id,
        md_status=md_status,
        md_verdict=md_verdict,
        thresholds=thresholds,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _fail(
    md_job_id: str,
    reasons: List[GatingReason],
    warnings: List[GatingReason],
    thresholds: Dict[str, Any],
    key: str,
    message: str,
    *,
    details: Optional[Dict[str, Any]] = None,
    md_status: Optional[str] = None,
    md_verdict: Optional[str] = None,
    resolved_sampling_plan: Optional[SamplingPlan] = None,
) -> GatingResult:
    if key not in REASON_KEYS:
        # Defensive: keep machine keys curated. Don't crash, just log.
        logger.warning("Unknown gating reason key: %s", key)
    reasons.append(GatingReason(
        key=key,
        message=message,
        severity="blocker",
        details=details or {},
    ))
    return GatingResult(
        can_run=False,
        reasons=reasons,
        warnings=warnings,
        resolved_sampling_plan=resolved_sampling_plan,
        md_job_id=md_job_id,
        md_status=md_status,
        md_verdict=md_verdict,
        thresholds=thresholds,
    )


def _count_unreadable_frames(
    frame_paths: List[Path],
    indices: List[int],
    sample_cap: int = 20,
) -> int:
    """Spot-check up to `sample_cap` frames for readability (cheap heuristic)."""
    if not indices:
        return 0
    step = max(1, len(indices) // sample_cap)
    bad = 0
    checked = 0
    for i in indices[::step]:
        if i >= len(frame_paths):
            continue
        try:
            text = frame_paths[i].read_text(encoding="utf-8", errors="ignore")
            if not any(line.startswith(("ATOM", "HETATM")) for line in text.splitlines()):
                bad += 1
        except Exception:
            bad += 1
        checked += 1
        if checked >= sample_cap:
            break
    return bad
