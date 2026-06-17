"""
Sampling-plan computation.

Pure logic: given (n_total_frames, last_fraction, stride, dt) → SamplingPlan.
No filesystem, no MD summary parsing — that lives in gating.py. Splitting
keeps the stride/window arithmetic unit-testable in isolation.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import List, Optional

from . import (
    DEFAULT_LAST_FRACTION,
    DEFAULT_MAX_WINDOW_FRAMES_TARGET,
    LAST_FRACTION_MAX,
    LAST_FRACTION_MIN,
)


@dataclass
class SamplingPlan:
    """Resolved plan for which MD frames will be fed to the free-energy step."""
    window_start_frame: int          # inclusive
    window_end_frame: int            # exclusive
    window_last_fraction: float      # what we actually used after clamping
    stride: int
    stride_was_auto: bool
    n_frames_in_window: int          # window size before stride
    n_frames_sampled: int            # after stride
    frame_indices: List[int]
    times_ps: Optional[List[float]] = field(default=None)
    snapshot_every_ps: Optional[float] = field(default=None)

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


def compute_sampling_plan(
    *,
    n_total_frames: int,
    last_fraction: Optional[float] = None,
    stride: Optional[int] = None,
    snapshot_every_ps: Optional[float] = None,
    max_window_frames_target: int = DEFAULT_MAX_WINDOW_FRAMES_TARGET,
) -> SamplingPlan:
    """
    Build a SamplingPlan from frame totals + user knobs.

    Rules:
      - last_fraction is clamped to [LAST_FRACTION_MIN, LAST_FRACTION_MAX].
      - If stride is None / <= 0 it is chosen automatically so the post-stride
        count does not exceed `max_window_frames_target`.
      - n_frames_sampled is the count *after* stride application.

    This function never raises on "too few frames" — that policy decision
    belongs to gating.validate(). It will, however, raise ValueError on
    unambiguously broken inputs (negative totals, etc.).
    """
    if n_total_frames < 0:
        raise ValueError(f"n_total_frames must be >= 0, got {n_total_frames}")
    if stride is not None and stride < 0:
        raise ValueError(f"stride must be >= 0 or None, got {stride}")

    lf = last_fraction if last_fraction is not None else DEFAULT_LAST_FRACTION
    if not math.isfinite(lf):
        raise ValueError(f"last_fraction must be finite, got {lf}")
    lf = max(LAST_FRACTION_MIN, min(LAST_FRACTION_MAX, lf))

    if n_total_frames == 0:
        return SamplingPlan(
            window_start_frame=0,
            window_end_frame=0,
            window_last_fraction=lf,
            stride=1,
            stride_was_auto=stride in (None, 0),
            n_frames_in_window=0,
            n_frames_sampled=0,
            frame_indices=[],
            times_ps=None,
            snapshot_every_ps=snapshot_every_ps,
        )

    # Take the last `lf` fraction of frames, but always include at least one.
    window_size = max(1, int(round(n_total_frames * lf)))
    window_start = max(0, n_total_frames - window_size)
    window_end = n_total_frames
    n_window = window_end - window_start

    # Auto-stride: keep n_sampled <= max_window_frames_target.
    auto = stride is None or stride == 0
    if auto:
        eff_stride = max(1, math.ceil(n_window / max_window_frames_target))
    else:
        eff_stride = max(1, int(stride))

    frame_indices = list(range(window_start, window_end, eff_stride))
    n_sampled = len(frame_indices)

    times_ps: Optional[List[float]] = None
    if snapshot_every_ps is not None and snapshot_every_ps > 0:
        times_ps = [i * float(snapshot_every_ps) for i in frame_indices]

    return SamplingPlan(
        window_start_frame=window_start,
        window_end_frame=window_end,
        window_last_fraction=lf,
        stride=eff_stride,
        stride_was_auto=auto,
        n_frames_in_window=n_window,
        n_frames_sampled=n_sampled,
        frame_indices=frame_indices,
        times_ps=times_ps,
        snapshot_every_ps=snapshot_every_ps,
    )
