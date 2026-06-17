"""
Unit tests for free_energy.sampling.compute_sampling_plan and
free_energy.gating.validate. Synthetic MD fixtures are created in tmp_path;
no real MD artifacts are touched.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from backend.app.free_energy.sampling import compute_sampling_plan
from backend.app.free_energy.gating import validate


# ---------------------------------------------------------------------------
# Fixture: a minimal but realistic MD job tree on disk
# ---------------------------------------------------------------------------
def _write_md_fixture(
    jobs_dir: Path,
    md_job_id: str,
    *,
    n_frames: int,
    status: str = "completed",
    verdict: str = "stable",
    snapshot_every_ps: float = 5.0,
    rmsd_ligand_pose_max_a: float | None = None,
    rmsd_backbone_final_a: float | None = None,
    engine_kind: str = "openmm_full",
    bad_frame_idx: tuple[int, ...] = (),
) -> Path:
    """Create jobs/<md_job_id>/md/summary.json + frames/frame_NNN.pdb.

    Q6b: parameter is `rmsd_ligand_pose_max_a` (pose, the primary metric
    soft gate E reads). The pre-Q6b legacy key `rmsd_ligand_max_a` has a
    dedicated regression test below that still exercises gating.py's
    fallback path for old summaries in the wild.
    """
    md_dir = jobs_dir / md_job_id / "md"
    md_dir.mkdir(parents=True, exist_ok=True)
    frames = md_dir / "frames"
    frames.mkdir(exist_ok=True)

    summary = {
        "schema_version": "1.1.0",
        "md_job_id": md_job_id,
        "status": status,
        "verdict": verdict,
        "n_frames": n_frames,
        "settings": {"snapshot_every_ps": snapshot_every_ps},
        "engine": {"kind": engine_kind},
        "metrics": {
            "rmsd_ligand_pose_max_a": rmsd_ligand_pose_max_a,
            "rmsd_backbone_final_a": rmsd_backbone_final_a,
        },
    }
    (md_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

    # One ATOM record per frame is enough for gating's readability check.
    for i in range(n_frames):
        content = (
            "" if i in bad_frame_idx
            else "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n"
        )
        (frames / f"frame_{i:03d}.pdb").write_text(content, encoding="utf-8")

    return md_dir


# ===========================================================================
# Sampling plan
# ===========================================================================
class TestSamplingPlan:
    def test_default_30pct_window_auto_stride(self):
        plan = compute_sampling_plan(n_total_frames=300)
        # last 30% = frames 210..299 (90 frames). 90 <= 200 → stride 1.
        assert plan.window_start_frame == 210
        assert plan.window_end_frame == 300
        assert plan.window_last_fraction == pytest.approx(0.30)
        assert plan.stride == 1
        assert plan.stride_was_auto is True
        assert plan.n_frames_in_window == 90
        assert plan.n_frames_sampled == 90
        assert plan.frame_indices[0] == 210
        assert plan.frame_indices[-1] == 299

    def test_large_window_auto_stride_caps_sampled_count(self):
        # 5000 frames, last 30% = 1500. Default target = 200.
        plan = compute_sampling_plan(n_total_frames=5000)
        assert plan.stride == max(1, -(-1500 // 200))  # ceil(1500/200) = 8
        assert plan.n_frames_sampled <= 200
        assert plan.n_frames_sampled >= 50

    def test_manual_stride_honored(self):
        plan = compute_sampling_plan(n_total_frames=300, stride=5)
        assert plan.stride == 5
        assert plan.stride_was_auto is False
        # window has 90 frames, stride 5 → ceil(90/5) = 18
        assert plan.n_frames_sampled == 18

    def test_last_fraction_clamped_low(self):
        plan = compute_sampling_plan(n_total_frames=1000, last_fraction=0.01)
        assert plan.window_last_fraction == pytest.approx(0.10)

    def test_last_fraction_clamped_high(self):
        plan = compute_sampling_plan(n_total_frames=1000, last_fraction=0.99)
        assert plan.window_last_fraction == pytest.approx(0.60)

    def test_dt_propagates_times(self):
        plan = compute_sampling_plan(
            n_total_frames=300, last_fraction=0.30,
            snapshot_every_ps=5.0,
        )
        assert plan.times_ps is not None
        assert plan.times_ps[0] == pytest.approx(210 * 5.0)
        assert plan.times_ps[-1] == pytest.approx(299 * 5.0)

    def test_empty_trajectory(self):
        plan = compute_sampling_plan(n_total_frames=0)
        assert plan.n_frames_in_window == 0
        assert plan.n_frames_sampled == 0
        assert plan.frame_indices == []

    def test_negative_total_rejected(self):
        with pytest.raises(ValueError):
            compute_sampling_plan(n_total_frames=-1)


# ===========================================================================
# Gating — hard gates
# ===========================================================================
class TestGatingHardGates:
    def test_missing_summary(self, tmp_path):
        result = validate("nonexistent", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_SUMMARY_MISSING"

    def test_md_not_complete(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=300, status="running")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_NOT_COMPLETE"
        assert result.md_status == "running"

    def test_md_not_stable(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=300, verdict="drifting")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_NOT_STABLE"
        assert result.md_verdict == "drifting"

    def test_frames_missing(self, tmp_path):
        # Build summary but skip writing frames dir. Must include the
        # engine block since 2026-06-05 added a B' hard gate on
        # engine.kind == "openmm_full" that fires before frame checks.
        md = tmp_path / "md1" / "md"
        md.mkdir(parents=True)
        (md / "summary.json").write_text(json.dumps({
            "status": "completed",
            "verdict": "stable",
            "engine": {"kind": "openmm_full"},
        }), encoding="utf-8")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_FRAMES_MISSING"

    def test_frames_empty(self, tmp_path):
        md = tmp_path / "md1" / "md"
        (md / "frames").mkdir(parents=True)
        (md / "summary.json").write_text(json.dumps({
            "status": "completed",
            "verdict": "stable",
            "engine": {"kind": "openmm_full"},
        }), encoding="utf-8")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_FRAMES_EMPTY"

    def test_first_frame_unreadable(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=300, bad_frame_idx=(0,))
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_FRAME_UNREADABLE"

    def test_too_few_total_frames(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=100)
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "INSUFFICIENT_FRAMES"
        assert result.reasons[0].details["n_total"] == 100

    def test_too_few_window_frames_after_stride(self, tmp_path):
        # 250 frames, last 30% = 75 raw frames; stride=10 → 8 sampled < 50.
        _write_md_fixture(tmp_path, "md1", n_frames=250)
        result = validate("md1", jobs_dir=tmp_path, stride=10)
        assert not result.can_run
        assert result.reasons[0].key == "INSUFFICIENT_FRAMES"

    def test_huge_ligand_drift_blocks(self, tmp_path):
        _write_md_fixture(
            tmp_path, "md1", n_frames=400,
            rmsd_ligand_pose_max_a=15.0,
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        # We reuse MD_NOT_STABLE for "ligand has clearly unbound".
        assert result.reasons[0].key == "MD_NOT_STABLE"

    def test_surrogate_engine_blocked_even_when_stable(self, tmp_path):
        # 2026-06-05 curation pass: even a verdict=='stable' surrogate
        # trajectory must be rejected by free-energy gating. The
        # surrogate's pre-aligned + rigid-receptor construction makes
        # pose RMSD ≡ internal RMSD by construction, so any free-energy
        # estimate from such a trajectory is meaningless.
        _write_md_fixture(
            tmp_path, "md1", n_frames=400,
            verdict="stable",
            rmsd_ligand_pose_max_a=1.5, rmsd_backbone_final_a=1.5,
            engine_kind="surrogate",
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_NOT_OPENMM_FULL"
        assert result.reasons[0].details["engine_kind"] == "surrogate"

    def test_missing_engine_kind_is_blocked(self, tmp_path):
        # Defensive: a pre-engine-field summary (engine block missing or
        # engine.kind == None) cannot be confirmed as openmm_full, so the
        # gate blocks. Operator must re-run on openmm_full to enable FE.
        md_dir = tmp_path / "md1" / "md"
        md_dir.mkdir(parents=True)
        (md_dir / "frames").mkdir()
        for i in range(400):
            (md_dir / "frames" / f"frame_{i:03d}.pdb").write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n"
            )
        summary = {
            "schema_version": "1.1.0",
            "md_job_id": "md1",
            "status": "completed",
            "verdict": "stable",
            "n_frames": 400,
            "settings": {"snapshot_every_ps": 5.0},
            # NO `engine` block at all
            "metrics": {
                "rmsd_ligand_pose_max_a": 1.5,
                "rmsd_backbone_final_a": 1.5,
            },
        }
        (md_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_NOT_OPENMM_FULL"
        assert result.reasons[0].details["engine_kind"] is None

    def test_pre_q6b_legacy_key_still_read(self, tmp_path):
        # Defensive regression: a pre-Q6b summary.json (schema 1.0.0) emits
        # rmsd_ligand_max_a (ligand-on-ligand), not rmsd_ligand_pose_max_a.
        # The gate's fallback path must keep reading the legacy key so
        # pre-Q6b jobs that haven't been re-analyzed don't crash the gate.
        # Write the summary by hand because _write_md_fixture now emits
        # the new key only.
        md_dir = tmp_path / "md1" / "md"
        md_dir.mkdir(parents=True)
        (md_dir / "frames").mkdir()
        for i in range(400):
            (md_dir / "frames" / f"frame_{i:03d}.pdb").write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n"
            )
        summary = {
            "schema_version": "1.0.0",
            "md_job_id": "md1",
            "status": "completed",
            "verdict": "stable",
            "n_frames": 400,
            "settings": {"snapshot_every_ps": 5.0},
            "engine": {"kind": "openmm_full"},
            "metrics": {
                "rmsd_ligand_max_a": 15.0,   # pre-Q6b key
                "rmsd_backbone_final_a": 1.5,
            },
        }
        (md_dir / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        result = validate("md1", jobs_dir=tmp_path)
        assert not result.can_run
        assert result.reasons[0].key == "MD_NOT_STABLE"
        # And the details payload uses the legacy key name, not pose, so
        # downstream consumers can tell which path fired.
        assert "rmsd_ligand_max_a" in result.reasons[0].details


# ===========================================================================
# Gating — happy path + warnings
# ===========================================================================
class TestGatingPasses:
    def test_basic_pass(self, tmp_path):
        _write_md_fixture(
            tmp_path, "md1", n_frames=400,
            rmsd_ligand_pose_max_a=2.0, rmsd_backbone_final_a=1.5,
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert result.can_run is True
        assert result.reasons == []
        assert result.md_status == "completed"
        assert result.md_verdict == "stable"
        plan = result.resolved_sampling_plan
        assert plan is not None
        assert plan.n_frames_sampled >= 50

    def test_warning_short_window_time(self, tmp_path):
        # 400 frames * 0.5 ps = 200 ps total, window 30% = 60 ps << 1 ns.
        _write_md_fixture(
            tmp_path, "md1", n_frames=400,
            snapshot_every_ps=0.5,
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert result.can_run is True
        keys = {w.key for w in result.warnings}
        assert "SHORT_WINDOW_TIME" in keys

    def test_warning_ligand_drift(self, tmp_path):
        _write_md_fixture(
            tmp_path, "md1", n_frames=400, rmsd_ligand_pose_max_a=7.0,
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert result.can_run is True
        keys = {w.key for w in result.warnings}
        assert "LIGAND_DRIFT" in keys

    def test_warning_backbone_not_equilibrated(self, tmp_path):
        _write_md_fixture(
            tmp_path, "md1", n_frames=400, rmsd_backbone_final_a=3.5,
        )
        result = validate("md1", jobs_dir=tmp_path)
        assert result.can_run is True
        keys = {w.key for w in result.warnings}
        assert "BACKBONE_NOT_EQUILIBRATED" in keys

    def test_manual_stride_warning(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=1000)
        result = validate("md1", jobs_dir=tmp_path, stride=2)
        assert result.can_run is True
        keys = {w.key for w in result.warnings}
        assert "STRIDE_MANUAL" in keys

    def test_experimental_engine_warning_always_present(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=400)
        result = validate("md1", jobs_dir=tmp_path)
        keys = {w.key for w in result.warnings}
        assert "EXPERIMENTAL_ENGINE" in keys


# ===========================================================================
# Window override behavior
# ===========================================================================
class TestWindowOverrides:
    def test_user_can_widen_window_to_pass(self, tmp_path):
        # 250 frames, default 30% window = 75 raw, stride auto = 1 → 75 sampled.
        # That should pass. But with a manual stride of 2 → 38 sampled → fail.
        _write_md_fixture(tmp_path, "md1", n_frames=250)
        fail = validate("md1", jobs_dir=tmp_path, stride=2)
        assert not fail.can_run
        # Widening last_fraction to 0.60 → 150 raw, stride 2 → 75 sampled → pass.
        ok = validate("md1", jobs_dir=tmp_path, last_fraction=0.60, stride=2)
        assert ok.can_run is True

    def test_min_window_threshold_override(self, tmp_path):
        _write_md_fixture(tmp_path, "md1", n_frames=250)
        # With min_window_frames=20, the stride=2 case (38 sampled) should pass.
        ok = validate(
            "md1", jobs_dir=tmp_path,
            stride=2, min_window_frames=20,
        )
        assert ok.can_run is True
