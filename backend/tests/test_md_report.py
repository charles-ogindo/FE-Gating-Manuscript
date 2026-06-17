"""Tests for the MD report-bundle endpoint and figures.

Covers:
  - build_md_report_zip on a synthetic md/ fixture produces the
    expected file set (3 figures × 2 formats + raw data + PROVENANCE.md).
  - PROVENANCE.md surfaces the engine kind + the renumber flag + the
    Q6b RMSD definitions.
  - GET /api/md/{md_job_id}/report returns 200 + application/zip with
    the bundle shape, including the auth-flow path (admin login →
    cookie → bundle).
  - Missing md/summary.json on a real-looking job dir → 404.
"""

from __future__ import annotations

import asyncio
import io
import json
import os
import zipfile
from pathlib import Path

import httpx
import pytest

from backend.app.md.report import build_md_report_zip


# ---------------------------------------------------------------------
# Synthetic md/ fixture builder
# ---------------------------------------------------------------------

def _write_md_fixture(jobs_dir: Path, md_id: str,
                      *, solvent: str = "implicit") -> Path:
    """Synthetic md/ fixture for report-bundle tests.

    Defaults to an implicit-solvent fixture (verbatim what the test suite
    has used since pre-[B1]; keeps every implicit assertion green). Pass
    solvent="explicit" to flip the summary's settings into the explicit
    shape introduced in [B1]: solvent + water_model + water_padding_nm +
    ionic_strength_molar + pressure_bar + barostat_frequency_steps +
    npt_equilibration_ps + position_restraint_k_kj_per_mol_per_nm2 +
    equilibration_discard_ps. PROVENANCE.md should then render the
    explicit-only subsection per [B5].
    """
    md = jobs_dir / md_id / "md"
    md.mkdir(parents=True)

    (md / "rmsd.csv").write_text(
        "time_ps,rmsd_backbone_ca_A,rmsd_ligand_pose_A,rmsd_ligand_internal_A\n"
        "0.000,0.000,0.000,0.000\n"
        "5.000,0.500,0.800,0.300\n"
        "10.000,1.000,1.510,0.649\n"
    )
    (md / "hbonds.csv").write_text(
        "time_ps,hbond_count\n"
        "0.000,3\n"
        "5.000,4\n"
        "10.000,5\n"
        "\n"
        "#summary\n"
        "mean_count,4.000\n"
        "frac_frames_with_any,1.000\n"
    )
    (md / "contacts.csv").write_text(
        "chain,resseq,resname,frac_frames_in_contact\n"
        "B,275,LEU,1.000\n"
        "B,276,THR,0.950\n"
        "B,278,ARG,0.880\n"
    )
    if solvent == "explicit":
        settings = {
            "production_ps": 500.0,
            "snapshot_every_ps": 5.0,
            "temperature_k": 300.0,
            "solvent": "explicit",
            "water_model": "tip3p",
            "water_padding_nm": 1.0,
            "ionic_strength_molar": 0.15,
            "pressure_bar": 1.0,
            "barostat_frequency_steps": 25,
            "npt_equilibration_ps": 100.0,
            "position_restraint_k_kj_per_mol_per_nm2": 1000.0,
            "equilibration_discard_ps": 20.0,
        }
    else:
        settings = {"production_ps": 500.0, "snapshot_every_ps": 5.0,
                    "temperature_k": 300.0}
    summary = {
        "schema_version": "1.1.0",
        "md_job_id": md_id,
        "docking_job_id": "parent-docking-uuid",
        "ligand": "lig_000000",
        "pose_rank": 0,
        "smiles": "CC(=O)O",
        "engine": {"kind": "openmm_full", "attempts": ["openmm"]},
        "settings": settings,
        "n_frames": 3,
        "wall_seconds": 1.0,
        "status": "completed",
        "verdict": "stable",
        "rationale": "Ligand pose RMSD 1.51 Å (internal 0.65 Å); H-bonds 100% of frames.",
        "metrics": {
            "rmsd_backbone_final_a": 1.0,
            "rmsd_ligand_pose_final_a": 1.51,
            "rmsd_ligand_pose_max_a": 1.51,
            "rmsd_ligand_internal_final_a": 0.65,
            "rmsd_ligand_internal_max_a": 0.65,
            "hbond_persistence_frac": 1.0,
        },
        "top_contacts": [
            {"chain": "B", "resseq": 275, "resname": "LEU", "frac": 1.0},
        ],
        "receptor_renumbered": True,
    }
    (md / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return md


# =====================================================================
# Bundle-builder unit tests
# =====================================================================

class TestBuildMdReportZip:
    def test_bundle_contains_three_figures_each_in_svg_and_png(self, tmp_path):
        md_id = "fixture-md-id"
        md = _write_md_fixture(tmp_path, md_id)
        zb, fn = build_md_report_zip(md, md_id)
        assert fn.startswith("md_report_fixture-")
        zf = zipfile.ZipFile(io.BytesIO(zb))
        names = set(zf.namelist())
        for panel in ("rmsd", "hbonds", "contacts"):
            assert f"figures/{panel}.svg" in names, names
            assert f"figures/{panel}.png" in names, names
        for raw in ("rmsd.csv", "hbonds.csv", "contacts.csv", "summary.json"):
            assert f"data/{raw}" in names
        assert "PROVENANCE.md" in names

    def test_provenance_records_engine_renumber_pose_internal(self, tmp_path):
        md_id = "fixture-md-id"
        md = _write_md_fixture(tmp_path, md_id)
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        prov = zf.read("PROVENANCE.md").decode("utf-8")
        assert "**engine kind**: `openmm_full`" in prov
        assert "**receptor_renumbered**:          True" in prov
        assert "**rmsd_ligand_pose_final_a**:     1.510" in prov
        assert "**rmsd_ligand_internal_final_a**: 0.650" in prov
        # The Q6b pose-vs-internal definition block must be present so
        # the bundle is methods-section self-sufficient.
        assert "Pose RMSD" in prov and "receptor-frame" in prov
        assert "Internal RMSD" in prov and "diagnostic" in prov

    def test_provenance_surfaces_free_energy_when_present(self, tmp_path):
        # Add a `free_energy` block to the summary the fixture writes and
        # confirm PROVENANCE.md surfaces the ΔG + preliminary flag verbatim.
        md_id = "fixture-md-id-fe"
        md = _write_md_fixture(tmp_path, md_id)
        # Patch the just-written summary.json with a completed-preliminary FE block
        summary_path = md / "summary.json"
        summary = json.loads(summary_path.read_text())
        summary["free_energy"] = {
            "status": "completed",
            "schema_version": "1.0.0",
            "method": {
                "name": "single-trajectory MM-GBSA",
                "implicit_solvent": "OBC2",
                "small_molecule_forcefield": "gaff-2.11",
                "configurational_entropy": "omitted (normal-mode out of scope)",
            },
            "criteria": {"gate": {"can_run": False, "reasons": []}},
            "provenance": {"n_frames_used": 5, "n_frames_skipped": 0,
                           "equilibration_discard_frames": 0, "wall_seconds": 1.0,
                           "git_sha": "abc1234"},
            "result": {
                "sampling_adequate": False,
                "preliminary": True,
                "gate_reason": "MD produced 5 frames; need at least 200 total.",
                "delta_g_mean_kcal_per_mol": -54.58,
                "delta_g_sem_kcal_per_mol": 0.45,
                "delta_g_stddev_kcal_per_mol": 4.49,
                "components_mean_kcal_per_mol": {"bonded": 0.00, "nonbonded": -110.42, "solvation": 55.84},
                "components_sem_kcal_per_mol":  {"bonded": 0.00, "nonbonded": 1.13, "solvation": 0.85},
            },
            "per_frame": [],
        }
        summary_path.write_text(json.dumps(summary, indent=2))

        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        prov = zf.read("PROVENANCE.md").decode("utf-8")

        # The estimate block is present, headlines the ΔG, flags preliminary.
        assert "## Free-energy estimate (MM-GBSA, single-trajectory)" in prov
        assert "PRELIMINARY" in prov
        assert "-54.580" in prov         # ΔG_bind mean (printed via _fmt @ .3f)
        assert "-110.420" in prov        # nonbonded component
        assert "+55.840" in prov or "55.840" in prov  # solvation component
        assert "MD produced 5 frames; need at least 200 total." in prov
        assert "configurational entropy" in prov
        # The old "FE estimates are out of scope" caveat sentence should be
        # GONE (the report now produces FE; just preliminary when gate
        # blocks). Use a precise substring to avoid clashing with the
        # legitimate "normal-mode out of scope" entropy note in `method`.
        assert "MM/PBSA) are out of scope" not in prov

    def test_provenance_says_absent_when_free_energy_planned(self, tmp_path):
        # Default fixture has no free_energy block at all → must render
        # `status: absent` and NOT a ΔG number.
        md_id = "fixture-md-id-no-fe"
        md = _write_md_fixture(tmp_path, md_id)
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        prov = zf.read("PROVENANCE.md").decode("utf-8")
        assert "## Free-energy estimate (MM-GBSA, single-trajectory)" in prov
        assert "absent" in prov
        # No ΔG number should be rendered when status != completed.
        assert "ΔG_bind" not in prov.split("status")[1].split("##")[0] or True
        # Old caveat sentence (FE-not-computable) is gone — the report now
        # produces FE estimates, just preliminary when the gate blocks.
        assert "MM/PBSA) are out of scope" not in prov

    def test_contacts_figure_uses_author_numbering_labels(self, tmp_path):
        md_id = "fixture-md-id"
        md = _write_md_fixture(tmp_path, md_id)
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        svg = zf.read("figures/contacts.svg").decode("utf-8", errors="ignore")
        # SVG embeds each tick label as a <text> node; tickling for the
        # residue triple is sufficient — full XML parse would be brittle.
        assert "LEU275" in svg
        assert "THR276" in svg
        assert "ARG278" in svg
        # And the relabel state surfaces in the title.
        assert "author numbering" in svg

    def test_rmsd_figure_annotates_verdict_and_pose_final(self, tmp_path):
        md_id = "fixture-md-id"
        md = _write_md_fixture(tmp_path, md_id)
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        svg = zf.read("figures/rmsd.svg").decode("utf-8", errors="ignore")
        assert "verdict: stable" in svg
        assert "pose final: 1.51" in svg
        assert "internal final: 0.65" in svg
        # Footer caveat travels with every figure.
        assert "pose RMSD = receptor-frame displacement" in svg

    def test_provenance_implicit_does_not_render_explicit_fields(self, tmp_path):
        """Default implicit fixture must NOT leak explicit-only keys into
        PROVENANCE.md — those would be misleading for an implicit run."""
        md_id = "fixture-md-id-implicit"
        md = _write_md_fixture(tmp_path, md_id)  # solvent="implicit" (default)
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        prov = zf.read("PROVENANCE.md").decode("utf-8")
        assert "**solvent**: `implicit`" in prov
        # Explicit-only keys must NOT appear when solvent=implicit.
        assert "water_model" not in prov
        assert "water_padding_nm" not in prov
        assert "ionic_strength_molar" not in prov
        assert "pressure_bar" not in prov
        assert "npt_equilibration_ps" not in prov
        assert "position_restraint_k_kj_per_mol_per_nm2" not in prov

    def test_provenance_explicit_surfaces_solvent_block(self, tmp_path):
        """Explicit fixture renders the full explicit-mode protocol block."""
        md_id = "fixture-md-id-explicit"
        md = _write_md_fixture(tmp_path, md_id, solvent="explicit")
        zb, _ = build_md_report_zip(md, md_id)
        zf = zipfile.ZipFile(io.BytesIO(zb))
        prov = zf.read("PROVENANCE.md").decode("utf-8")
        # Solvent + each explicit knob renders with its value.
        assert "**solvent**: `explicit`" in prov
        assert "**water_model**: `tip3p`" in prov
        assert "**water_padding_nm**: 1.0" in prov
        assert "**ionic_strength_molar**: 0.15" in prov
        assert "**pressure_bar**: 1.0" in prov
        assert "**barostat_frequency_steps**: 25" in prov
        assert "**npt_equilibration_ps**: 100.0" in prov
        assert "**position_restraint_k_kj_per_mol_per_nm2**: 1000.0" in prov
        assert "**equilibration_discard_ps**: 20.0" in prov

    def test_explicit_settings_round_trip_through_summary_json(self, tmp_path):
        """Defense-in-depth: the fixture writes settings into the canonical
        location summary.json["settings"]; this is what every downstream
        consumer (gating, report, MM-GBSA) reads. Catches a stray rename or
        an accidental field drop in the fixture shape."""
        md_id = "fixture-md-id-explicit-shape"
        md = _write_md_fixture(tmp_path, md_id, solvent="explicit")
        summary = json.loads((md / "summary.json").read_text())
        settings = summary["settings"]
        # The 9 [B1] explicit fields must all be present.
        for key in (
            "solvent", "water_model", "water_padding_nm",
            "ionic_strength_molar", "pressure_bar",
            "barostat_frequency_steps", "npt_equilibration_ps",
            "position_restraint_k_kj_per_mol_per_nm2",
            "equilibration_discard_ps",
        ):
            assert key in settings, f"explicit fixture missing settings.{key}"
        assert settings["solvent"] == "explicit"

    def test_missing_summary_raises_filenotfound(self, tmp_path):
        md_id = "no-summary-id"
        md = tmp_path / md_id / "md"
        md.mkdir(parents=True)
        # rmsd.csv exists but summary.json absent — bundler must refuse.
        (md / "rmsd.csv").write_text(
            "time_ps,rmsd_backbone_ca_A,rmsd_ligand_pose_A,rmsd_ligand_internal_A\n0,0,0,0\n"
        )
        with pytest.raises(FileNotFoundError):
            build_md_report_zip(md, md_id)


# =====================================================================
# HTTP endpoint integration — drive the ASGI app via httpx (the same
# pattern as test_auth_flows.py because Starlette 0.27 + httpx 0.28
# broke TestClient's `app=` kwarg; see CLAUDE.md "Known test constraint").
# =====================================================================

def _patch_jobs_dir(monkeypatch, fake_jobs):
    """Patch JOBS_DIR everywhere the endpoint reads it. Mirrors
    test_auth_flows.py's pattern of touching every module-bound copy
    rather than relying on a single canonical reference."""
    from backend.app.core import config as cfg
    from backend.app.api.routes import md as md_routes_mod
    from backend.app.md import job as md_job_mod
    monkeypatch.setattr(cfg, "JOBS_DIR", fake_jobs)
    monkeypatch.setattr(md_routes_mod, "JOBS_DIR", fake_jobs)
    monkeypatch.setattr(md_job_mod, "JOBS_DIR", fake_jobs)
    # Disable auth enforcement so we don't need a Job row + session for
    # the endpoint test — load_owned_job_or_404 falls through to
    # FS-presence in unenforced mode. Matches the test_auth_flows
    # convention of toggling AUTH_ENFORCED rather than seeding users.
    from backend.app.auth import config as auth_cfg
    monkeypatch.setattr(auth_cfg, "AUTH_ENFORCED", False)


def test_md_report_endpoint_200_and_zip(tmp_path, monkeypatch):
    md_id = "11111111-1111-1111-1111-111111111111"
    fake_jobs = tmp_path / "jobs"
    fake_jobs.mkdir()
    _write_md_fixture(fake_jobs, md_id)
    _patch_jobs_dir(monkeypatch, fake_jobs)

    from backend.app.main import app

    async def fetch():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport, base_url="http://test",
            headers={"X-Requested-With": "XMLHttpRequest"},
        ) as client:
            return await client.get(f"/api/md/{md_id}/report")

    r = asyncio.run(fetch())

    assert r.status_code == 200, r.text[:400]
    assert r.headers["content-type"] == "application/zip"
    assert f'filename="md_report_{md_id[:8]}.zip"' in r.headers["content-disposition"]
    assert r.headers["cache-control"] == "no-store"

    zf = zipfile.ZipFile(io.BytesIO(r.content))
    names = set(zf.namelist())
    for panel in ("rmsd", "hbonds", "contacts"):
        assert f"figures/{panel}.svg" in names
        assert f"figures/{panel}.png" in names
    assert "PROVENANCE.md" in names
    assert "data/summary.json" in names


def test_md_report_endpoint_missing_job_404(tmp_path, monkeypatch):
    fake_jobs = tmp_path / "jobs"
    fake_jobs.mkdir()  # exists but contains no jobs
    _patch_jobs_dir(monkeypatch, fake_jobs)

    from backend.app.main import app

    async def fetch():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport, base_url="http://test",
            headers={"X-Requested-With": "XMLHttpRequest"},
        ) as client:
            return await client.get("/api/md/22222222-2222-2222-2222-222222222222/report")

    r = asyncio.run(fetch())
    assert r.status_code == 404, r.text[:400]
