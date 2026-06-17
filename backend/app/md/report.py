"""Bundle a completed MD job's data + publication figures into a .zip.

Companion to `backend/app/md/figures.py` (which renders the panels) —
this module assembles the zip a single endpoint returns to the operator.

Contents of the zip (when all sources are present):

  figures/
    rmsd.svg, rmsd.png         — pose vs. internal RMSD time series
    hbonds.svg, hbonds.png     — H-bond count time series
    contacts.svg, contacts.png — top per-residue contacts (author-numbered
                                 when receptor_renumbered=True)
  data/
    rmsd.csv                   — raw time series (4 columns post-Q6b)
    hbonds.csv                 — H-bond count + persistence footer
    contacts.csv               — per-residue frequencies
    summary.json               — the full md summary artifact
  PROVENANCE.md                — provenance + caveats; the methods-section
                                 prose for any figure leaving this bundle.

Missing source files degrade gracefully (the matching artifact is
omitted) — the report still bundles whatever IS available. The
PROVENANCE.md always renders so a reader has the context to interpret
whatever did make it in.
"""

from __future__ import annotations

import io
import json
import os
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

from .figures import CAVEAT_FOOTER, render_all_md_figures


def _git_sha(repo_root: Path) -> Optional[str]:
    """Best-effort short SHA of the backend tree. Returns None on any
    failure (no git binary, repo missing, etc.) — provenance still
    renders, the SHA line just reads "unknown"."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(repo_root),
            capture_output=True, text=True, check=False, timeout=5,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip()
    except (FileNotFoundError, subprocess.SubprocessError, OSError):
        pass
    return None


def _render_provenance_md(
    summary: Dict,
    md_job_id: str,
    git_sha: Optional[str],
    fe_gate: Optional[Dict],
) -> str:
    """The README that travels with every report bundle.

    Lists every method/parameter a methods section would cite, the
    mechanics-fixture caveat, the pose-vs-internal RMSD definition,
    the engine kind + free-energy gating status, and the git SHA so the
    figures can be re-derived if the source data is rebuilt.
    """
    engine = (summary.get("engine") or {})
    metrics = (summary.get("metrics") or {})
    settings = (summary.get("settings") or {})
    rmsd_pose_final = metrics.get("rmsd_ligand_pose_final_a")
    rmsd_pose_max = metrics.get("rmsd_ligand_pose_max_a")
    rmsd_internal_final = metrics.get("rmsd_ligand_internal_final_a")
    rmsd_internal_max = metrics.get("rmsd_ligand_internal_max_a")
    rmsd_bb_final = metrics.get("rmsd_backbone_final_a")
    hbond_frac = metrics.get("hbond_persistence_frac")

    lines = []
    lines.append(f"# MD report — {md_job_id[:8]}…")
    lines.append("")
    lines.append(f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}.")
    lines.append("")
    lines.append("## Identity")
    lines.append("")
    lines.append(f"- **md_job_id**: `{md_job_id}`")
    lines.append(f"- **parent docking job**: `{summary.get('docking_job_id') or '—'}`")
    lines.append(f"- **ligand**: `{summary.get('ligand') or '—'}`  (pose_rank={summary.get('pose_rank')})")
    lines.append(f"- **SMILES**: `{summary.get('smiles') or '—'}`")
    lines.append(f"- **git SHA (backend)**: `{git_sha or 'unknown'}`")
    lines.append(f"- **summary schema**: `{summary.get('schema_version') or '—'}`")
    lines.append("")
    lines.append("## Engine + protocol")
    lines.append("")
    lines.append(f"- **engine kind**: `{engine.get('kind') or 'unknown'}`")
    if engine.get("attempts"):
        lines.append(f"- **engine attempts**: {engine.get('attempts')}")
    lines.append(f"- **production_ps**: {settings.get('production_ps')}")
    lines.append(f"- **snapshot_every_ps**: {settings.get('snapshot_every_ps')}")
    lines.append(f"- **temperature_k**: {settings.get('temperature_k')}")
    # Solvent disclosure — implicit (OBC2, NVT) vs explicit (TIP3P, PME,
    # staged NPT). Pre-[B1] summaries default to implicit.
    solvent_mode = settings.get("solvent", "implicit")
    lines.append(f"- **solvent**: `{solvent_mode}`")
    if solvent_mode == "explicit":
        lines.append(f"  - **water_model**: `{settings.get('water_model')}`")
        lines.append(f"  - **water_padding_nm**: {settings.get('water_padding_nm')}")
        lines.append(f"  - **ionic_strength_molar**: {settings.get('ionic_strength_molar')}")
        lines.append(f"  - **pressure_bar**: {settings.get('pressure_bar')}")
        lines.append(f"  - **barostat_frequency_steps**: "
                     f"{settings.get('barostat_frequency_steps')}")
        lines.append(f"  - **npt_equilibration_ps**: "
                     f"{settings.get('npt_equilibration_ps')}")
        lines.append(f"  - **position_restraint_k_kj_per_mol_per_nm2**: "
                     f"{settings.get('position_restraint_k_kj_per_mol_per_nm2')}")
        eq_discard = settings.get("equilibration_discard_ps")
        if eq_discard is not None:
            lines.append(f"  - **equilibration_discard_ps**: {eq_discard}")
    lines.append(f"- **n_frames**: {summary.get('n_frames')}")
    lines.append(f"- **wall_seconds**: {summary.get('wall_seconds')}")
    lines.append("")
    lines.append("## Verdict + metrics")
    lines.append("")
    lines.append(f"- **verdict**: `{summary.get('verdict') or '—'}`")
    if summary.get("rationale"):
        lines.append(f"- **rationale**: {summary['rationale']}")
    lines.append(f"- **rmsd_ligand_pose_final_a**:     {_fmt(rmsd_pose_final)}  Å  ← primary, gates the verdict")
    lines.append(f"- **rmsd_ligand_pose_max_a**:       {_fmt(rmsd_pose_max)}  Å")
    lines.append(f"- **rmsd_ligand_internal_final_a**: {_fmt(rmsd_internal_final)}  Å  ← diagnostic, pre-Q6b semantics")
    lines.append(f"- **rmsd_ligand_internal_max_a**:   {_fmt(rmsd_internal_max)}  Å")
    lines.append(f"- **rmsd_backbone_final_a**:        {_fmt(rmsd_bb_final)}  Å")
    lines.append(f"- **hbond_persistence_frac**:       {_fmt(hbond_frac)}")
    lines.append(f"- **receptor_renumbered**:          {summary.get('receptor_renumbered')}  (Q6b PART 2 author-numbering relabel)")
    lines.append("")
    lines.append("## Free-energy gating status")
    lines.append("")
    if fe_gate is None:
        lines.append("- (not computed — see /api/free-energy/gating for the live decision)")
    elif fe_gate.get("can_run"):
        lines.append("- **passes hard gates** (verdict==stable, engine==openmm_full, frames sufficient).")
        warns = fe_gate.get("warnings") or []
        if warns:
            lines.append(f"- soft-gate warnings: {[w.get('key') for w in warns]}")
    else:
        reasons = fe_gate.get("reasons") or []
        lines.append("- **blocked.** Hard-gate reasons:")
        for r in reasons:
            lines.append(f"  - `{r.get('key')}` — {r.get('message')}")
    lines.append("")

    # Free-energy estimate (when present) — surface the actual ΔG number,
    # not just the gate decision. The two are independent: gate-blocked
    # estimates ARE computed but flagged preliminary, and a passing-gate
    # estimate that's somehow not in summary still has to surface as
    # "not yet computed". This block reads summary.free_energy directly;
    # if the block is missing or status=='planned' we say so verbatim.
    fe = summary.get("free_energy") or {}
    fe_status = fe.get("status")
    lines.append("## Free-energy estimate (MM-GBSA, single-trajectory)")
    lines.append("")
    if fe_status != "completed":
        lines.append(f"- **status**: `{fe_status or 'absent'}` — no ΔG_bind to report.")
    else:
        method = fe.get("method") or {}
        result = fe.get("result") or {}
        prov = fe.get("provenance") or {}
        comps = result.get("components_mean_kcal_per_mol") or {}
        comps_sem = result.get("components_sem_kcal_per_mol") or {}
        preliminary = bool(result.get("preliminary"))
        lines.append(f"- **status**: `completed`{'  (PRELIMINARY — see gate reason)' if preliminary else ''}")
        dg = result.get("delta_g_mean_kcal_per_mol")
        sem = result.get("delta_g_sem_kcal_per_mol")
        std = result.get("delta_g_stddev_kcal_per_mol")
        lines.append(f"- **ΔG_bind**:    {_fmt(dg)} ± {_fmt(sem)} kcal/mol   (stddev {_fmt(std)})")
        lines.append(f"- **components (mean ± SEM, kcal/mol):**")
        for k in ("bonded", "nonbonded", "solvation"):
            lines.append(f"  - {k:<10s}: {_fmt(comps.get(k))} ± {_fmt(comps_sem.get(k))}")
        if preliminary and result.get("gate_reason"):
            lines.append(f"- **preliminary because**: {result['gate_reason']}")
        lines.append(f"- **method**: {method.get('name', '?')} · {method.get('implicit_solvent', '?')} · "
                     f"{method.get('small_molecule_forcefield', '?')}")
        if method.get("configurational_entropy"):
            lines.append(f"- **configurational entropy**: {method['configurational_entropy']}")
        lines.append(f"- **frames used**: {prov.get('n_frames_used')}  (skipped {prov.get('n_frames_skipped', 0)}; "
                     f"equilibration discard {prov.get('equilibration_discard_frames', 0)})")
        if prov.get("wall_seconds") is not None:
            lines.append(f"- **wall_seconds (estimator)**: {prov['wall_seconds']}")
    lines.append("")
    lines.append("## RMSD definitions (post-Q6b, 2026-06-04)")
    lines.append("")
    lines.append("- **Pose RMSD** = receptor-frame ligand displacement. For each "
                 "frame, superpose its backbone Cα onto the reference's via "
                 "Kabsch, apply that same (R, t) to the ligand heavy atoms, "
                 "and report plain RMSD vs. the reference-pose ligand. This is "
                 "the metric that actually answers \"did the ligand stay in "
                 "the pocket?\". It is what the verdict is keyed on.")
    lines.append("- **Internal RMSD** = ligand-on-ligand Kabsch RMSD — the "
                 "pre-Q6b metric, kept as a diagnostic. Captures conformational "
                 "change of the ligand itself (sliding torsions, bond rotations) "
                 "but is blind to pocket displacement.")
    lines.append("- **Verdict bands**: ≤ 2 Å stable, 2–4 Å drifting, > 4 Å "
                 "unstable. Bands are MD-convention defaults; a calibration "
                 "pass against a labeled bound/unbound set is the right "
                 "follow-up.")
    lines.append("")
    lines.append("## Author numbering (Q6b PART 2)")
    lines.append("")
    if summary.get("receptor_renumbered"):
        lines.append("Contacts panel labels are in the docking receptor's "
                     "**author numbering** (e.g., LEU275/THR276/ARG278 on "
                     "β-tubulin's taxane site). The MD prep renumbers chain B "
                     "to 1..N contiguous; the analysis-time relabel maps it "
                     "back to author numbering by per-chain sequence position. "
                     "See `backend/app/md/analyze.py:build_md_to_docking_resseq_map`.")
    else:
        lines.append("The contacts panel labels are in MD numbering (no "
                     "relabel ran — typically because the docking-receptor "
                     "and MD-frame chain lengths or residue names disagreed; "
                     "see Q6b notes for the per-chain safety rule).")
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append(f"- {CAVEAT_FOOTER}")
    lines.append("- Free-energy estimates in this report are **single-trajectory "
                 "MM-GBSA** (no explicit ligand-alone or receptor-alone MD, no "
                 "MM/PBSA, no entropy). The implicit-solvent model is OBC2, "
                 "consistent with the MD that produced the frames. The number "
                 "is preliminary whenever the gating reasons block above lists "
                 "`reasons` — most commonly INSUFFICIENT_FRAMES below the "
                 "200-frame minimum.")
    lines.append("- All thresholds, engine attempts, and the renumber-map "
                 "policy are self-describing in `data/summary.json` and "
                 "`data/rmsd.csv` so this report can be re-derived from raw "
                 "artifacts without the figures module.")
    return "\n".join(lines) + "\n"


def _fmt(v) -> str:
    if v is None:
        return "—"
    try:
        f = float(v)
        return f"{f:.3f}"
    except (TypeError, ValueError):
        return str(v)


def build_md_report_zip(
    md_dir: Path,
    md_job_id: str,
    *,
    repo_root: Optional[Path] = None,
    fe_gate: Optional[Dict] = None,
) -> Tuple[bytes, str]:
    """Render the full report bundle. Returns (zip_bytes, filename).

    Raises FileNotFoundError if md_dir or summary.json is missing — the
    route layer catches this and returns 404.
    """
    summary_path = md_dir / "summary.json"
    if not md_dir.is_dir() or not summary_path.is_file():
        raise FileNotFoundError(
            f"MD job {md_job_id!r} has no md/summary.json at {summary_path}"
        )

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    figures = render_all_md_figures(md_dir, summary)

    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[3]  # backend/app/md → repo root
    git_sha = _git_sha(repo_root)

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        # Figures — SVG + PNG per panel.
        for key, (svg, png) in figures.items():
            zf.writestr(f"figures/{key}.svg", svg)
            zf.writestr(f"figures/{key}.png", png)

        # Raw data: copy bytes verbatim so the operator can re-render
        # off the same source the figures were built from.
        for fname in ("rmsd.csv", "hbonds.csv", "contacts.csv", "summary.json"):
            p = md_dir / fname
            if p.is_file():
                zf.writestr(f"data/{fname}", p.read_bytes())

        zf.writestr(
            "PROVENANCE.md",
            _render_provenance_md(summary, md_job_id, git_sha, fe_gate),
        )

    filename = f"md_report_{md_job_id[:8]}.zip"
    return buf.getvalue(), filename
