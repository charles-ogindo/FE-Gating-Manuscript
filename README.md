# free-energy-gating

Reviewer / reproduction repository for the **stability- and free-energy-gated
docking → MD → MM-GBSA** pipeline described in the accompanying manuscript
(β-tubulin / taxane-site case study).

This repo contains **only the evidence needed to read, run, and verify the
scientific logic** — the gating criteria, the MD setup, the MM-GBSA scoring,
the convergence/stability diagnostics, the replicate protocol, and a hermetic
test suite. The full interactive application (web frontend, Mol\* viewer, API
server, auth, job database, deployment) is **deliberately excluded** — it is
product, not evidence.

---

## What's here

```
backend/app/
  free_energy/        gating criteria, MM-GBSA scoring, convergence sampling, artifact schema
    gating.py           stability + free-energy gate logic (validate(), GatingResult)
    sampling.py         convergence / sampling-plan logic
    mmgbsa.py           MM-GBSA energy decomposition
    mmgbsa_runner.py    runs MM-GBSA over MD frames
    protocol.py         end-to-end FE protocol descriptor
    __init__.py         ARTIFACT_SCHEMA_VERSION + shared dataclasses
  md/
    engine_openmm.py    explicit-solvent MD engine (solvation, force field, box shape, chain truncation)
    engine_surrogate.py RDKit fallback engine (used by tests; no GPU/MD needed)
    analyze.py          RMSD / H-bond / contact analysis + classify_stability (stability gate)
    io.py               docked-pose resolution (the docking → MD link)
    report.py           MD report bundle
    figures.py          MD figures
    job.py              run_md_job orchestrator (workflow engine: docking → MD → analysis → scoring)
    __init__.py         ARTIFACT_SCHEMA_VERSION + verdict constants (stable/drifting/unstable/failed)
  docking/
    md_receptor_prep.py MD-grade receptor prep (PDBFixer, cofactors, ion templates, E-hook truncation)
  utils/
    rmsd.py             Kabsch RMSD
  core/
    config.py           TRIMMED shim — provides only JOBS_DIR (see note below)
backend/tests/          hermetic test suite (synthetic trajectories; no DB, no network, no GPU)
scripts/
  run_replicates_taxol.py     the replicate MD sweep (corrected octahedron + E-hook truncation protocol)
  aggregate_replicates_taxol.py replicate aggregator
  stability_diagnostics.py    stability diagnostics
  all_md_diagnostics.py       per-run MD diagnostics
  corrected_fe_gate.py        corrected free-energy gate
  fe_gate_all_runs.py         free-energy gate across runs
  run_fe_qualifying.py        MM-GBSA on gate-qualifying runs
  pocket_aligned_rmsd.py      pocket-frame RMSD diagnostic
docs/
  meeko_env.lock.txt    exact conda environment (conda list --explicit)
examples/beta_tubulin/  minimal example inputs (SMILES, docking box, MD settings; receptor = PDB 1JFF)
```

### Deliberately excluded (product, not evidence)
Web frontend & Mol\* viewer · FastAPI server & API routes · authentication/session code ·
the Postgres job catalog (`db/`, migrations, `repository.py`) and its **contents** ·
deployment / `docker-compose` / server config · internal backfill & admin tooling ·
credentials · large trajectory files and job artifacts.

---

## 1. Recreate the environment

The scientific stack is conda-managed. Recreate it from the exact lockfile:

```bash
conda create --name fe-gating --file docs/meeko_env.lock.txt
conda activate fe-gating
```

(The lockfile is `conda list --explicit` from the original `meeko_env`:
Python 3.10, OpenMM 8.5, AmberTools 23.6, openmmforcefields, OpenFF toolkit,
RDKit, PDBFixer, NumPy/SciPy, pytest, etc.)

## 2. Verify the gating logic (hermetic — this is the runnable verification)

The test suite runs against **synthetic trajectories and fixtures** — no GPU, no
database, no network. Run it from the repo root with the repo on `PYTHONPATH`:

```bash
PYTHONPATH=. python -m pytest backend/tests/ -v
```

What it covers:

| Test | Verifies |
|---|---|
| `test_free_energy_gating.py` | free-energy gate `validate()` + convergence `compute_sampling_plan()` |
| `test_mmgbsa.py`, `test_mmgbsa_strip.py` | MM-GBSA scoring + frame-position collection / solvent stripping |
| `test_engine_truncate.py` | chain-range truncation + OXT capping (the A439–451 E-hook cut) |
| `test_md_analyze_waters.py`, `test_q6b_pose_rmsd.py` | MD frame parsing, pose-frame RMSD series |
| `test_md_report.py` | MD report bundling |
| `test_rmsd.py` | Kabsch RMSD |

## 3. Reproduce the science (docking → MD → MM-GBSA)

Inputs for the worked example are in `examples/beta_tubulin/` (fetch receptor
PDB **1JFF** yourself — see that folder's README). The protocol:

1. **Dock** the taxane-site box (`examples/beta_tubulin/box.json`) against 1JFF.
2. **MD** each top pose with the explicit-solvent settings in
   `examples/beta_tubulin/md_settings.json` (octahedron TIP3P box, 1.0 nm
   padding, 0.15 M, **A439–451 α-tubulin E-hook truncation**). The replicate
   sweep driver is:
   ```bash
   PYTHONPATH=. python scripts/run_replicates_taxol.py
   ```
   3 poses × 3 seeded replicates; `random_seed` is the only per-replicate
   difference. Results are aggregated by:
   ```bash
   PYTHONPATH=. python scripts/aggregate_replicates_taxol.py
   ```
3. **Gate + score**: stability gate (`md/analyze.py::classify_stability`),
   free-energy gate (`free_energy/gating.py::validate`), then MM-GBSA on
   qualifying runs:
   ```bash
   PYTHONPATH=. python scripts/run_fe_qualifying.py
   PYTHONPATH=. python scripts/fe_gate_all_runs.py
   ```

> **Runnable boundary.** The **hermetic test suite (step 2) runs standalone** in
> this repo and is the self-contained verification of the gating/scoring logic.
> The **orchestration scripts in step 3** (`run_replicates_taxol.py`,
> `aggregate_*`, `*_fe_gate*`, diagnostics) and `md/job.py` reference the full
> application's **job database layer** (`backend.app.db`, `core.job_artifacts`),
> which is part of the private product and not included here. They are included
> as the **authoritative, unedited record of the exact protocol and commands**;
> running them end-to-end against real jobs requires the full application.

## Notes

- **`backend/app/core/config.py` is a trimmed shim.** The real config (scenario
  rules, safety weights, database URL, CORS) is product config. The copied
  scientific modules only need `JOBS_DIR`, so the shim provides just that.
- **No receptor or trajectory files are committed.** Fetch 1JFF from the RCSB
  (see `examples/beta_tubulin/README.md`); `.gitignore` excludes structures,
  trajectories, databases, and job artifacts.
- **Artifact schema** versions live in `backend/app/md/__init__.py`
  (`ARTIFACT_SCHEMA_VERSION`, verdict labels) and
  `backend/app/free_energy/__init__.py`.
- **License:** see `LICENSE` — **TODO**, choose before public release.
