# All MD diagnostics — every run ever queued (READ-ONLY)

Sweep of the full MD catalogue (DB md_runs ∪ jobs/<id>/md/) at this branch state. 33 MD jobs total. Pocket residue set derived ONCE from the canonical reference run (5bc61f59-834f-4e71-a492-d32ddfdc7326, explicit taxol pose 0) frame 0 — 25 residues within 5 Å of any ligand heavy atom — and intersected with each run's own receptor for the deep per-run analysis.

Honest reporting: every run is listed, including failed / empty / surrogate / smoke-test / wrong-pocket. Caveat per row in the Notes column. `Converged?` = yes when |1st-half − 2nd-half mean| ≤ 0.3 Å AND |last-half slope| ≤ 1e-3 Å/ps; ≈ if borderline; no otherwise.

**FE gate verdicts are deliberately NOT applied here.** Just the raw metrics, grouped for readability.

## Table

| job_id (8) | compound | engine | sol | dur(ps) | frames | pose-pkt fin/max | COM fin/max | conv | hb/contact persist | in-pkt fin | buried fin | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `183b71bc` | taxol | surrogate | explicit | 1000.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `1f01da83` | taxol | openmm_full | implicit | 500.00 | 101 | 1.57 / 2.15 | 0.94 / 1.47 | yes | 31%  (of 60 init contacts) | Y | 0.79 | plateaued |
| `2a854872` | taxol | surrogate | explicit | 10.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction); production_ps=10.0 |
| `3fac203e` | taxol | surrogate | explicit | 1000.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `5817afff` | taxol | openmm_full | implicit | 1000.00 | 201 | 4.28 / 4.80 | 3.32 / 3.59 | no | 19%  (of 80 init contacts) | Y | 0.60 | still drifting at end |
| `5bc61f59` | taxol | openmm_full | explicit | 1000.00 | 201 | 2.36 / 2.89 | 1.29 / 2.25 | yes | 23%  (of 85 init contacts) | Y | 0.79 | plateaued |
| `613a0fb0` | taxol | surrogate | explicit | 1000.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `6b581508` | taxol | surrogate | implicit | 10.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction); production_ps=10.0 |
| `6cfae070` | taxol | openmm_full | explicit | 10.00 | 11 | 0.96 / 1.18 | 0.44 / 0.69 | yes | 53%  (of 91 init contacts) | Y | 0.82 | smoke-test scale; production_ps=10.0; plateaued |
| `6d4115e4` | taxol | surrogate | implicit | 10.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction); production_ps=10.0 |
| `ac0a7952` | taxol | surrogate | explicit | 10.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction); production_ps=10.0 |
| `dba42e7d` | taxol | openmm_full | implicit | 10.00 | 11 | 0.79 / 1.03 | 0.41 / 0.76 | no | 45%  (of 67 init contacts) | Y | 0.65 | smoke-test scale; production_ps=10.0; still drifting at end |
| `e06fdb95` | taxol | surrogate | implicit | 10.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction); production_ps=10.0 |
| `e19003b6` | taxol | surrogate | explicit | 1000.00 | 41 | — / — | — / — | — | hbond 100% | — | — | RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `f9c29ee9` | taxol | openmm_full | implicit | 10.00 | 11 | 1.88 / 2.07 | 1.12 / 1.20 | no | 25%  (of 42 init contacts) | Y | 0.65 | smoke-test scale; production_ps=10.0; still drifting at end |
| `34840aa1` | taxol | openmm_full | explicit | 1000.00 | 201 | 1.68 / 2.30 | 0.85 / 1.87 | ≈ | 28%  (of 115 init contacts) | Y | 0.61 | — |
| `80e53d8a` | taxol | openmm_full | explicit | 1000.00 | 201 | 1.96 / 2.43 | 1.01 / 1.82 | ≈ | 20%  (of 102 init contacts) | Y | 0.73 | — |
| `a0b04941` | Juliprosopine | openmm_full | explicit | 1000.00 | 201 | 4.72 / 5.24 | 2.13 / 2.86 | no | 0%  (of 42 init contacts) | Y | 0.95 | still drifting at end |
| `00dda37f` | Primaquine | openmm_full | explicit | 1000.00 | 201 | 2.85 / 4.10 | 1.66 / 2.37 | yes | 8%  (of 42 init contacts) | Y | 0.84 | plateaued |
| `11fa6090` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `27467317` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `524febf8` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `592ef3e0` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `6ef436cd` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `86e5195b` | taxol (wrong-pocket) | openmm_full | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started) |
| `987322c1` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `c2284960` | taxol (wrong-pocket) | surrogate | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started); RDKit surrogate (rigid receptor; pose-RMSD ≡ internal-RMSD by construction) |
| `d34b991f` | taxol (wrong-pocket) | openmm_full | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started) |
| `ff0484cb` | taxol (wrong-pocket) | — | implicit | — | 0 | — / — | — / — | — | — | — | — | WRONG POCKET (d96d719a lineage; do not use for binding analysis); no summary.json (DB-only / never started) |
| `9fc6ad18` | taxol | — | implicit | 1000.00 | 0 | — / — | — / — | — | — | — | — | no summary.json (DB-only / never started) |
| `a691297a` | taxol | — | implicit | 10.00 | 0 | — / — | — / — | — | — | — | — | no summary.json (DB-only / never started); production_ps=10.0 |
| `ba9e876e` | taxol | — | implicit | 10.00 | 0 | — / — | — / — | — | — | — | — | no summary.json (DB-only / never started); production_ps=10.0 |
| `d2f06998` | taxol | — | implicit | 10.00 | 0 | — / — | — / — | — | — | — | — | no summary.json (DB-only / never started); production_ps=10.0 |

## Columns

- **engine**: `openmm_full` (real MD), `surrogate` (RDKit MMFF wiggle in rigid receptor — pose-RMSD ≡ internal-RMSD by construction), or `—` (DB row only, no run output).
- **sol**: implicit (OBC2) or explicit (TIP3P NPT octahedron per the [B7] / [B9] protocol).
- **pose-pkt fin/max**: ligand pose RMSD in pocket-Cα-aligned frame (5 Å pocket Cα fit), final + trajectory max (Å). The primary binding-mode drift metric.
- **COM fin/max**: ligand center-of-mass displacement after pocket alignment, final + max (Å).
- **conv**: yes / ≈ / no per the convergence rule above.
- **hb/contact persist**: for openmm_full runs with frames, the fraction of the frame-0 (ligand × protein heavy atom ≤ 4 Å) contacts retained per frame, mean over the equilibrated window (last 30 %). For surrogate / short runs, the engine's own `hbond_persistence_frac` is reported instead (different metric).
- **in-pkt fin**: Y if the ligand COM in the last frame is < 5 Å from frame-0 COM after pocket alignment.
- **buried fin**: fraction of ligand heavy atoms with any protein heavy atom within 4.5 Å, in the last frame.
- **notes**: short caveat list; full per-row JSON in `docs/all_md_diagnostics.csv`.

## Grouping

Order: taxol on the verified taxane pocket (4a37bb0c), grouped by pose_rank; controls (Juliprosopine, Primaquine); other / unidentified runs; wrong-pocket lineage (d96d719a — receptor is a different construct, do NOT use for any binding analysis); DB-only rows (md_run entry exists but no summary.json on disk; typically a queued job that crashed before any output).