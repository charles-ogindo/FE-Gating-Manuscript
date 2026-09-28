# Fourteen-receptor benchmark

Results supporting Section 9 of the manuscript. The workflow was run on
fourteen protein-ligand series from the OpenFF Protein-Ligand Benchmark
collection, comprising 359 compounds and 443 production trajectories.

## Layout

    screening/<receptor>.jsonl     one line per screening run
    calibration/<receptor>.json    the core contacts derived for that receptor
    stage9_threshold_sweep.json    the retention threshold sweep

## Where each result comes from

| Paper item | File and fields |
|---|---|
| Table 1, calibration | `calibration/*.json`: `calibration_ligand`, `core_set_size`, `core_retentions` |
| Table 2, ranking performance | `screening/*.jsonl`: `dG_pred_kcal`, `dG_exp_kcal`, `core_gated` |
| Table 3 and Figure 2, threshold | `stage9_threshold_sweep.json` |
| Figure 1, ranking across receptors | `../../figures/figures_pgf.py` |
| Section 9.4, the two TNKS2 cases | `screening/tnks2.jsonl`: compounds `lig_1a` and `lig_7` |
| Section 9.6, PDE2 | `screening/pde2.jsonl` |

## Screening record fields

Each line is one molecular dynamics run with its endpoint estimate.

    core_retention    mean fraction of core contacts held over the converged window
    landed            whether the pose RMSD settled
    core_gated        the gate's verdict
    dG_pred_kcal      the endpoint estimate
    dG_exp_kcal       the measured affinity
    n_eff             independent samples in the converged window

## The threshold sweep

`../../scripts/stage9_threshold_sweep.py` produced `stage9_threshold_sweep.json`.
It re-derives each receptor's core at 0.60, 0.70 and 0.80 from the calibration
replicates, re-applies the retention condition to every screening run, and
recomputes the correlations. It runs no simulation, though it parses the stored
trajectory frames, so it requires the full application and the trajectory
archive. It is included as the record of how the reported sweep was performed.
