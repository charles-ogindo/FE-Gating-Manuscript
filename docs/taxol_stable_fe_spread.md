# Taxol stable-set operative-window ΔG — between-replicate spread

MM-GBSA operative window (t0_energy→end; NOT last-30%) over the **20 stable** taxol 1 ns
runs (Stage-1c structural gate). `pose0_rep5` lands very late (t0_energy 0.90 ns → 21-frame
window, **N_eff 8.8 < floor 10**) and is the only stable run below the sampling-sufficiency
minimum. This report gives the spread WITH and WITHOUT it. **Provisional (1 ns < 3 ns); read-only; not committed.**

## Per-pose ΔG (kcal/mol)

| set | pose | n | mean ΔG | SD | SEM | t0_energy range |
|---|---|---|---|---|---|---|
| all 20 | 0 | 5 | -52.04 | 6.36 | 2.84 | 0.21–0.90 |
| all 20 | 1 | 8 | -51.55 | 4.94 | 1.75 | 0.13–0.80 |
| all 20 | 2 | 7 | -48.71 | 4.32 | 1.63 | 0.04–0.73 (+1 no-land) |
| excl pose0_rep5 (n=19) | 0 | 4 | -51.56 | 7.24 | 3.62 | 0.21–0.59 |
| excl pose0_rep5 (n=19) | 1 | 8 | -51.55 | 4.94 | 1.75 | 0.13–0.80 |
| excl pose0_rep5 (n=19) | 2 | 7 | -48.71 | 4.32 | 1.63 | 0.04–0.73 (+1 no-land) |

## Spread summary

| metric | all 20 | excl pose0_rep5 (n=19) |
|---|---|---|
| overall mean ΔG | -50.68 | -50.50 |
| overall SD | 5.07 | 5.14 |
| overall SEM | 1.13 | 1.18 |
| **pooled within-pose SD** | **5.12** | **5.25** |
| within-run median cSEM | 0.606 | 0.571 |
| **between-replicate SD ÷ within-run median cSEM** | **8.4×** | **9.2×** |

**Defensible taxol estimate (excl. under-sampled pose0_rep5): ΔG = -50.5 ± 5.2 kcal/mol** (between-replicate SD), SEM 1.18 on the 19-run mean.

Excluding the one below-floor run barely moves the estimate (mean -50.7→-50.5, pooled SD 5.1→5.2), confirming the result is not driven by the under-sampled window. The single-run corrected SEM still understates the true uncertainty by ~9× — the real error bar is the between-replicate σ, not the intra-trajectory SEM.

## Per-run detail (all 20; ⚠ = below N_eff floor 10)

| run | pose | t0_energy (ns) | window frames | ΔG | cSEM | N_eff |
|---|---|---|---|---|---|---|
| pose0_rep1 | 0 | 0.29 | 143 | -56.38 | 0.403 | 109.4 |
| pose0_rep2 | 0 | 0.21 | 158 | -51.05 | 0.353 | 104.1 |
| pose0_rep3 | 0 | 0.59 | 83 | -41.52 | 0.420 | 58.1 |
| pose0_rep4 | 0 | 0.38 | 125 | -57.28 | 0.567 | 48.1 |
| pose0_rep5 ⚠ | 0 | 0.90 | 21 | -53.95 | 1.393 | 8.8 |
| pose1_rep1 | 1 | 0.60 | 81 | -50.91 | 0.466 | 50.2 |
| pose1_rep2 | 1 | 0.13 | 175 | -50.87 | 1.112 | 15.0 |
| pose1_rep3 | 1 | 0.17 | 168 | -48.87 | 0.973 | 14.8 |
| pose1_rep4 | 1 | 0.15 | 171 | -60.38 | 0.571 | 42.7 |
| pose1_rep5 | 1 | 0.80 | 42 | -55.43 | 0.641 | 38.7 |
| pose1_rep6 | 1 | 0.28 | 145 | -47.46 | 0.493 | 45.8 |
| pose1_rep7 | 1 | 0.76 | 49 | -53.86 | 1.096 | 15.9 |
| pose1_rep8 | 1 | 0.44 | 113 | -44.62 | 0.839 | 20.0 |
| pose2_rep1 | 2 | 0.73 | 55 | -44.46 | 0.426 | 49.7 |
| pose2_rep2 | 2 | 0.56 | 88 | -49.99 | 0.773 | 16.8 |
| pose2_rep4 | 2 | 0.04 | 192 | -52.02 | 0.482 | 50.3 |
| pose2_rep5 | 2 | 0.69 | 64 | -47.63 | 1.466 | 11.1 |
| pose2_rep6 | 2 | — (no land) | 201 | -42.09 | 1.072 | 16.7 |
| pose2_rep7 | 2 | 0.23 | 155 | -50.19 | 0.317 | 118.1 |
| pose2_rep8 | 2 | 0.68 | 66 | -54.57 | 0.886 | 18.3 |
