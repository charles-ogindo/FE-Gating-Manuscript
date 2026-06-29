# Replicate sweep aggregation — 3 reps × 3 taxol poses

Per-replicate MM-GBSA ΔG_bind (kcal/mol) aggregated into between-
replicate σ + pooled bootstrap 95% CI. Replaces the misleading intra-
run SEM from the single-replicate sweep (`docs/free_energy_qualifying.md`).

## Per-replicate table

| Pose | Rep | seed | md_id | gate | Δ(begin−end) | top-5 % | ΔG_bind | σ_intra | n_FE |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 20260617 | `6d53c1af` | ✅ | 0.08 Å | 100.0% | **-56.51** | 4.49 | — |
| 0 | 2 | 20260618 | `127ee13d` | ✅ | 0.07 Å | 99.7% | **-51.77** | 3.57 | — |
| 0 | 3 | 20260619 | `dbd62fa3` | ✅ | 0.22 Å | 98.3% | **-41.65** | 3.07 | — |
| 0 | 4 | 20260620 | `c2380fc0` | ✅ | 0.07 Å | 100.0% | **-58.28** | 4.07 | — |
| 0 | 5 | 20260621 | `deaec49f` | ✅ | 0.03 Å | 99.7% | **-51.26** | 4.20 | — |
| 0 | 6 | 20260622 | `fad90b00` | ❌ | 0.67 Å | 99.3% | **-55.08** | 3.83 | — |
| 0 | 7 | 20260623 | `65741d10` | ❌ | 1.46 Å | 99.3% | **-39.27** | 3.23 | — |
| 0 | 8 | 20260624 | `e45b2dd3` | ❌ | 0.68 Å | 98.3% | **-45.44** | 5.23 | — |
| 1 | 1 | 20261617 | `3c4d38e6` | ✅ | 0.32 Å | 100.0% | **-50.97** | 3.19 | — |
| 1 | 2 | 20261618 | `4aaa8ab8` | ✅ | 0.14 Å | 100.0% | **-50.47** | 3.43 | — |
| 1 | 3 | 20261619 | `9ac9849a` | ✅ | 0.02 Å | 99.7% | **-50.53** | 2.85 | — |
| 1 | 4 | 20261620 | `09b77a40` | ✅ | 0.10 Å | 100.0% | **-58.98** | 3.18 | — |
| 1 | 5 | 20261621 | `536ba129` | ✅ | 0.04 Å | 100.0% | **-57.62** | 5.16 | — |
| 1 | 6 | 20261622 | `d596130c` | ✅ | 0.24 Å | 100.0% | **-47.68** | 3.03 | — |
| 1 | 7 | 20261623 | `849118c3` | ✅ | 0.06 Å | 100.0% | **-54.07** | 4.09 | — |
| 1 | 8 | 20261624 | `289b578e` | ✅ | 0.12 Å | 100.0% | **-44.97** | 2.99 | — |
| 2 | 1 | 20262617 | `f95fe87d` | ✅ | 0.11 Å | 99.7% | **-45.00** | 3.45 | — |
| 2 | 2 | 20262618 | `c366e126` | ✅ | 0.02 Å | 99.7% | **-50.55** | 3.09 | — |
| 2 | 3 | 20262619 | `7b459365` | ✅ | 0.35 Å | 99.7% | **-45.10** | 3.44 | — |
| 2 | 4 | 20262620 | `4e12de8c` | ✅ | 0.11 Å | 100.0% | **-51.16** | 2.79 | — |
| 2 | 5 | 20262621 | `1a820507` | ✅ | 0.02 Å | 99.3% | **-47.66** | 4.92 | — |
| 2 | 6 | 20262622 | `ebe9f241` | ✅ | 0.18 Å | 98.7% | **-41.46** | 4.00 | — |
| 2 | 7 | 20262623 | `34c5cbdf` | ❌ | 0.41 Å | 100.0% | **-49.52** | 3.09 | — |
| 2 | 8 | 20262624 | `61af92c2` | ✅ | 0.17 Å | 98.3% | **-54.53** | 3.87 | — |

## Per-pose aggregation

| Pose | n_replicates | mean ΔG_bind | σ_between_replicates |
|---|---|---|---|
| 0 | 8 | **-49.91 kcal/mol** | 7.04 kcal/mol |
| 1 | 8 | **-51.91 kcal/mol** | 4.75 kcal/mol |
| 2 | 8 | **-48.12 kcal/mol** | 4.17 kcal/mol |

## Pooled bootstrap — 9 taxol replicates

Pooled mean ΔG_bind (n=24 replicates): **-49.98 kcal/mol**
95% bootstrap CI (n_resamples=10000, seed=20260616): [**-52.10**, **-47.81**] kcal/mol

## Interpretation

- **σ_between_replicates** is the headline replication uncertainty — replaces the misleading intra-run SEM (which only captures the per-frame variance within a single trajectory and dramatically understates the true uncertainty of single-trajectory MM-GBSA).
- **Pooled bootstrap 95% CI** is the supportable confidence interval for paclitaxel's ΔG_bind on this pocket under this force-field stack, over the population of 9 independent 1-ns explicit-TIP3P replicates across 3 docked poses.
- The single-replicate FE numbers in `docs/free_energy_qualifying.md` (5bc61f59 / 34840aa1 / 80e53d8a) are NOT in this aggregation — they're separate reference runs that established the qualifying baseline. The 9 replicates here are deliberately independent draws from the same protocol.