# MM-GBSA on FE-gate-qualifying MD runs

Single-trajectory MM-GBSA computed on the four explicit-TIP3P MD runs that cleared the final FE gate (eq-window begin/end Δ < 0.4 Å, in-pocket ≥ 0.95, top-5 residue persistence ≥ 65 %, ≥ 50 eq frames). Estimator: `backend.app.free_energy.mmgbsa.estimate_mmgbsa` (single-trajectory; bonded + nonbonded + GBSA-OBC2 solvation; configurational entropy omitted). Implicit-lineage run `02f30602` deliberately NOT included.

**Caveat — bundled gate vs corrected gate.** The bundled `free_energy.gating.validate` keys on `summary.verdict='stable'` and marks all four runs preliminary (verdict='drifting'). The corrected gate at `docs/all_md_corrected_gate.md` is the authoritative qualification source — it ignores `verdict` and judges convergence + current state. The MM-GBSA number is computed either way; `sampling_adequate=false / preliminary=true` in the per-run JSON is the bundled-gate's verdict-based flag, not a corrected-gate signal.

## Per-run results

| Run | Compound | Pose | Sol | ΔG_bind | bonded / nonbonded / solvation | sampling_adequate | wall (s) |
|---|---|---|---|---|---|---|---|
| `5bc61f59` | taxol | 0 | explicit | — | — | ⚠ preliminary | 629 |
| `34840aa1` | taxol | 1 | explicit | — | — | ⚠ preliminary | 481 |
| `80e53d8a` | taxol | 2 | explicit | — | — | ⚠ preliminary | 500 |
| `a0b04941` | Juliprosopine | 0 | explicit | — | — | ⚠ preliminary | 461 |
| `02f30602` | taxol (Run A extended) | 0 | implicit | — | — | ⚠ preliminary | 631 |

## Method

- Force-field stack: amber14-all + amber14/tip3p ion templates + implicit/obc2 + gaff-2.11 (reuses the MD parameterization).
- Nonbonded: CutoffNonPeriodic, cutoff 1.0 nm (same NonbondedForce class as the MD).
- Per-frame component split via OpenMM force-group dispatch.
- Equilibration discard: auto 20 ps for explicit MD (waters relax around the unrestrained solute post-restraint-release).
- Configurational entropy: computed but EXCLUDED from ranking. An Interaction-Entropy −TΔS term (Gaussian second-order and full exponential-average backends) is available, but its uncertainty swamps the enthalpic signal: across 30 MM-GBSA runs the −TΔS SEM is ~19× the ΔH SEM (median 8.4 vs 0.44 kcal/mol), with −TΔS magnitudes of 23–172 kcal/mol and the full exponential-average estimate one-frame-dominated in 29/30 runs. The decisive point is reproducibility: for 29 replicate/pose MD runs of a SINGLE compound (paclitaxel) in one correctly-prepared taxane pocket the −TΔS term is not self-consistent (23–172 kcal/mol run-to-run), and adding it collapses the clean ΔH separation between paclitaxel and the weaker-binding juliprosopine control. Ranking therefore keys on entropy-free ΔH, and each free-energy record carries an `entropy.ranking.unusable_for_ranking` flag. The reported number is enthalpic + solvation only — the standard MM-GBSA quantity. (Normal-mode / quasi-harmonic backends remain out of scope.)
- Single-trajectory subtraction: ΔG_bind = ⟨E_complex⟩ − ⟨E_receptor⟩ − ⟨E_ligand⟩.

Per-run artifacts: `jobs/<md_id>/free_energy/summary.json` (the full self-describing FE block written by `mmgbsa_runner.compute_md_fe`).