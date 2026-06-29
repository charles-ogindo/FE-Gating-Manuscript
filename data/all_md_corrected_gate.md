# FE gate (final 2026-06-16) — every openmm_full MD run

C1 collapses to a single eq-window begin-vs-end check — slope (a) and visual trend (c) are gone; both were measuring the same drift less directly than the begin/end means. The full gate is now 4 criteria. C2 (self-stable) remains removed (redundant). `summary.verdict` ignored; implicit-solvent runs flagged.

## Criteria

**1. CONVERGED** — **|eq_end_mean − eq_begin_mean| < 0.4 Å**. eq_begin = mean pocket-aligned RMSD over the first **10 %** of the eq-window frames; eq_end = mean over the last **10 %**. For a 60-frame eq window that's the first and last 6 frames. Asks: once the ligand entered its settled window, did RMSD stop changing?

**2. SELF-STABLE — REMOVED** (redundant once C1 + C3 + C4 hold).

**3. IN-POCKET** — in-pocket fraction in eq window ≥ **0.95** AND COM displacement at the final frame ≤ **3.0 Å**.
**4. INTERACTIONS** — **top-5 residue contact persistence** (≥ **65 %**). The top-5 are the 5 protein residues with the highest per-frame contact fraction across the whole trajectory (any-heavy-atom ≤ **4.0 Å**). Per eq-window frame: count how many of those 5 are in contact. C4 metric = mean of (count ÷ 5). Captures “key interactions held” while tolerating the all-pair reshuffle.
**5. EQ-WINDOW FRAMES** — ≥ **50** frames in the last 30 % of the trajectory.

## Result

**5 of 11 openmm_full runs QUALIFY for FE computation.**

**Solvent-strategy note (2026-06-16):** explicit TIP3P is the only supported FE-feeder. Implicit-solvent OBC2 runs are kept in this table as historical evidence — even if they qualify on the metric criteria they are NOT FE candidates; the **Notes** column flags them as `implicit lineage — declined strategy`.

## Summary table

| Job ID (8) | Compound | Sol | Δ(begin−end) Å | C3 in-pocket | C4 top-5 % | C5 eq-frames | **QUALIFIES** | Notes |
|---|---|---|---|---|---|---|---|---|
| `02f30602` | taxol | implicit | ✅ (0.28 Å) | ✅ (frac 1.00, COM 1.39 Å) | ✅ (86.33%) | ✅ (60/50) | ✅ **YES** | **implicit lineage** — declined strategy (explicit TIP3P is the only supported FE-feeder; kept for evidence, not an FE candidate) · verdict.json said `drifting` (ignored) |
| `5bc61f59` | taxol | explicit | ✅ (0.07 Å) | ✅ (frac 1.00, COM 1.29 Å) | ✅ (99.33%) | ✅ (60/50) | ✅ **YES** | verdict.json said `drifting` (ignored) |
| `34840aa1` | taxol | explicit | ✅ (0.40 Å) | ✅ (frac 1.00, COM 0.85 Å) | ✅ (99.67%) | ✅ (60/50) | ✅ **YES** | verdict.json said `drifting` (ignored) |
| `80e53d8a` | taxol | explicit | ✅ (0.17 Å) | ✅ (frac 1.00, COM 1.01 Å) | ✅ (100.00%) | ✅ (60/50) | ✅ **YES** | verdict.json said `drifting` (ignored) |
| `a0b04941` | Juliprosopine | explicit | ✅ (0.29 Å) | ✅ (frac 1.00, COM 2.13 Å) | ✅ (81.33%) | ✅ (60/50) | ✅ **YES** | verdict.json said `unstable` (ignored) |
| `1f01da83` | taxol | — | ✅ (0.17 Å) | ✅ (frac 1.00, COM 0.94 Å) | ✅ (98.00%) | ❌ (30/50) | ❌ no | verdict.json said `stable` (ignored) |
| `5817afff` | taxol | — | ❌ (0.79 Å) | ❌ (frac 1.00, COM 3.32 Å) | ✅ (97.33%) | ✅ (60/50) | ❌ no | verdict.json said `drifting` (ignored) |
| `6cfae070` | taxol | explicit | ✅ (0.14 Å) | ✅ (frac 1.00, COM 0.44 Å) | ✅ (100.00%) | ❌ (3/50) | ❌ no | verdict.json said `stable` (ignored) |
| `dba42e7d` | taxol | — | ✅ (0.13 Å) | ✅ (frac 1.00, COM 0.41 Å) | ✅ (100.00%) | ❌ (3/50) | ❌ no | verdict.json said `stable` (ignored) |
| `f9c29ee9` | taxol | — | ✅ (0.02 Å) | ✅ (frac 1.00, COM 1.12 Å) | ✅ (100.00%) | ❌ (3/50) | ❌ no | verdict.json said `stable` (ignored) |
| `00dda37f` | Primaquine | explicit | ❌ (0.85 Å) | ✅ (frac 1.00, COM 1.66 Å) | ✅ (94.67%) | ✅ (60/50) | ❌ no | verdict.json said `drifting` (ignored) |
| `86e5195b` | taxol (wrong-pocket) | — | — | — | — | — | — | skipped (no summary.json) |
| `d34b991f` | taxol (wrong-pocket) | — | — | — | — | — | — | skipped (no summary.json) |

## Detailed per-run blocks

### `02f30602-8624-4c0a-863b-9c130d357364` — taxol (pose 0, solvent=implicit, 1000.00 ps, 201 frames)

```
Job: 02f30602-8624-4c0a-863b-9c130d357364 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 3.59 vs 3.87 Å, Δ(begin−end) = 0.28 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 1.39 Å → PASS
Interactions: top-5 residue persistence 86.3% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): PRO B:271 (97%); ARG B:358 (93%); ARG B:281 (85%); HIS B:226 (84%); ARG B:275 (81%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: YES** | Notes: **implicit lineage** — declined strategy. verdict.json said `drifting` (ignored)
```

### `5bc61f59-834f-4e71-a492-d32ddfdc7326` — taxol (pose 0, solvent=explicit, 1000.00 ps, 201 frames)

```
Job: 5bc61f59-834f-4e71-a492-d32ddfdc7326 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 2.27 vs 2.34 Å, Δ(begin−end) = 0.07 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 1.29 Å → PASS
Interactions: top-5 residue persistence 99.3% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): THR B:273 (100%); PRO B:271 (100%); VAL B:22 (100%); LEU B:360 (99%); HIS B:226 (99%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: YES** | Notes: verdict.json said `drifting` (ignored)
```

### `34840aa1-cbf5-4c6b-a665-ea4f52110f5d` — taxol (pose 1, solvent=explicit, 1000.00 ps, 201 frames)

```
Job: 34840aa1-cbf5-4c6b-a665-ea4f52110f5d | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 2.05 vs 1.65 Å, Δ(begin−end) = 0.40 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 0.85 Å → PASS
Interactions: top-5 residue persistence 99.7% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): ASP B:25 (100%); ARG B:358 (100%); HIS B:226 (100%); ARG B:275 (100%); LEU B:272 (100%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: YES** | Notes: verdict.json said `drifting` (ignored)
```

### `80e53d8a-1926-4525-b2e1-55cb1e30eedd` — taxol (pose 2, solvent=explicit, 1000.00 ps, 201 frames)

```
Job: 80e53d8a-1926-4525-b2e1-55cb1e30eedd | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 1.85 vs 2.02 Å, Δ(begin−end) = 0.17 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 1.01 Å → PASS
Interactions: top-5 residue persistence 100.0% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): HIS B:226 (100%); PRO B:357 (100%); ALA B:230 (100%); THR B:273 (100%); VAL B:22 (99%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: YES** | Notes: verdict.json said `drifting` (ignored)
```

### `a0b04941-e4b0-40ef-9459-67becac4a61c` — Juliprosopine (pose 0, solvent=explicit, 1000.00 ps, 201 frames)

```
Job: a0b04941-e4b0-40ef-9459-67becac4a61c | Compound: Juliprosopine | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 5.01 vs 4.73 Å, Δ(begin−end) = 0.29 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 2.13 Å → PASS
Interactions: top-5 residue persistence 81.3% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): ARG B:358 (94%); PRO B:357 (92%); LEU B:360 (80%); ALA B:230 (61%); PRO B:271 (60%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: YES** | Notes: verdict.json said `unstable` (ignored)
```

### `1f01da83-ac18-4103-8631-2053873967c0` — taxol (pose 0, solvent=None, 500.00 ps, 101 frames)

```
Job: 1f01da83-ac18-4103-8631-2053873967c0 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (3 frames each): 1.46 vs 1.64 Å, Δ(begin−end) = 0.17 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 0.94 Å → PASS
Interactions: top-5 residue persistence 98.0% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): THR B:273 (100%); LEU B:272 (100%); VAL B:22 (100%); ARG B:275 (98%); PRO B:271 (98%)
Equilibrated frames: 30 (threshold 50) → FAIL

**QUALIFIES: NO** | Notes: verdict.json said `stable` (ignored)
```

### `5817afff-8e6b-4130-9997-46c9965379b2` — taxol (pose 0, solvent=None, 1000.00 ps, 201 frames)

```
Job: 5817afff-8e6b-4130-9997-46c9965379b2 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 3.28 vs 4.08 Å, Δ(begin−end) = 0.79 Å (threshold < 0.4 Å) → FAIL

In-pocket: fraction 1.00, COM 3.32 Å → FAIL
Interactions: top-5 residue persistence 97.3% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): ARG B:275 (100%); THR B:273 (100%); LEU B:272 (98%); PRO B:271 (98%); HIS B:226 (98%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: NO** | Notes: verdict.json said `drifting` (ignored)
```

### `6cfae070-eddf-4cbb-923d-679cfeb305e3` — taxol (pose 0, solvent=explicit, 10.00 ps, 11 frames)

```
Job: 6cfae070-eddf-4cbb-923d-679cfeb305e3 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (1 frames each): 0.82 vs 0.96 Å, Δ(begin−end) = 0.14 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 0.44 Å → PASS
Interactions: top-5 residue persistence 100.0% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): LEU B:272 (100%); GLU B:26 (100%); VAL B:22 (100%); GLU B:21 (100%); LEU B:360 (100%)
Equilibrated frames: 3 (threshold 50) → FAIL

**QUALIFIES: NO** | Notes: verdict.json said `stable` (ignored)
```

### `dba42e7d-71a0-4188-b00a-515e1c58388a` — taxol (pose 0, solvent=None, 10.00 ps, 11 frames)

```
Job: dba42e7d-71a0-4188-b00a-515e1c58388a | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (1 frames each): 0.92 vs 0.79 Å, Δ(begin−end) = 0.13 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 0.41 Å → PASS
Interactions: top-5 residue persistence 100.0% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): LEU B:360 (100%); PRO B:271 (100%); LEU B:272 (100%); THR B:273 (100%); ARG B:275 (100%)
Equilibrated frames: 3 (threshold 50) → FAIL

**QUALIFIES: NO** | Notes: verdict.json said `stable` (ignored)
```

### `f9c29ee9-7596-496d-b907-0969bb203a88` — taxol (pose 0, solvent=None, 10.00 ps, 11 frames)

```
Job: f9c29ee9-7596-496d-b907-0969bb203a88 | Compound: taxol | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (1 frames each): 1.90 vs 1.88 Å, Δ(begin−end) = 0.02 Å (threshold < 0.4 Å) → PASS

In-pocket: fraction 1.00, COM 1.12 Å → PASS
Interactions: top-5 residue persistence 100.0% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): PRO B:271 (100%); ARG B:281 (100%); ARG B:358 (91%); HIS B:226 (91%); THR B:273 (91%)
Equilibrated frames: 3 (threshold 50) → FAIL

**QUALIFIES: NO** | Notes: verdict.json said `stable` (ignored)
```

### `00dda37f-e218-49ce-8ac0-a69bfff6851d` — Primaquine (pose 0, solvent=explicit, 1000.00 ps, 201 frames)

```
Job: 00dda37f-e218-49ce-8ac0-a69bfff6851d | Compound: Primaquine | Engine: openmm_full
──────────────────────────────────────────────────────────────────────────
C1 Converged: eq-window begin vs end (6 frames each): 3.86 vs 3.01 Å, Δ(begin−end) = 0.85 Å (threshold < 0.4 Å) → FAIL

In-pocket: fraction 1.00, COM 1.66 Å → PASS
Interactions: top-5 residue persistence 94.7% (threshold 65%) → PASS
  Top-5 residues (whole-trajectory contact frac): ARG B:317 (100%); GLU B:26 (98%); PRO B:357 (96%); PHE B:269 (96%); VAL B:22 (95%)
Equilibrated frames: 60 (threshold 50) → PASS

**QUALIFIES: NO** | Notes: verdict.json said `drifting` (ignored)
```

### `86e5195b` — taxol (wrong-pocket)  _(skipped: no summary.json)_

### `d34b991f` — taxol (wrong-pocket)  _(skipped: no summary.json)_

## QUALIFYING runs

- `02f30602-8624-4c0a-863b-9c130d357364` — **taxol** pose 0, solvent=implicit, 1000.00 ps, 201 frames
- `5bc61f59-834f-4e71-a492-d32ddfdc7326` — **taxol** pose 0, solvent=explicit, 1000.00 ps, 201 frames
- `34840aa1-cbf5-4c6b-a665-ea4f52110f5d` — **taxol** pose 1, solvent=explicit, 1000.00 ps, 201 frames
- `80e53d8a-1926-4525-b2e1-55cb1e30eedd` — **taxol** pose 2, solvent=explicit, 1000.00 ps, 201 frames
- `a0b04941-e4b0-40ef-9459-67becac4a61c` — **Juliprosopine** pose 0, solvent=explicit, 1000.00 ps, 201 frames