# FE gate applied to every openmm_full MD run

Loops `backend.app.free_energy.gating.validate(md_job_id)` over every MD run flagged engine=`openmm_full` in the catalogue (12 runs). Sourced from db.MDRun ∪ filesystem summary.json (engine).

## Gate thresholds (verbatim from `backend.app.free_energy.__init__`)

- `DEFAULT_MIN_TOTAL_FRAMES`   = **200**    (hard gate — ≥ 200 frames before any windowing)
- `DEFAULT_MIN_WINDOW_FRAMES`  = **50**     (hard gate — ≥ 50 frames in the sampled window)
- `LIGAND_MAX_RMSD_WARN_A`     = **5.0** Å   (soft gate — warns above this)
- `LIGAND_MAX_RMSD_BLOCK_A`    = **10.0** Å   (soft gate — *blocks* if ligand max RMSD ≥ this — "unbound")
- `BACKBONE_FINAL_RMSD_WARN_A` = **2.5** Å  (soft gate — backbone-not-equilibrated warning)
- `WINDOW_DURATION_WARN_PS`    = **1000.0** ps   (soft gate — short post-eq window warning)

The gate's hard gates (must pass all):
  A `MD_SUMMARY_MISSING` — md/summary.json present + parseable
  B `MD_NOT_COMPLETE`    — summary.status == 'completed'
  C `MD_NOT_STABLE`      — summary.verdict == 'stable'
  D `MD_NOT_OPENMM_FULL` — summary.engine.kind == 'openmm_full'
  E `MD_FRAMES_MISSING` / `MD_FRAMES_EMPTY` / `MD_FRAME_UNREADABLE`
  F `INSUFFICIENT_FRAMES` — n_total ≥ 200 AND n_sampled ≥ 50

Soft gates emit warnings but do NOT block: `LIGAND_DRIFT` (2-5 Å warn), `BACKBONE_NOT_EQUILIBRATED` (BB Cα RMSD final ≥ 2.5 Å), `OUTLIER_FRAMES`, `SHORT_WINDOW_TIME` (< 1 ns sampled), `STRIDE_MANUAL`, `EXPERIMENTAL_ENGINE` (standing flag — every passing run gets it; this is the "experimental / approximate" chip the dashboard shows on the FE estimate).

## Result

Out of 12 openmm_full runs, **0 pass** the FE gate.

| job_id (8) | compound | pose | dur(ps) | n_frm | ligand pose RMSD (fin/max) | BB final | gate | hard-gate reasons | soft-gate warnings |
|---|---|---|---|---|---|---|---|---|---|
| `1f01da83` | taxol | 0 | 500.00 | 101 | 1.51 / 2.51 | 3.14 | ❌ decline | INSUFFICIENT_FRAMES | — |
| `5817afff` | taxol | 0 | 1000.00 | 201 | 3.84 / 4.76 | 4.08 | ❌ decline | MD_NOT_STABLE | — |
| `5bc61f59` | taxol | 0 | 1000.00 | 201 | 2.63 / 3.23 | 2.18 | ❌ decline | MD_NOT_STABLE | — |
| `6cfae070` | taxol | 0 | 10.00 | 11 | 1.69 / 1.69 | 1.37 | ❌ decline | INSUFFICIENT_FRAMES | — |
| `dba42e7d` | taxol | 0 | 10.00 | 11 | 0.82 / 1.30 | 1.34 | ❌ decline | INSUFFICIENT_FRAMES | — |
| `f9c29ee9` | taxol | 0 | 10.00 | 11 | 1.75 / 1.94 | 1.26 | ❌ decline | INSUFFICIENT_FRAMES | — |
| `34840aa1` | taxol | 1 | 1000.00 | 201 | 2.53 / 3.23 | 2.92 | ❌ decline | MD_NOT_STABLE | — |
| `80e53d8a` | taxol | 2 | 1000.00 | 201 | 2.31 / 3.60 | 2.66 | ❌ decline | MD_NOT_STABLE | — |
| `a0b04941` | Juliprosopine | 0 | 1000.00 | 201 | 5.11 / 5.70 | 3.04 | ❌ decline | MD_NOT_STABLE | — |
| `00dda37f` | Primaquine | 0 | 1000.00 | 201 | 2.23 / 3.99 | 2.68 | ❌ decline | MD_NOT_STABLE | — |
| `86e5195b` | taxol (wrong-pocket) | 0 | — | — | — / — | — | ❌ decline | MD_SUMMARY_MISSING | — |
| `d34b991f` | taxol (wrong-pocket) | 0 | — | — | — / — | — | ❌ decline | MD_SUMMARY_MISSING | — |

## Per-run gate-message detail (the why)

- **`1f01da83`** — taxol pose 0 (500.00 ps, 101 frames, verdict=`stable`): ❌ **DECLINE**  hard=[`INSUFFICIENT_FRAMES`], soft=[—]
  - hard-gate msgs: INSUFFICIENT_FRAMES: MD produced 101 frames; need at least 200 total. Re-run MD with longer production (or smaller snapshot interval).
- **`5817afff`** — taxol pose 0 (1000.00 ps, 201 frames, verdict=`drifting`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'drifting'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`5bc61f59`** — taxol pose 0 (1000.00 ps, 201 frames, verdict=`drifting`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'drifting'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`6cfae070`** — taxol pose 0 (10.00 ps, 11 frames, verdict=`stable`): ❌ **DECLINE**  hard=[`INSUFFICIENT_FRAMES`], soft=[—]
  - hard-gate msgs: INSUFFICIENT_FRAMES: MD produced 11 frames; need at least 200 total. Re-run MD with longer production (or smaller snapshot interval).
- **`dba42e7d`** — taxol pose 0 (10.00 ps, 11 frames, verdict=`stable`): ❌ **DECLINE**  hard=[`INSUFFICIENT_FRAMES`], soft=[—]
  - hard-gate msgs: INSUFFICIENT_FRAMES: MD produced 11 frames; need at least 200 total. Re-run MD with longer production (or smaller snapshot interval).
- **`f9c29ee9`** — taxol pose 0 (10.00 ps, 11 frames, verdict=`stable`): ❌ **DECLINE**  hard=[`INSUFFICIENT_FRAMES`], soft=[—]
  - hard-gate msgs: INSUFFICIENT_FRAMES: MD produced 11 frames; need at least 200 total. Re-run MD with longer production (or smaller snapshot interval).
- **`34840aa1`** — taxol pose 1 (1000.00 ps, 201 frames, verdict=`drifting`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'drifting'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`80e53d8a`** — taxol pose 2 (1000.00 ps, 201 frames, verdict=`drifting`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'drifting'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`a0b04941`** — Juliprosopine pose 0 (1000.00 ps, 201 frames, verdict=`unstable`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'unstable'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`00dda37f`** — Primaquine pose 0 (1000.00 ps, 201 frames, verdict=`drifting`): ❌ **DECLINE**  hard=[`MD_NOT_STABLE`], soft=[—]
  - hard-gate msgs: MD_NOT_STABLE: MD stability verdict is 'drifting'; need 'stable'. Inspect MD results before requesting a free-energy run.
- **`86e5195b`** — taxol (wrong-pocket) pose 0 (— ps, — frames, verdict=`—`): ❌ **DECLINE**  hard=[`MD_SUMMARY_MISSING`], soft=[—]
  - hard-gate msgs: MD_SUMMARY_MISSING: No MD summary at jobs/86e5195b-5f5c-4ecc-a221-199731ddb4bd/md/summary.json
- **`d34b991f`** — taxol (wrong-pocket) pose 0 (— ps, — frames, verdict=`—`): ❌ **DECLINE**  hard=[`MD_SUMMARY_MISSING`], soft=[—]
  - hard-gate msgs: MD_SUMMARY_MISSING: No MD summary at jobs/d34b991f-f4a3-4105-b58c-b2c3c8e3cf49/md/summary.json