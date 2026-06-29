# Changelog

All notable changes to this reviewer / reproduction repository are recorded here.

## 2026-06-29 — Repository consolidated and corrected

This repository was reorganized into a single canonical reviewer package and its
history was rewritten. If you cloned an earlier version, please re-clone or reset
to the current `main`.

Changes in this revision:

- Corrected the gate description to match the manuscript: two structural
  conditions (core-contact retention; pose convergence) plus a separate
  sampling-adequacy check (N_eff ≥ 10). The earlier "0.75 Å, mean+3SD"
  description was superseded.
- Added the sanitized method source, the Table 4 reproduction artifacts,
  figures, structures, and notes into one repository.
- Removed local filesystem paths from `structures/ligand_pose0.pdb` and the
  `data/` gate tables (coordinate and result data unchanged).
