# Beta-tubulin / taxol example inputs

Minimal inputs to reproduce the taxane-site docking + explicit-solvent MD +
MM-GBSA replicate study. **No receptor coordinates are committed here** — use
the public crystal structure instead.

## Receptor

Use **PDB [1JFF](https://www.rcsb.org/structure/1JFF)** — the 3.5 Å Zn-induced
αβ-tubulin sheet with bound **taxol (ligand TA1)**, GTP, GDP, and Mg²⁺. Download
it yourself:

```bash
wget https://files.rcsb.org/download/1JFF.pdb
```

The docked system models chains A (α-tubulin) and B (β-tubulin) with the GTP /
GDP / Mg²⁺ / Zn²⁺ cofactors; taxol binds the β-tubulin taxane pocket.

## Files

| File | What it is |
|---|---|
| `ligands.smi` | The three docked ligands (SMILES + name). `taxol` is the reference / `lig_000000` simulated in the replicate sweep. |
| `box.json` | Docking search box for the taxane site, **center + size in Å**, derived from the reference taxol (TA1) ligand identity (the corrected ligand-identity box logic, not the old "smallest chain" heuristic). |
| `md_settings.json` | Explicit-solvent MD settings (octahedron box, 1.0 nm TIP3P padding, 0.15 M, and the A439–451 α-tubulin E-hook truncation). Maps to `backend.app.md.job.MdRunSettings`. |

## Notes

- The box was anchored on TA1 in the holo reference complex; against an apo
  model you would re-derive it from the taxane pocket.
- `truncate_chain_ranges` removes the disordered α-tubulin C-terminal E-hook
  (residues 439–451, `VEGEGEEEGEEY`) and OXT-caps A438. This is required to
  match the reference ~191k-atom solvated system; leaving it in inflates the
  octahedron toward ~520k atoms.
