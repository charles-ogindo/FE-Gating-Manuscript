# Publication notes — `docking_app2` (process / milestone register)

Companion to `publication_record.md`. Where `publication_record.md` is the
**claims** register (paper-shaped, defensible scientific statements with
sources and caveats), this file is the **process** register: a dated log
of milestones with the methods, parameters, versions, decisions, and
figures/data each one touched.

Both files are updated on the same cadence; commits land with a footer
line `(publication-notes: YYYY-MM-DD)` pointing back to the entry that
records them.

**Scope of "MD" in this register (2026-06-05 curation pass).** Reported
MD = full-OpenMM trajectories only; surrogate runs and abandoned-setup /
mis-specified-pocket docking directories archived under
`jobs/archive/` (recoverable, excluded from analysis scans + free-energy
gating); non-binder controls retained as true negatives. See
`jobs/archive/MANIFEST.md` for the per-job classification.

Entry block format (use verbatim):

```markdown
## YYYY-MM-DD — <Phase or item label> — <one-line summary>
**Commit(s):** `<sha>`
**Scope:** backend | frontend | both
**Supports claim(s):** <date(s) of related publication_record.md entries, or —>
**What changed:** <1–3 sentences>
**Methods / params / versions touched:** <thresholds, force fields, tool
  versions, runtime knobs — exactly what a methods section would state>
**Decisions + rationale:** <locked-in choices + why; "we chose X because Y">
**Figures / data produced:** <viewer URL w/ full params; artifact paths; job IDs>
**Caveat / known limitation:** <if any>
```

---

## 2026-06-03 — Phase 1A — analysis/interactions submodule (detection + thresholds)
**Commit(s):** `e2762f6`
**Scope:** backend
**Supports claim(s):** —
**What changed:** Added the analysis/interactions/ submodule (12 files) detecting
9 protein–ligand interaction types per pose and emitting them to a per-pose
contacts artifact.
**Methods / params / versions touched:** 9 types (H-bond, weak C–H···O,
hydrophobic, salt bridge, π-stacking, π-cation, halogen, water bridge, metal).
Per-type geometric cutoffs live in analysis/interactions/thresholds.py with
literature citations (IUPAC, PLIP, Auffinger, Desiraju) and are serialized
verbatim into every contacts.json as a `thresholds_applied` block — this is the
methods-section source of truth for the detection criteria.
**Decisions + rationale:** Thresholds are data, not code constants — one cited
module, emitted into each artifact, so every figure/table is self-describing and
reproducible without reading source.
**Figures / data produced:** contacts artifact schema (interactions[],
thresholds_applied), consumed by the viewer and exports.
**Caveat / known limitation:** π-stacking detection incomplete when
ligand_smiles is null (no aromatic-ring perception).

## 2026-06-03 — Phase 1B — contacts.json persistence + backfill
**Commit(s):** `e3bac1f`
**Scope:** backend
**Supports claim(s):** —
**What changed:** Wired the interactions submodule into dock-time export (each
pose writes contacts_NN.json) and backfilled existing jobs (30/30 present).
**Methods / params / versions touched:** No new chemistry; uses Phase 1A
detection + thresholds. Artifact path docking/viewer/<ligand>/contacts_NN.json
(rank-padded).
**Decisions + rationale:** Backfill over compute-on-demand so historical jobs are
immediately publication-ready and the artifact set is complete and auditable.
**Figures / data produced:** 30/30 contacts_NN.json across existing jobs.
**Caveat / known limitation:** Scientific-stack versions that generated these
artifacts are not pinned in git — record the conda env if these specific
artifacts back a figure (see environment-snapshot.md).

## 2026-06-04 — Item 23 (v1) — analysis/consensus_pharmacophore/ backend engine
**Commit(s):** `b12a8d1`
**Scope:** backend
**Supports claim(s):** —
**What changed:** New backend submodule `analysis/consensus_pharmacophore/` that
extracts a recurrent-feature consensus pharmacophore across many `(ligand, pose)`
inputs for a single receptor. Generic in inputs; thin `assemble_from_job(job_id)`
adapter packages the canonical 3-ligand fixture. Emits
`analysis/consensus_pharmacophore/<run_id>/consensus_pharmacophore.json` with the
self-describing-artifact + criteria-block discipline established by
`analysis/interactions/`.
**Methods / params / versions touched:** Four locked decisions, recorded verbatim
in every artifact's `criteria` block:
- **Feature primitive = contacts.** Reads `contacts_NN.json` interactions and
  projects each via a single `feature_map.map_interaction_type()` rule onto one
  of 9 pharmacophore types. H-bond direction (`protein_donor` vs `ligand_donor`)
  maps to the LIGAND-anchored role (acceptor / donor); weak C–H···O subtype →
  donor; missing direction → acceptor with `direction_missing` flag. The
  schema-reserved `aromatic` type is NOT produced by v1 (deferred Meeko
  cross-check pass).
- **Clustering = (feature_type, chain, residue_number) + receptor-atom-set
  overlap, union-find.** Transitively unioned within a bucket so "same sub-region
  of the same residue" spans interactions that overlap pairwise but not all-pairs.
- **Consensus = N-of-M majority, N = ceil(M / 2) distinct-ligand support.**
  `n_override` accepted; default majority puts the fixture at 2-of-3.
- **Pose-validation gate = OFF (v1).** Configurable flag; default off. ON
  requires the Q6 / Q7 / F12 work in the saved roadmap to land first.

Position = centroid of contributing ligand-side positions. Tolerance = max
distance of any contributing position from the centroid (Å). Activity-weighting
hook (`weight_fn(ligand_id) -> float`) wired through but unused in v1 (uniform).
Single choke point `extract.normalize_feature()` for `contacts.json` field
access — schema drift is a one-function fix.
**Decisions + rationale:** Contacts-only as the feature primitive (not Meeko or
both) because the contacts artifact is already cited, validated, and per-pose
self-describing — the Meeko ligand-side cross-check is an additive future pass,
not a v1 dependency. Distinct-ligand support count over raw-interaction count so
two interactions from the same ligand at the same residue don't inflate
consensus. Receptor-atom-set overlap as the clustering rule (rather than spatial
clustering on positions) because it's the cheapest correct rule that handles the
ND1/NE2 distinction within HIS without arbitrary distance cutoffs. Pose-gate
default OFF surfaced explicitly in every artifact's `criteria` block so a future
ON-by-default switch is a traceable policy decision.
**Figures / data produced:** Fixture smoke on job
`4a37bb0c-7655-411e-a238-1ff5b7bab910` (paclitaxel + juliprosopine + primaquine
vs β-tubulin) produced 11 consensus features at
`jobs/4a37bb0c-…/analysis/consensus_pharmacophore/smoke-2026-06-04/consensus_pharmacophore.json`
from 93 raw features → 43 candidate clusters → 32 dropped (singleton) → 11 kept.
Per-type: 1 hbond_acceptor (GLY370 backbone N, 2-of-3), 9 hydrophobic across the
canonical taxane site (VAL23, GLU27, HIS229, PHE272 ×2 distinct ring faces,
LEU275, PRO360, ARG369), 1 pi_cation (HIS229 imidazole, 2-of-3). pytest 55/55
PASS (32 existing + 23 new tests) in 40.26s.
**Caveat / known limitation:** **MECHANICS demonstration only — NOT a scientific
pharmacophore claim.** The 3 fixture ligands are not a curated actives panel: one
known binder + two designed-as-non-binder controls (juliprosopine, primaquine —
per `publication_record.md`'s 2026-06-02 multi-ligand selectivity entry). All
three Vina-dock into the same search box by construction, so spatial neighborhood
is forced; "consensus" here means "all three happened to put a heavy atom near
the same receptor residue," not "all three engage the same chemistry." Scientific
validity requires (a) a curated actives panel against a single receptor, (b)
measured activities for the weighting hook, (c) MD-validated poses with the gate
ON — which in turn requires Q6 (backbone-Cα-RMSD investigation), Q7 (free-energy
gating on artifactual RMSD), and F12 (MD-1c root-cause fix at source) from the
saved roadmap. Each artifact emits a `caveats` block surfacing this gap.
No frontend render this pass (Item 23 frontend deferred until the Phase 1C/1D
browser verifications V1–V4 land); the Item 23 frontend will reuse the
`useMolstarContacts` rendering pattern + the `contactsExport.normalizeContact`
choke-point philosophy with a new `consensusExport.js`.

**Update (2026-06-04 post-Item-25): hydrophobic dominance fixed.** After the
Item 25 classifier landed in commit `3af9025` and the 30-artifact backfill ran,
re-running this fixture's consensus on the same 3-ligand input shrinks the
hydrophobic-dominance shape dramatically. Run id: `post-item-25-2026-06-04`.

| | Pre (`smoke-2026-06-04`, this entry's original numbers) | Post (`post-item-25-2026-06-04`) |
|---|---:|---:|
| raw features | 93 | 29 |
| candidate clusters | 43 | 21 |
| consensus features kept | 11 | 4 |
| hbond_acceptor | 1 | 1 |
| hydrophobic | 9 | **2** |
| pi_cation | 1 | 1 |

The 2 surviving hydrophobic features are at canonical taxane-site residues
(LEU275 CD1, PHE272 CZ). The 7 dropped were all at residues outside the
sidechain-7 hydrophobic set (importance=0) OR were non-primary contacts within a
qualifying residue (one primary per residue per the Item 25 rule). Hydrophobic
dominance dropped from 82% (9/11) of the consensus to 50% (2/4); specific
interactions (H-bond + π-cation) went from 18% (2/11) to 50% (2/4). Signal-to-
noise improvement matches the design intent. Mechanics-only caveat above STILL
applies — the four kept features all sit at canonical taxane-site residues,
which is the expected pattern when one of three input ligands is a known taxol
binder; it remains NOT a scientific pharmacophore claim until a curated actives
panel + MD-validated poses are wired.

## 2026-06-04 — Item 25 backend — hydrophobic importance classifier + schema 1.1.0 + backfill
**Commit(s):** `3af9025`
**Scope:** backend
**Supports claim(s):** —
**What changed:** New `analysis/interactions/hydrophobic_classifier.py` (single
function `apply_hydrophobic_classification`) annotates every hydrophobic
interaction in `contacts_NN.json` with `importance` (per-residue qualifying-
contact count), `is_primary` (one per qualifying residue, smallest distance), and
`buried_surface_area_a2: null` (v1 placeholder). Schema bump
`contacts_NN.json: 1.0.0 → 1.1.0` with a new top-level `criteria` key carrying
the `hydrophobic_classification` block. New backfill script
`scripts/backfill_hydrophobic_classification.py` mirrors the Phase 1B pattern;
30/30 existing artifacts updated; idempotency confirmed on re-run.
`analysis/consensus_pharmacophore/extract.py` now filters non-primary
hydrophobics via a new `skipped_hydrophobic_non_primary` stat — fixes the
Item 23 v1 9-of-11 hydrophobic dominance.
**Methods / params / versions touched:** Three locked decisions, recorded verbatim
in every artifact's `criteria.hydrophobic_classification` block:
- **Importance metric** = per-residue count of qualifying ligand-atom ↔
  hydrophobic-sidechain-carbon contacts where the receptor residue is in
  `{LEU, ILE, VAL, PHE, MET, TRP, ALA}`. Non-qualifying residues get `importance=0`.
- **`is_primary` rule** = the single smallest-distance contact at each qualifying
  residue; non-qualifying residues have no primary. Tie-break by stable first-seen
  order in the artifact.
- **`buried_surface_area_a2`** = `null` in v1 (forward-compatible field; BSA
  refinement deferred to v2).

Backbone-carbon exclusion is INHERITED via the upstream
`atom_typing.is_hydrophobic_carbon` (applied by the Item 14 detector when the
hydrophobic interactions are first emitted into `contacts_NN.json`). The
classifier does NOT re-check backbone exclusion — single-source-of-truth.

Single function shared between dock-time orchestrator emission and the backfill
script. `is_already_classified()` idempotency probe keys off
`criteria.hydrophobic_classification` presence, NOT `schema_version` — so a
future 1.2.0 bump unrelated to this classifier won't falsely report as
already-classified.

The new top-level `criteria` key is separate from the existing
`thresholds_applied` block — `thresholds_applied` = per-type GEOMETRIC cutoffs
(D-A distance, angle), `criteria` = downstream CLASSIFICATION rules. Clean
separation; consumers read whichever applies.
**Decisions + rationale:** Classifier as a SHARED capability (one source, three
consumers — Items 14, 23, 26) is the architectural invariant Item 25 exists to
establish; the 9-of-11 hydrophobic dominance in the Item 23 v1 smoke is the
demonstrated failure mode of re-deriving per consumer. Per-residue count over
per-pair count for `importance` — chemically meaningful unit. Sidechain-7 set
chosen for sidechain-carbon richness without polar/charged functional groups;
glycine intentionally excluded (no sidechain carbon at all). Smallest-distance
tie-break for `is_primary` because the closest contact is the most chemically
informative representative. BSA deferred because it adds an SASA dependency
(MSMS / FreeSASA / RDKit) without changing the qualitative signal in v1; the
field is reserved in the schema so adding it later is backward-compatible. The
new top-level `criteria` block (not extending `thresholds_applied`) preserves
the clean separation of detection cutoffs from downstream gating, matching the
`criteria` key naming already used by `analysis/consensus_pharmacophore/`.
**Figures / data produced:** 30 contacts_*.json artifacts updated under
`jobs/4a37bb0c-…` and `jobs/d96d719a-…`. Across all 30: `hydrophobic_total=687`,
`qualifying_residues=102`, `non_qualifying_residues=132`, `primary_count=102`
(`qualifying_contacts=342`, `non_qualifying_contacts=345`). Re-run consensus for
the canonical 3-ligand job at
`jobs/4a37bb0c-…/analysis/consensus_pharmacophore/post-item-25-2026-06-04/consensus_pharmacophore.json`
(see the Item 23 entry above for the pre/post comparison table). pytest 73/73
PASS (55 prior + 18 new tests) in 22.58s.
**Caveat / known limitation:** Frontend pieces of Item 25 (glyph encoding,
legend overlay, low-weight arc rendering for hydrophobic-primary) are
deferred until the Phase 1C/1D V1–V4 browser-verify lands — same gate that holds
the Item 23 frontend. Pre-1.1.0 artifacts (none remain after this backfill, but
defensive code is wired for them) pass non-primary hydrophobics through
unchanged because `is_primary is None` is treated as permissive. The
sidechain-7 set is the conservative choice; if the publication corpus has a
target receptor with consistently strong CYS / HIS hydrophobic contacts in the
literature, an expansion to a sidechain-N set is a one-line constant edit (also
naturally records into the artifact's `qualifying_residues` field for
reproducibility).

## 2026-06-04 — bugfix — endpoint xyz in contacts.json (schema 1.2.0)
**Commit(s):** backend `77d5a94`, frontend `a940d05`
**Scope:** both
**Supports claim(s):** —
**What changed:** Browser verification of the Phase 1C/1D contacts view (V1–V4
in the saved roadmap) surfaced two visible rendering bugs that share a single
root cause. Fixed at the backend layer per the spec's decision rule (artifact
lacked endpoint coords → emit them at detection time, not at render time).
Schema bump `1.1.0 → 1.2.0`; backfill the 30 existing artifacts via the
existing Phase 1B script with `--overwrite`.

**Bug 1 (fan-out)**: every Mol* cylinder emanated from a small set of shared
ligand-side points instead of running between the actual interacting atoms.

**Bug 2 (contactsLigand no-op)**: per-ligand URLs (`?contactsLigand=lig_000001`)
rendered indistinguishably from the all-ligands view.

**Root cause (shared)**: `contacts_NN.json` carried atom NAMES only — never
indices, never coordinates. The frontend `buildLigandIndexes()` in
`useMolstarContacts.js` built `Map<atomName, xyz>` from the loaded ligand pose
and silently dropped duplicates. Because obabel-converted ligand poses use
generic atom names — every C is named "C", every O is "O" — the Map only
stored ONE coord per element letter (whichever atom was parsed first). Every
interaction with `ligand_atom: "O"` collapsed to the same xyz; every "C" to
another shared xyz; every "N" to a third. This is the fan-out (Bug 1). Bug 2 is
a symptom of the same flaw: all three fixture ligands dock into the same
search box, so each ligand's "first C / first O" coords land in roughly the
same neighborhood; switching contactsLigand produced visually indistinguishable
cylinders even though the per-ligand contacts.json data was correctly
different. The backend HAD the precise atom indices at detection time (every
detector iterates with concrete pandas indices and computes distances from
resolved xyz) but discarded them before emitting.

**Methods / params / versions touched:** Schema bump `1.1.0 → 1.2.0`. New
fields under `atoms`:
- `receptor_xyz: List[List[float]]` — one entry per receptor atom name, OR a
  single-element centroid for centroid-based detectors (salt_bridge,
  pi_stacking, pi_cation/protein_cation).
- `ligand_xyz: List[List[float]]` — same shape on the ligand side.
- `water_xyz: List[List[float]]` — only present for water_bridge interactions.

New helpers `xyz_of_row(row) / xyz_of_position(vec)` in
`analysis/interactions/atom_typing.py` produce rounded 3-decimal lists. All 7
detectors updated (hbond strong + weak_ch_o, hydrophobic, salt_bridge,
pi_stacking, pi_cation × 2, halogen, water_bridge, metal) to emit the new
fields at their existing dict-build site.

Frontend: `buildCylinderShape()` reads `atoms.receptor_xyz` /
`atoms.ligand_xyz` and array-centroids via a new `_arrayCentroid()` helper
(no-op on length 1; ring-centroid on length N). The 73-line
`buildLigandIndexes()` and 19-line `selectLigandStructure()` functions are
deleted — both were dead code under the new design AND were the root cause of
the fan-out bug. The receptor atom index (`buildReceptorAtomIndex`) stays
because residue labels still anchor at Cα. Pre-1.2.0 artifacts are now an
error condition surfaced via a single `console.warn` pointing at the backfill
command — deliberate, because keeping the legacy fallback would re-introduce
the fan-out on any artifact that slipped past the backfill.

**Decisions + rationale:** Backend emission (not frontend repair) because the
backend has the resolved atom indices at detection time AND the artifact is
the durable contract — emitting coords makes every downstream consumer
(Mol* viewer, 2D diagram Item 26, consensus engine Item 23) read the same
authoritative endpoints. Centroid emission for centroid-based detectors
(rather than per-atom arrays of all ring atoms) because the detector ALREADY
computes the centroid as its distance reference; emitting it directly is the
cheapest correct representation. Killing the legacy atom-name fallback in the
frontend rather than keeping it as a "soft" path because the bug ONLY existed
because the soft path looked correct on cursory inspection but produced the
wrong geometry; a hard `console.warn` + missing-coord-skip is louder and
points the operator at the fix.

**Figures / data produced:** 30 contacts_*.json artifacts refreshed at
schema 1.2.0 under `jobs/4a37bb0c-…` and `jobs/d96d719a-…`. Artifact-level
verification: lig_000002's three PHE272 hydrophobic contacts carry three
DISTINCT `ligand_xyz` values (`(0.56,-14.94,10.07)`, `(5.84,-9.96,12.04)`,
`(1.82,-14.10,12.17)`); same atom name "C" across taxol/juliprosopine/primaquine
resolves to three distinct positions in the binding pocket. Backend pytest
73/73 PASS in 29.20s (one existing schema-version assertion updated to read
the `SCHEMA_VERSION` constant instead of a hardcoded literal so future bumps
don't re-break it). Frontend Vite production build clean (40.38s).

**Caveat / known limitation:** Browser verification pending on the user — per-
ligand URLs should now render visibly distinct cylinder geometry, and
individual cylinders should run between actual interacting atom positions
rather than fanning out from element-representative points. This pass is
geometry + filtering only; the Item 25 frontend pieces (glyphs, legend, low-
weight hydrophobic-primary arcs) still wait on the Phase 1C/1D V1–V4
browser-verify gate. Pre-1.2.0 contacts.json artifacts (none should exist on
this machine after the backfill, but any imported from elsewhere or generated
by an older backend in the future) will surface as missing cylinders + a
single console.warn — by design.
## 2026-06-04 — Q6b — ligand RMSD reference-frame fix + MD residue renumber (schema 1.1.0)
**Commit(s):** `e79bfea`
**Scope:** backend
**Supports claim(s):** — (analysis-correctness fix; verdict reclassification of existing MD jobs is descriptive, not a new scientific claim)
**What changed:** Two-part MD-analysis correction pass with no MD re-runs.
PART 1 (ligand RMSD frame): the pre-Q6b ligand RMSD was
`kabsch_rmsd(ref_lig, lig)` — a Kabsch SUPERPOSITION of the ligand on
itself, which removes the rigid-body translation/rotation that "did it
leave the pocket?" depends on. The 0.65 Å reported for taxol back-dock
was the ligand's INTERNAL conformational RMSD (sliding torsions, bond
rotations), not pocket retention. The new metric superposes each frame's
receptor backbone Cα onto the reference's, applies that SAME (R, t) to
the ligand heavy atoms, and computes plain RMSD vs the reference ligand
— this is POSE RMSD = pocket displacement. New
`utils/rmsd.kabsch_fit(ref, mob) -> (R, t, rmsd)` returns the transform
(the old scalar `kabsch_rmsd` and its callers are untouched).
`compute_rmsd_series` now returns 4-tuples
`(time_ps, bb_ca, lig_internal, lig_pose)` and `write_rmsd_csv` writes
the new `rmsd_ligand_pose_A` column ahead of `rmsd_ligand_internal_A`.
`classify_stability` is keyed on the POSE metric (thresholds unchanged:
≤2 Å stable / 2-4 drifting / >4 unstable); internal stays in the
artifact as a reported diagnostic only. `summary.json.metrics` adds
`rmsd_ligand_pose_{final,max}_a` as the primaries and renames the legacy
ligand metrics to `rmsd_ligand_internal_{final,max}_a`. Free-energy soft
gate E (`free_energy/gating.py`) reads `rmsd_ligand_pose_max_a` with a
fallback to the pre-Q6b key so pre-Q6b summaries that haven't been
re-analyzed don't error.
PART 2 (residue renumber): the MD prep (PDBFixer / OpenMM normalization)
drops resseq gap markers and renumbers chains 1..N contiguous; the
docking pipeline preserves AUTHOR numbering with the explicit gaps
(β-tubulin chain B has missing-density gaps at resseq 45-46 and
361-368). The contacts/consensus pipeline + the tubulin literature both
label by author numbering — so MD top_contacts at "LEU272 / THR273 /
ARG275" actually describe LEU275 / THR276 / ARG278 in the canonical
naming, off by a NON-constant chain-dependent gap-skip offset. New
`build_md_to_docking_resseq_map(md_pdb, docking_pdb)` builds a per-chain
`(chain, md_resseq) → docking_resseq` map (PER-CHAIN safety: a length
or residue-name disagreement on ONE chain drops that chain from the map
but lets other chains still relabel; rows for the unaligned chain pass
through unchanged in `relabel_contacts`). Called from `md/job.py` after
`per_residue_contact_frequency`, gated on `pose.receptor_pdb` being
available. Summary surfaces `receptor_renumbered: bool` provenance.
**Methods / params / versions touched:** Schema bump `1.0.0 → 1.1.0`
(additive — `summary.json.metrics` gains `rmsd_ligand_pose_*` and
renames `rmsd_ligand_*` to `rmsd_ligand_internal_*`; consumers of the
old keys need to migrate). New artifact field `receptor_renumbered: bool`
+ `reanalysis` provenance block written by the re-analyze driver. CSV
column reorder in `rmsd.csv` — `rmsd_ligand_pose_A` before
`rmsd_ligand_internal_A`. Verdict thresholds unchanged.
**Decisions + rationale:** Pose RMSD (not internal) as the verdict key
because the question MD answers is "did the ligand stay in the pocket?"
— internal RMSD is geometry-only and blind to the answer. Receptor
backbone Cα as the alignment selection because backbone is what stays
put if anything does; sidechains move too much under thermal motion to
fix the frame. Same (R, t) APPLIED (not re-fit) to the ligand because a
second alignment of the ligand would simply re-derive internal RMSD —
the whole point is to MEASURE the displacement that the first alignment
exposes. Per-chain relabel safety (not all-or-nothing) so a single
chain-A length mismatch (e.g., MD-prep truncated 39 N-term Cα residues
in some surrogate-engine fixtures) doesn't block the chain-B relabel
that contains the binding-site contacts. Renumber map built from PDB Cα
positions alone (no biopython/Bio.Align dependency) because the prep
preserves residue NAMES at corresponding positions — a name-match check
at every position serves the same role as a sequence alignment with
gap penalty ∞. Fallback in gating.py (pose key → legacy key) so
free-energy gating keeps working on pre-Q6b summaries without forcing a
re-analysis as a hard prerequisite.
**Figures / data produced (CORRECTED 2026-06-05):** The earlier framing of
this entry — "8 of 10 jobs flip verdict = the metric gaining discriminating
power" — was overstated. The corrected reading:

ONE clean demonstration of the Q6b pose metric:

  - **d34b991f (OpenMM full MD on d96d719a):** post-Q6b *internal*
    RMSD final 1.46 Å (which would be classified `stable` under the
    pre-Q6b internal-keyed rule, since 1.46 < 2.0) versus *pose* RMSD
    final 3.13 Å (`drifting` under the new pose-keyed rule). The ligand
    stayed conformationally tight (internal small) but slid ~3 Å in
    the pocket — exactly the displacement internal RMSD is blind to and
    pose RMSD captures. THE headline-positive evidence that the metric
    swap is doing real work.

Headline-positive confirmation:

  - **1f01da83 (OpenMM full MD on 4a37bb0c, taxol back-dock):** verdict
    `stable → stable`, pose RMSD 1.51 Å (matches the user's "real pose
    RMSD ~1.5 Å" prediction), internal RMSD preserved at 0.65 Å, and
    top_contacts now read PHE / LEU / THR / ARG numbered in author
    space (LEU275, THR276, ARG278, PRO274, ASP26, LEU371) — the same
    biology the pre-Q6b labels LEU272 / THR273 / ARG275 described but
    in the canonical numbering that matches the contacts/consensus
    pipeline and the β-tubulin taxane-site literature.

  - **86e5195b (OpenMM full MD, 11-frame short run on d96d719a):**
    stable → stable, pose RMSD 1.95 Å. Boundary case but on the
    correct side of the 2.0 Å threshold.

DEGENERATE (excluded from "Q6b discrimination" claim, archived under
`jobs/archive/surrogate-md/` 2026-06-05):

  - **11fa6090, 27467317, 524febf8, 592ef3e0, 6ef436cd, 987322c1,
    c2284960:** all 7 are `engine.kind="surrogate"` (RDKit-based
    fallback when OpenMM is unavailable). The surrogate pre-aligns
    each conformer to the docked pose AND holds the receptor rigid —
    removing by construction exactly the pose displacement the Q6b
    pose metric measures. Their pose-final values therefore EQUAL
    their internal-final values (3.72=3.72, 4.47=4.47) and their
    verdict flips under Q6b are NOT pose-discrimination demonstrations
    — they are da3221b catch-up. The pre-Q6b cofactor-exclusion fix
    (da3221b) only ever re-ran taxol (1f01da83); the other 9 MD jobs
    were carrying the pre-da3221b lumped cofactor RMSD until THIS
    re-analysis, so their pre→post numbers reflect both da3221b + Q6b
    arriving at the same time. With internal-only keying these 7 would
    flip identically (3.72 > 2.0 → drifting; 4.47 > 4.0 → unstable);
    Q6b adds zero discrimination on them.

  - **d34b991f's pre→post number swing** is the most extreme example
    of the same catch-up effect: pre-Q6b `lig_final=8.39 Å` was a
    pre-da3221b lumped-cofactor artifact, NOT a real ligand drift.
    The post-Q6b numbers (internal 1.46, pose 3.13) reflect da3221b's
    cofactor exclusion AND Q6b's frame swap landing at once. The
    1.46→3.13 split IS the clean pose-vs-internal demo; the
    8.39→1.46 internal reduction is da3221b finally reaching this job.

Quantitative summary (RE-FRAMED). On the 3 real-OpenMM jobs:
| Job       | PRE      | POST     | pose_final | pose_max | int_final | int_max | Q6b adds discrimination? |
|-----------|----------|----------|-----------:|---------:|----------:|--------:|---|
| 1f01da83  | stable   | stable   | 1.51       | 2.51     | 0.65      | 1.30    | Confirmation; pose≠internal but both under threshold |
| 86e5195b  | stable   | stable   | 1.95       | 1.99     | 1.13      | 1.13    | Pose-keyed result agrees with internal-keyed |
| d34b991f  | unstable | drifting | 3.13       | 3.87     | 1.46      | 2.02    | YES — internal-keyed stable, pose-keyed drifting |

Tooling outcomes: each prior summary backed up as `summary.json.pre-q6b`.
Backend pytest 88/88 PASS in 15.80s (5 new Q6b tests under
`backend/tests/test_q6b_pose_rmsd.py`: rigid-translation control proves
internal ≈ 0 AND pose ≈ |Δ|; joint-translation control proves both
metrics cancel; kabsch_fit returns identity transform on equal coords;
recovers pure translation; recovers pure 90° rotation). Free-energy
gating tests 25/25 PASS — the pose key + legacy fallback are exercised
by the existing fixtures.
**Caveat / known limitation:** PART 2's renumber map is robust to
per-chain gap-skip renumbering (the actual MD-prep behavior observed
here) but does NOT handle MID-CHAIN insertions or deletions — those
would need a real sequence-alignment with gap penalties (e.g.,
Bio.Align). The per-chain length-match + residue-name-match guard
correctly skips relabel for chains where this is the case (chain A on
11fa6090/27467317/524febf8: MD has 412 Cα, docking has 451 — those
chains' contacts retain MD numbering; the artifact's
`receptor_renumbered=True` reflects that AT LEAST ONE chain aligned).
PART 1's verdict thresholds (≤2 / 2-4 / >4 Å) are unchanged from the
pre-Q6b classifier — they were calibrated against the internal metric
but happen to be reasonable for pose RMSD too (1.5 Å is well within
"stable in pocket", 3-4 Å is the classic sliding boundary). A future
calibration pass against a labeled bound/unbound dataset is the right
follow-up; treat the current thresholds as MD-convention defaults, not
literature-derived constants. Q3 (MM-GBSA) consuming the pose-keyed
verdict can now claim "free-energy estimate gated on the ligand staying
in the pocket" — previously it could only claim "gated on the ligand
not changing internal conformation," which is a much weaker statement.

---

## 2026-06-04 — Q9 — ligand-aromatic ring perception via SMILES template; π-stacking enabled (schema 1.3.0)
**Commit(s):** `2dd4c01`
**Scope:** backend
**Supports claim(s):** — (mechanics-only enablement; see caveat)
**What changed:** Threaded the authoritative ligand SMILES (the upload-time
SMILES recorded in `docking/scores.json` / `dock_scores.json` by the docking
pipeline, read back via `read_ligand_smiles`) into the interactions detector,
combined with the obabel-converted `.pdb` of each pose (which carries explicit
CONECT records). RDKit's `AllChem.AssignBondOrdersFromTemplate` reads the
docked 3D geometry, assigns bond orders + aromaticity from the SMILES
template, and feeds the resulting bonded mol into `parse_ligand_rings`. The
two π detectors (`detect_pi_stacking`, `detect_pi_cation`) consume the
perceived rings list once per pose. New top-level `ligand_smiles_status`
field (`"authoritative" | "missing" | "failed"`) + new
`criteria.smiles_assignment` block (method + rdkit_version +
aromaticity_model + status + citation) make the path self-describing.
**Methods / params / versions touched:** Schema bump `1.2.0 → 1.3.0`
(additive — no existing field shape changes). New orchestrator kwarg
`ligand_pdb_bonded_text: Optional[str]` — when provided, used as the bonded
geometry source for ring perception (the source PDBQT lacks CONECT records
and trips RDKit's proximity-bonding over-valence check on complex ligands
like taxol; this was the Q9 root cause). When absent, the orchestrator falls
back to `ligand_text` and may yield `ligand_smiles_status="failed"` on
non-trivial scaffolds. RDKit version captured per artifact:
`2024.03.5`. Aromaticity model: `rdkit_default`
(`Chem.AROMATICITY_DEFAULT`). Method: `AllChem.AssignBondOrdersFromTemplate`.
Fail-safe contract: any template-assignment failure (atom-count mismatch,
invalid SMILES, tautomer mismatch, sanitisation error) returns
`([], "failed")` — NEVER fabricate rings from coordinates. π-stacking +
the protein-cation↔ligand-ring direction of π-cation stay empty on
"missing" / "failed". `pi.py` signature already takes `ligand_rings:
Optional[List[Ring]] = None` as its third positional (no dead SMILES param);
the orchestrator binds rings correctly.
**Decisions + rationale:** SMILES + template assignment (not coordinate-only
ring perception) because the upload-time SMILES is the only source that
KNOWS the molecule's bond orders and aromaticity — re-deriving from 3D pose
coordinates loses tautomer + protonation information that obabel's PDBQT
conversion discarded. CONECT-record requirement for the bonded mol because
RDKit's `MolFromPDBBlock` without CONECT falls back to proximity bonding,
which over-valences sp2 carbons on strained polycyclic scaffolds (taxol
fails this way). The obabel-converted `.pdb` written at viewer-export time
already carries ~620 CONECT records per taxol pose — the orchestrator now
reads that file as the bonded-geometry source. Single SMILES source: the
upload `ligands.smi` flows through `batch_docking.py` into `scores.json`
and is read by `read_ligand_smiles`; both the dock-time export path
(`core/job_artifacts.py`) and the backfill path
(`scripts/backfill_contacts_json.py`) use the same helper, so SMILES
provenance is one chain end-to-end. Ring perception centralised in the
orchestrator (not duplicated in both pi detectors) so the
`ligand_smiles_status` provenance can be stamped on the artifact in one
place. Coord-based mapping from RDKit ring atoms → ligand_df rows (3 dp
xyz match) rather than positional rank, because Vina/obabel-prepared poses
do NOT preserve SMILES atom order; the prior positional mapping was the
silent fragility this Q9 work also closes.
**Figures / data produced:** 30/30 contacts artifacts backfilled at schema
1.3.0 across `jobs/4a37bb0c-…` (15 poses) + `jobs/d96d719a-…` (15 poses),
all `ligand_smiles_status="authoritative"`. π-stacking enablement evidence:
2 poses (lig_000000 rank 02 + 04 on 4a37bb0c) surfaced 1 π-stacking finding
each — 0 across all 30 artifacts pre-Q9. Regression check on the
4a37bb0c lig_000000 rank 00 fixture: all 7 non-π keys identical pre→post
(hbond=3, weak_ch_o=2, hydrophobic=36, salt_bridge=1, halogen=0,
water_bridge=0, metal=0). Consensus re-run on 4a37bb0c pose-rank-0 best
poses (artifact at `analysis/consensus_pharmacophore/post-q9-2026-06-04/`)
identical to the pre-Q9 baseline: 4 features
(hbond_acceptor=1, hydrophobic=2, pi_cation=1, pi_stacking=0),
candidate_clusters=21, dropped_below_threshold=17, raw_feature_count=29,
per-ligand extraction.total = 43/37/15. The 4 surviving features anchor
at GLY370 (H-bond acceptor), LEU275 + PHE272 (hydrophobic), HIS229
(π-cation) — same set + same anchors. Inputs read contacts at schema 1.3.0
(vs 1.1.0 pre-Q9). Backend pytest 83/83 PASS in 20.54s (10 new Q9 tests
under `backend/tests/test_pi_smiles_template.py` covering: template
assignment yields ≥1 aromatic ring on benzene, missing SMILES → "missing"
status, parallel benzene-over-PHE control yields exactly 1 π-stacking
finding at the orchestrator level, atom-count mismatch → empty + "failed"
flag with no fabricated rings, invalid SMILES → "failed", schema version
asserts SCHEMA_VERSION, artifact carries `ligand_smiles_status` +
`criteria.smiles_assignment`).
**Caveat / known limitation:** The mechanics-only caveat from Item 23 (v1)
+ Item 25 backend extends to π-stacking: π-stacking detection is now
enabled but is POSE-DEPENDENT — it surfaces only when a ligand's aromatic
ring sits 3.5–5.5 Å from a protein PHE/TYR/TRP/HIS centroid with the right
plane-plane angle. On this fixture's best-pose set (rank 00 per ligand,
what consensus reads by default), no such geometry exists, which is why
the consensus output is unchanged. The 2 taxol non-best-pose findings
demonstrate the detector is now alive; the consensus delta would land if
a future fixture's best poses contain a stacked configuration, or if the
consensus's `pose_rank` parameter is broadened to scan multiple ranks.
The downstream Item 26 (2D ligand interaction diagram) is unaffected by
Q9 — it's the next-pass consumer of the same artifact and will inherit
the perceived rings automatically. Q9 also exercises the
`pi_stacking` `line_disc` Mol* glyph (added in Item 25 frontend); next
browser-verify pass should see those glyphs render on the lig_000000
rank 02 / 04 poses.

---

## 2026-06-16 — MM-GBSA FE sweep on the 5 corrected-gate-qualifying MD runs
**Commit(s):** `117e0cd` (engine fixes), this commit (sweep + publication notes)
**Scope:** backend (FE engine), docs (sweep output + publication registers)
**Decisions + rationale:** The bundled MM-GBSA engine
(`backend/app/free_energy/mmgbsa.py` + `mmgbsa_runner.py`) was fully
implemented + synthetic-tested but had never been exercised end-to-end on
real MD output from this repo (per `test_mmgbsa.py:10-15` — "real-data
taxol smoke is NOT a unit test"). The 5-run sweep exposed four blockers,
all fixed in `117e0cd`:

  (a) Topology source — `_build_complex_modeller` now loads `frame_000.pdb`
      directly instead of reconstructing from `docking/viewer/receptor.pdb`
      + cofactor sidecars + the ligand-with-H sidecar. The reconstruction
      route drifted from what `engine_openmm` wrote into the trajectory
      frames (N-terminal MET `H1` vs `H`; cofactor/ligand chain IDs A/B/X
      from sidecars vs C/D/E/F from the frame writer), so the
      `(chain, resname, resseq, atom_name)` strip mask mismatched 177-6500
      atoms per run. Sourcing the topology from frame 0 closes the entire
      writer-vs-builder convention-drift class by construction.
  (b) GAFF cache reuse — `_make_system_generator` now threads the MD
      pipeline's per-docking-job `_gaff_template_cache.json` through
      `SystemGenerator(cache=…)`. This sidesteps the AM1-BCC charge
      backend (AmberTools/OpenEye), which is env-dependent and was not
      registered with openff-toolkit in the runtime env. Cold ligands
      without an MD cache still need the wrapper; this is the dominant
      ~1-4 hour `sqm` step the engine's docstring warns about.
  (c) Shared SystemGenerator — the cache is delta-encoded (entry 1
      defines all GAFF atom types, later entries only define
      not-yet-registered ones), so 3 separate SystemGenerators meant the
      ligand-only build hit a fresh ForceField and KeyError'd on the
      first unknown type (`c3`, sp3 aliphatic C). One generator threaded
      through complex → receptor → ligand accumulates atom types on the
      complex build; the next two reuse the populated ForceField.
  (d) Eq-window selection — `compute_md_fe` previously fed ALL
      post-discard frames into the estimator (197/201 for explicit;
      201/201 for implicit). Aligning with the corrected gate's last-30 %
      window (`docs/all_md_corrected_gate.md`), the runner now slices
      `positions[int(0.7 * n_post_discard):]` so MM-GBSA integrates the
      same 60-61 frames the gate qualified the run on, not the
      unequilibrated head. Persisted as
      `provenance.eq_window_last_fraction` + `eq_window_skipped_frames`.

  Cache invalidation: stale `_mmgbsa_cache/` directories written by the
  pre-(a) topology code path were wiped before the final sweep (one-time
  manual `rm -rf jobs/<id>/md/_mmgbsa_cache/` per qualifying MD job);
  the runner has no cache version field, so a future schema bump should
  add one.

**Figures / data produced:** Sweep output at
`docs/free_energy_qualifying.md`; per-MD JSONs at
`jobs/<md_id>/free_energy/summary.json` (5 files, schema 1.0.0; not
committed — `jobs/*` gitignored). Methodology: single-trajectory MM-GBSA
on OBC2 implicit solvent at the energy-eval stage, force-group dispatch
(bonded / nonbonded / GBSAOBCForce solvation), force-field stack
amber14-all + amber14/tip3p ion templates + implicit/obc2 + gaff-2.11,
last-30 % eq window, configurational entropy omitted (single-trajectory
MM-GBSA convention).

ΔG_bind (kcal/mol, mean ± SEM, n=60-61 frames):

  - 5bc61f59  taxol pose 0       (explicit, 1 ns)   −58.55 ± 0.39  σ=3.00
                                        nonbonded −106.17, solvation +47.62
  - 34840aa1  taxol pose 1       (explicit, 1 ns)   −56.58 ± 0.53  σ=4.07
                                        nonbonded −120.76, solvation +64.19
  - 80e53d8a  taxol pose 2       (explicit, 1 ns)   −51.41 ± 0.70  σ=5.45
                                        nonbonded −104.03, solvation +52.62
  - a0b04941  Juliprosopine pose 0 (explicit, 1 ns) −20.32 ± 0.33  σ=2.57
                                        nonbonded  −41.95, solvation +21.63
  - 02f30602  taxol pose 0 Run A ext. (IMPLICIT, 1 ns) −43.40 ± 0.66  σ=5.12
                                        nonbonded  −87.69, solvation +44.28

Bonded ≈ 0 by construction in single-trajectory MM-GBSA (ligand bonded
energy is the same in the complex and isolated subsystems).

**Caveat / known limitation:** Each per-MD `summary.json` carries
`sampling_adequate: false` + `preliminary: true` — these flags come from
the bundled gate at `backend/app/free_energy/gating.py`, which keys on
`summary.verdict == "stable"` and so blocks every run with verdict
`drifting`. The corrected gate at `docs/all_md_corrected_gate.md`
(eq-window begin/end ≤ 0.4 Å, in-pocket ≥ 0.95, top-5 residue
persistence ≥ 65 %, ≥ 50 eq frames) is the authoritative qualification
source — it explicitly ignores `verdict` to avoid conflating drift from
the docked pose with current instability. Configurational entropy
(normal-mode / quasi-harmonic) is OMITTED; the reported ΔG is the
enthalpic + solvation part — the standard MM-GBSA quantity in the
literature. Implicit-solvent MD lineage (`02f30602`) is flagged as a
declined strategy in the gate report — its result is kept here as
evidence that explicit-vs-implicit shifts ΔG by ~15 kcal/mol on the same
pocket, NOT as an FE candidate. No replicate runs — single-replicate
MM-GBSA on a single MD trajectory; statistical uncertainty is intra-run
SEM only and doesn't bound seed-to-seed variability. AM1-BCC backend is
not env-registered, so cold-cache ligands without an existing MD-pipeline
GAFF cache would still hit the AM1-BCC `assign_partial_charges` error;
all 5 runs here hit the cache.
