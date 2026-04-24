# 04 — Limitations

## 1. Scope

### 1.1 PoC subset, not full CASP16 pharma

We ran 44 of 263 pharma-ligand targets (17% coverage). The deferred
set is the L3000 Autotaxin (219 targets) batch, excluded because:

- ~15 min/target × 219 = ~55 GPU hours ≈ $47 on L40S — too large for a
  PoC with outstanding correctness concerns.
- Autotaxin's receptor is a single sequence across all 219 targets;
  `clone_prepped.py` will cache-hit heavily, so the full run would be
  feasible once we decide to do it. Not a technical blocker, a
  priorities one.

### 1.2 14 of 44 PoC targets folded pre-fix

We Ctrl-C'd the pod after L2001 and before the full 44-target run
completed, because mounting evidence suggested the homodimer bug was
corrupting roughly half the results. So the comparison tables here
represent:

- 14 targets that actually ran inference.
- 11 of those scored by OST (3 graph-iso failures in both pipelines).
- 4 spot-test re-runs of catastrophic pre-fix failures (L2001, L4001,
  L4003, L4004).

A full 44-target re-run on a fresh pod is planned. Expected ~$1.50 and
~2 hours.

### 1.3 Single model, no ensemble

We ran `protenix_base_default_v1.0.0` only. The `protenix-v2` 464M
model is withheld. Other v1 variants
(`protenix_base_20250630_v1.0.0`, `protenix_mini_*`) weren't
evaluated. An ensemble comparison would strengthen the "Protenix vs
chai" conclusion — on its own this is a v1-vs-chai comparison with
sample-of-1 model architecture on each side.

## 2. Ground-truth caveats

### 2.1 Novel ligand SDFs without CCD canonicalization

CASP16 pharma ligands don't have CCD codes. Reference SDFs are RDKit-
generated from 3D coordinates in the crystal-ligand PDB files. Bond
orders are inferred, not canonical. This drives the OST graph-iso
failures on L1006/L1007/L1008 — a scorer-side issue independent of
the model.

### 2.2 Single crystal pose per ligand

Each reference is one crystal pose. Many of these drug candidates are
flexible with multiple plausible bindings modes; a model placing the
ligand correctly in an alternate-conformation pocket would still score
poorly against the single-crystal reference. We don't attempt to
discount this.

### 2.3 Protein structure from "protein_aligned.pdb"

The CASP16 organizers provide a "pre-aligned" receptor PDB per set —
chains are oriented consistently across all targets within a set to
make cross-target comparison easier. Our builder reads this file but
doesn't distinguish "crystal receptor" from "docked receptor." If the
alignment introduced any subtle coordinate shifts, they'd propagate to
our ref. We've not checked.

## 3. Model-side caveats

### 3.1 Training data cutoff

`protenix_base_default_v1.0.0` cuts off at 2021-09-30. CASP16 was held
in 2024. Most CASP16 pharma targets postdate training, but we haven't
audited for data leakage (e.g. very similar structures in the training
set). For Chymase and Cathepsin G, these are classic drug-discovery
receptors with decades of pre-2021 PDB coverage — so *templates* are
rich even when the exact complex isn't in training.

### 3.2 Template search DB is 2022-09-28

The `pdb_seqres_2022_09_28.fasta` cutoff means any PDB deposition
after September 2022 is invisible to template search. For this set,
Mpro drug-candidate complexes deposited during COVID era are mostly
pre-2022, so the constraint is minimal in practice.

### 3.3 Single inference configuration

Documented in [01_methods.md §4](01_methods.md#4-inference). We ran
5 seeds × 5 samples × 10 recycles × 200 diffusion steps with template
ON and TFG OFF. Other configurations (e.g. TFG on for chirality
refinement, more aggressive sampling) are untested here.

## 4. Scoring caveats

### 4.1 OST graph-isomorphism failure rate

3 of 14 (21%) CASP16 targets failed OST scoring. Same 3 also fail
under chai. This is a scorer-side issue but it means our aggregate
statistics are drawn from the scored subset, not the full PoC.

For alternatives:

- **Coordinate-only RMSD after manual atom correspondence** — loses
  the chirality/symmetry-awareness of BiSyRMSD but avoids graph-iso.
- **PoseBusters checks** — chemically plausibility metrics (planarity,
  bond lengths). Doesn't compare to a reference but validates the
  pose is physically sensible.
- **Visual inspection** — appropriate for small N but doesn't scale.

We have not implemented any of these fallbacks.

### 4.2 Best-of-25 reporting

Best-sample selection on 25 samples is the best-of-25 statistic. This
is biased upward relative to the mean, and the bias isn't comparable
across methods if they produce different numbers of samples. Chai's
chai-lab output is 5 models, so a fair comparison would be "best-of-5"
from Protenix too. We haven't done this head-to-head with matched
sample budget.

Alternatively, reporting per-sample distributions (mean + IQR) avoids
this — also not done here.

### 4.3 Ligand-level vs target-level scoring

For multi-ligand CASP16 targets (e.g. Mpro with its ~307-aa chain + a
drug candidate), LDDT_pli is computed per ligand but reported per
(target, ligand) pair. Our aggregates pick the "primary ligand"
(highest-count or first-encountered), which biases toward the most
prominent ligand. Targets with ambiguous primary ligands (none in this
PoC, but possible in L3000) would need a more careful treatment.

## 5. Head-to-head caveats

### 5.1 Different MSA sources

Chai used ColabFold's `api.colabfold.com`; Protenix used ByteDance's
mirror `protenix-server.com/api/msa`. Both are MMseqs2 backends but
with potentially different databases / versions. If one has more
recent deposits or better-curated data, that's a confound not
controlled for here.

### 5.2 Different sample budgets

Chai: 5 models. Protenix: 25 (5 seeds × 5 samples). See §4.2 above.

### 5.3 Different templates policy

Chai runs here had templates OFF (the default); Protenix runs had
them ON. This is a meaningful confound if chai's wins / Protenix's
wins are driven by template availability. An apples-to-apples rerun
with matched template policy is pending.

## 6. Things this work does NOT claim

- "Protenix-v1.0.0 is better than chai-lab on CASP16 pharma ligands."
  (4 wins vs 3 is not a statistical claim; confounds unresolved.)
- "Protenix-v1.0.0 achieves AF3-class performance on drug-discovery
  benchmarks." (Would require a direct AF3 comparison we haven't done.)
- Anything about the deferred L3000 Autotaxin set.
- Anything about `protenix-v2` (weights not public).
- Anything about how these results would translate to CASP17 or
  subsequent blind benchmarks.
- Any clinical / drug-discovery implication — this is a pipeline
  validation, not a drug-discovery claim.

## 7. Known technical debt

- OST graph-iso workaround (preprocessing SDFs) not implemented.
- L3000 Autotaxin batch not run.
- Full 25-multi-chain rerun (post-fix) not completed.
- `--use_tfg_guidance` not benchmarked against default sampling.
- No attempt at per-ligand confidence via `chain_plddt` scraping.
- No seed-variance analysis (we report best-of-25 only).
