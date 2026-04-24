# 04 — Limitations

## 1. Scope

### 1.1 We ran 12 of 25+ CASP15 ligand targets

Targets excluded from this work:

- **RNA-only targets** (R1117, R1117v2, R1126, R1136v1): Protenix
  supports RNA via `rnaSequence` entities, but the chai baseline
  excluded them, and we opted for a symmetric comparison over maximum
  coverage. The pipeline can handle them; re-running on the R* set is
  a ~1-hour follow-up.
- **Protein-ligand targets not in chai manifests** (T1170, T1181, T1186,
  H1114, H1171v1-v2, H1172v1-v4): excluded for the same symmetry reason.
  These include some of CASP15's harder targets and any claim of
  "Protenix performs well on CASP15" should state the subset.

### 1.2 One inference configuration

All results use:

- Model: `protenix_base_default_v1.0.0` (the newer `protenix-v2` is
  withheld by ByteDance).
- 5 seeds × 5 samples/seed = 25 structures per target.
- Templates enabled.
- No Training-Free Guidance (TFG) for ligand chirality/planarity
  refinement.

Alternatives not explored:

- More seeds (v1 paper shows log-linear improvement to 1000 seeds).
- `--use_tfg_guidance true` for physics-aware ligand sampling.
- `protenix_base_20250630_v1.0.0` (newer data cutoff).
- Ensemble selection across different models.

### 1.3 Single scoring metric

We use OST LDDT_pli + BiSyRMSD as the primary pose-quality signals.
Not included:

- **PoseBusters checks** (planarity, bond lengths, clash detection,
  inter-molecular geometry sanity).
- **Drug-discovery metrics** (Vina score, ΔG estimates).
- **Target-specific** metrics that community assessors at CASP would
  apply (e.g. for T1187 glycan: glycosidic torsion angles).

## 2. Ground-truth caveats

### 2.1 Reference SDFs are derived

`refs/<T>/lig_*.sdf` files are extracted from `<T>_lig.pdb` by
`../casp15_ligands/prep_references.py`, which reads HETATM records
through RDKit. A few gotchas to know:

- The resulting SDF's bond orders are inferred, not canonical. This
  drives the OST graph-isomorphism failures on T1158v3 / XH0.
- Partial occupancies in the crystal are collapsed to single positions.
- Alternative conformers (altlocs) default to the first one.

### 2.2 Crystal vs solution binding sites

CASP15 ligand assessors often evaluate against a single crystal
structure, but biological binding can be broader (e.g. a solvent-
exposed surface that binds under different conditions). A "wrong" pose
by RMSD may still be a biologically plausible alternative. We don't
attempt to distinguish these cases here.

### 2.3 Receptor structure inside the CIF

Protenix's model CIF contains both the receptor and the ligands. OST
uses the receptor alignment to define the protein-ligand interface for
LDDT_pli. Target-wide receptor mispredictions (rare in this set but
not zero) can drag the LDDT_pli down even when the ligand is
approximately placed.

## 3. Model-side caveats

### 3.1 Training data cutoff

`protenix_base_default_v1.0.0` is trained on data through
**2021-09-30**. CASP15 was held in **2022**. For most CASP15 targets
this is fine (the PDB depositions postdate training, so the model
hasn't seen the exact complex). For a few, closely-related structures
may have been in training — the "leakage" audit was not performed here.

### 3.2 Template DB cutoff

Template search uses `pdb_seqres_2022_09_28.fasta`. Any template newer
than September 2022 is invisible. For this CASP15 set it's not a
constraint in practice (the deposited complexes are either in the DB or
not yet deposited, both cases).

### 3.3 Confidence ≠ correctness

See [03_discussion.md §3](03_discussion.md#3-confidence-calibration).
`ranking_score` can be simultaneously high (0.94 on T1188) and the
ligands misplaced (two ions 30+ Å off). Per-chain `chain_plddt`
is more informative for ligand-specific confidence but wasn't used
in the primary best-sample selection.

## 4. Scoring caveats

### 4.1 Graph-isomorphism failure is a known OST issue

1 of 12 CASP15 and 3 of 14 CASP16 evaluated targets here produced
"0 ligand pair(s) matched" entries. Chai-lab produced the same errors
on the same targets, suggesting it's scorer-side (OST's RDKit-based
matching) not method-side.

Not a scientific limitation of Protenix — but it does mean these
targets contribute zero signal to the aggregate statistics, and the
aggregate benchmarks (e.g. "50% of targets near-native") are really
"50% of the targets we could score".

### 4.2 "Best sample" is a biased estimator

`--best-only` picks the highest-`ranking_score` sample and reports only
its score, ignoring the variance across the other 24 samples. For
reporting, this is standard (CASP uses best-of-N as well). For
statistical rigor a more careful treatment would report median + IQR
across seeds and flag cases where best and worst differ wildly.

`score_lddt_pli.py` without `--best-only` scores all 25 samples if
per-sample distributions are wanted.

## 5. Reproducibility caveats

### 5.1 MSA is stochastic-ish

Protenix's server-side MSA search hits a live MMseqs2 service whose
underlying databases are refreshed periodically. Reruns on different
dates can return slightly different hit sets, which can perturb
predictions. For this work all preps occurred within a single 24-hour
window.

### 5.2 Model checkpoint may disappear

`protenix_base_default_v1.0.0.pt` is hosted on ByteDance's Volcano
Engine TOS. The `protenix-v2.pt` URL already returns 403. If
`v1.0.0.pt` follows the same fate, this work won't be reproducible
without the user having cached the checkpoint.

**Mitigation:** `/workspace/protenix_data/checkpoint/` on the pod has
the checkpoint cached from our runs. Pulling this before tearing
down the pod preserves reproducibility. (~1 GB, worth keeping.)

### 5.3 OST version pinning

We used OpenStructure 2.x via the `ost` conda env. Different OST
versions handle `--substructure-match` with slightly different
tolerances. Results here were generated with a specific binary at
`~/anaconda3/envs/ost/bin/ost`; changes in OST may shift
LDDT_pli values by a few percent.

## 6. What this work does NOT claim

- "Protenix-v1 beats chai on CASP15." (We don't have a direct
  head-to-head here; see CASP16 head-to-head for that analysis, where
  results are mixed.)
- "Protenix-v1 matches AlphaFold3." (The v1 paper makes that claim
  against AF3 on a broader benchmark; we don't re-benchmark here.)
- Any generalization beyond the 12-target subset.
- Any claim about Protenix-v2 performance (weights are not public).
- Any clinical or drug-discovery implication — this is a methods
  validation, not a usability study.
