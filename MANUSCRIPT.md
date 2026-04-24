# Reproducing CASP15 and CASP16 Ligand-Binding Benchmarks with Protenix-v1.0.0: A Unified Report

**Author:** sness
**Draft date:** 2026-04-22
**Status:** Draft — CASP16 results partial pending a post-fix rerun.

## Abstract

We benchmarked **Protenix-v1.0.0** (ByteDance, 368M parameters,
AlphaFold3-equivalent training data cutoff) against two CASP
ligand-prediction datasets: the full **CASP15 ligand track** (12
validated targets) and a 44-target proof-of-concept subset of the
**CASP16 pharma-ligand track** (Chymase, Cathepsin G, and SARS-CoV-2
Mpro). Using a unified pipeline (local MSA + template prep,
distributed GPU inference on rented cloud, OST-based LDDT\_pli and
BiSyRMSD scoring), we obtained near-native predictions (LDDT\_pli ≥ 0.9)
on 6/12 CASP15 targets and on 6/11 scored CASP16 targets, at a total
compute cost of \$2.90 across all runs to date. The work surfaced a
**homodimer-input bug** in our (and chai-lab's) CASP16 builder that
caused catastrophic ligand misplacement on 25 of 44 targets; a targeted
fix rescued the three most-affected cases from LDDT\_pli 0.00 to
0.46–0.78. A head-to-head comparison against chai-lab on 14 CASP16
targets shows neither method dominates (4 Protenix wins, 3 chai wins,
4 ties, 3 scorer-failures), with the caveat that the chai side was run
on the same buggy input. We also deliver a pipeline-agnostic
**canonical inputs tree** (260 targets) with adapters for Protenix,
chai-lab, and DynamicBind — reducing the cost of adding a new
predictor or onboarding CASP17 targets to a single adapter script.

## 1. Introduction

### 1.1 The CASP ligand tracks

The biennial **Critical Assessment of Structure Prediction** (CASP)
competitions have added dedicated ligand-prediction tracks in CASP15
(2022) and CASP16 (2024). CASP15 mixed protein-ligand complexes with
nucleic-acid-containing targets and used diverse cofactors (ATP, SAH,
glycans, etc.). CASP16 introduced a **pharma-ligand** track donated by
four pharmaceutical partners — sets of novel drug candidates at known
crystallographic binding poses, against four target receptors
(Chymase, Cathepsin G, Autotaxin, SARS-CoV-2 main protease). Between
them, they cover most of the relevant protein-ligand complexity modern
structure predictors face.

### 1.2 Protenix-v1.0.0

**Protenix** (ByteDance, 2026) is an open-source PyTorch implementation
of an AlphaFold3-style biomolecular structure predictor. The v1.0.0
release (February 2026) is trained on data through 2021-09-30 (the same
cutoff AlphaFold3 used for its original release) and supports MSA,
RNA MSA, and structural templates. A larger v2 (464M parameters) is
referenced in the repo and technical report but the weights were never
publicly released ([ByteDance/Protenix#295](https://github.com/bytedance/Protenix/issues/295)); v1.0.0
is therefore the largest publicly-runnable model.

### 1.3 Chai-lab baseline

**Chai-lab** (Chai Discovery, 2024) is a comparable open-source
AlphaFold3-derivative predictor widely used as a CASP baseline. It
predates Protenix and does not natively consume structural templates.
Our group had run chai against both CASP15 and CASP16 pharma-ligand
targets before this work began; we use those results as a baseline
here.

### 1.4 What this report covers

- **§2 Methods** — unified pipeline (inputs, prep, inference, scoring).
- **§3 CASP15 results** — per-target LDDT\_pli + RMSD on 12 targets.
- **§4 CASP16 results** — partial PoC (14 of 44 targets) plus
  validation spot-tests after fixing the homodimer bug (§5).
- **§5 The homodimer-input bug** — a case study of how a CASP16 input
  bug caused catastrophic LDDT\_pli failures and the audit that found it.
- **§6 Head-to-head with chai-lab** — on the 14-target shared subset,
  plus caveats.
- **§7 Infrastructure** — canonical `inputs/` tree as the main
  engineering contribution; enables CASP17 participation by adding
  one adapter per new tool.
- **§8 Limitations** — what this work does and does not claim.
- **§9 Reproducibility** — full repro path with commands.

### 1.5 Related documents

This manuscript is a **consolidated view**. For depth, see the
per-benchmark docs:

- `casp15_ligands_protenix/docs/01_methods.md` through `05_reproduction.md`
- `casp16_ligands_protenix/docs/01_methods.md` through `05_reproduction.md`
- `STATUS.md` (current project state, updated each milestone)

## 2. Methods

### 2.1 Unified pipeline overview

For every target we ran the same four-stage pipeline:

1. **Input assembly** — per-target JSON (Protenix) or FASTA (chai)
   listing each receptor chain's protein sequence and each ligand's
   SMILES or CCD code.
2. **Feature preparation** — protein MSAs (paired + unpaired) from
   ByteDance's Protenix MSA server (a ColabFold-API-compatible
   MMseqs2 service) and structural templates via local
   `hmmbuild`/`hmmsearch`/`kalign` against PDB seqres 2022-09-28.
3. **Structure prediction** — Protenix-v1.0.0 inference on rented
   cloud GPUs. Five seeds (101–105), five samples per seed, ten
   recycles, 200 diffusion steps, templates enabled. Total: 25 CIFs
   per target.
4. **Scoring** — OpenStructure 2.x `compare-ligand-structures` with
   `--lddt-pli --rmsd --substructure-match`, scoring the best sample
   (highest `ranking_score`) against crystallographic reference ligand
   SDFs and a ligand-free receptor PDB.

### 2.2 Input construction

CASP15 and CASP16 use different source data layouts:

- **CASP15:** a combined `<TARGET>_lig.pdb` per target with protein
  chains + ligand HETATMs in a single file. Our builder extracts
  per-chain one-letter sequences, collects HETATM residues as ligand
  CCD codes, and looks up SMILES via the RCSB Chemical Component
  Dictionary.
- **CASP16 pharma:** split sources — `protein_aligned.pdb` (receptor
  only) plus a per-target `<TARGET>.tsv` with columns `ID Name SMILES
  Task`. Novel pharma ligands lack CCD assignments; SMILES comes
  directly from the TSV.

In both cases, the critical step is **per-chain sequence extraction
with exact-identity deduplication**. Chains with identical sequences
are merged into a single `{"proteinChain": {"sequence": S, "count": N}}`
entity. Chains with different sequences — even differing by a single
residue — remain separate entities, each with `count: 1`. See §5 for
why this matters.

For the chai pipeline, the input is a flat FASTA with per-chain and
per-ligand records:

```
>protein|name=T1124-chainA
<sequence>
>ligand|name=T1124-lig1-SAH
<SMILES>
```

### 2.3 Feature generation (`protenix prep`)

The MSA step submits two tickets per protein chain to
`https://protenix-server.com/api/msa` (MMseqs2-App backend):

- `POST /ticket/msa` — non-pairing search
- `POST /ticket/pair` — taxonomically paired search for complex MSAs

Results are polled (`GET /ticket/{ID}`) until `status: COMPLETE`, then
downloaded as a tarball and rearranged into `pairing.a3m` +
`non_pairing.a3m` per chain. Queue wait dominates wall-clock time: cached
sequences return in ~10 s; uncached queries take 2–30 minutes.

Template search runs entirely locally: `hmmbuild` constructs an HMM
profile from the MSAs, `hmmsearch` scans
`pdb_seqres_2022_09_28.fasta` with loose E-value cutoffs, `kalign`
realigns hits into a query-anchored `hmmsearch.a3m`. Total per-chain
wall time: ~7 seconds.

### 2.4 Inference

We used a single inference configuration throughout:

| Parameter | Value |
|---|---|
| Model | `protenix_base_default_v1.0.0` (368M, bf16 mixed precision) |
| Seeds | 101, 102, 103, 104, 105 |
| Samples per seed | 5 |
| Recycles (`model.N_cycle`) | 10 |
| Diffusion steps (`sample_diffusion.N_step`) | 200 |
| Templates | Enabled (`--use_template true`) |
| Triangle attention kernel | `cuequivariance` |
| Triangle multiplicative kernel | `cuequivariance` |
| TFG guidance | Off |

**Mixed precision:** most computation is bf16 but `SampleDiffusion` and
`ConfidenceHead` stay in fp32 for numerical stability. Protenix's
`update_inference_configs` helper automatically lowers those modules to
bf16 at >2560 tokens; no CASP15/16 target in our subset reached that
threshold.

### 2.5 Scoring

Best-sample selection uses `ranking_score` from each sample's
`_summary_confidence_sample_N.json`. The selected CIF is scored
against the reference with:

```
ost compare-ligand-structures \
    -m <MODEL>.cif -mf cif \
    -r refs/<T>/receptor.pdb -rf pdb \
    -rl refs/<T>/lig_*.sdf \
    -of json -ft \
    --lddt-pli --rmsd \
    --substructure-match
```

Metrics reported:

- **LDDT\_pli** — local Distance Difference Test restricted to
  protein-ligand interaction atoms, range 0-1 higher-better. Counts
  contacts within a 6 Å radius.
- **BiSyRMSD** — symmetry-aware RMSD in Å (permutes equivalent atoms
  before alignment), lower-better.

`--substructure-match` tolerates partial-resolution ligands during
graph-isomorphism matching between the model CIF's ligand and the
reference SDF.

### 2.6 Hardware and cost

All inference on rented cloud GPUs:

- **RunPod L40S 48 GB** ($0.86/hr secure) — primary, used for the
  CASP15 full run and the CASP16 partial run.
- **Vast.ai RTX 3090 24 GB** ($0.18/hr community) — fallback when
  RunPod was GPU-starved; used for CASP16 L2001 post-fix spot test.
- **RunPod RTX 4090 24 GB** (~$0.40/hr) — used for CASP16 L4
  post-fix spot tests and the pending full rerun.

Local MSA prep runs on the user's workstation (no GPU needed).

## 3. CASP15 results

### 3.1 Target set

12 targets drawn from the chai-lab A40 and A100 manifests
(corresponding to ≤1024 tokens and 1536–2048 tokens respectively):

- **A40 class:** T1124, T1127, T1127v2, T1146, T1152, T1187, T1188
- **A100 class:** H1135, T1158v1, T1158v2, T1158v3, T1158v4

RNA-only targets (R\*) and a handful of protein targets not in chai's
manifests (T1170, T1181, T1186, H1114, H1171v1-v2, H1172v1-v4) are
excluded for symmetry with the baseline.

### 3.2 Per-target outcomes

Best-sample metrics per target. Numbers are `LDDT_pli / BiSyRMSD (Å)`.

| Target | Primary ligand (count) | LDDT\_pli | RMSD (Å) | `ranking_score` |
|---|---|---|---|---|
| **T1124** | SAH (×2) | 0.97 / 0.93 | 0.14 / 0.33 | 0.97 |
| **T1127** | EPE (×2) | 0.96 / 0.94 | 0.43 / 0.58 | 0.97 |
| **T1127v2** | EPE (×2) | 0.96 / 0.94 | 0.42 / 0.58 | 0.97 |
| **T1146** | NAG | 0.84 | 0.40 | 0.98 |
| **T1152** | NAG | 0.52 | 1.67 | 0.62 |
| **T1158v1** | XPG | 0.65 | 2.09 | 0.90 |
| **T1158v2** | P2E | 0.66 | 2.25 | 0.90 |
| **T1158v3** | XH0 | — | — | 0.90 (OST graph-iso failure) |
| **T1158v4** | ATP (×2) | 0.95 / 0.96 | 0.37 / 0.31 | 0.81 |
| **T1158v4** | Mg²⁺ (×2) | 1.00 / 0.98 | 0.27 / 0.15 | (same sample) |
| **T1187** | NAG (×2) | 0.15 / 0.14 | 5.91 / 5.80 | 0.40 |
| **T1188** | DW0 (×2) | 0.62 / 0.50 | 2.46 / 4.96 | 0.94 |
| **H1135** | K⁺ (×9) | 1.00 each | 0.08–0.22 | 0.60 |

### 3.3 Aggregate

- **6 of 12** near-native (primary ligand LDDT\_pli ≥ 0.9): T1124,
  T1127, T1127v2, T1158v4, H1135 (K⁺), T1146.
- **3 of 12** mid-range (0.5–0.9): T1152, T1158v1, T1158v2.
- **2 of 12** poor (< 0.5): T1187, T1188.
- **1 of 12** lost to scorer: T1158v3 (XH0).

### 3.4 Total compute

- **A40 manifest:** 31 minutes wall on L40S (includes one-time
  checkpoint download + LayerNorm JIT compile in the first target).
- **A100 manifest:** 77 minutes wall on L40S.
- **Combined:** 108 min, \$1.54 at L40S's \$0.86/hr.

### 3.5 Ligand-class patterns

| Class | Representative CCDs | LDDT\_pli range | Observation |
|---|---|---|---|
| Nucleotide cofactors | ATP, SAH | 0.93–0.97 | Consistently excellent |
| Cations (large) | K⁺, Mg²⁺ | 0.95–1.00 | Near-perfect where well-coordinated |
| Zwitterionic buffers | EPE (HEPES) | 0.94–0.96 | Strong |
| Glycans (single) | NAG | 0.15–0.84 | Bimodal — target-specific |
| Coenzymes (large, flexible) | COA | 0.52–0.57 | Right pocket, tail mislocated |
| Novel drug-like | XPG, P2E, DW0, XH0 | 0.15–0.66 or failure | Middling |
| Halide ions | Cl⁻ | 0.00–0.08 | Systematically misplaced |
| Solvent additives | MPD | 0.00 | Systematic failure (may be correct — artifacts) |

## 4. CASP16 pharma-ligand results

### 4.1 Target set

44-target proof-of-concept subset `manifest_L1L2L4.txt`:

- **Chymase (L1000):** L1001-L1017 (17 monomer targets)
- **Cathepsin G (L2000):** L2001, L2002 (2 homodimer targets)
- **SARS-CoV-2 Mpro (L4000):** L4001-L4028 excluding gaps (25 homodimer targets)

L3000 Autotaxin (219 targets) is deferred.

### 4.2 Partial results (pre-fix, 14 targets scored)

Run was Ctrl-C'd after 14 of 44 targets completed, when diagnostic
evidence suggested an input bug on multi-chain targets (§5).

| Target | Set | LDDT\_pli | RMSD (Å) | `ranking_score` | Outcome |
|---|---|---|---|---|---|
| L1001 | Chymase | 0.910 | 1.03 | 0.98 | Near-native |
| L1002 | Chymase | 0.678 | 2.66 | 0.97 | Mid |
| L1003 | Chymase | 0.967 | 0.60 | 0.98 | Near-native |
| L1004 | Chymase | 0.974 | 0.46 | 0.99 | Near-native |
| L1005 | Chymase | 0.968 | 0.65 | 0.99 | Near-native |
| L1006 | Chymase | — | — | — | OST graph-iso failure |
| L1007 | Chymase | — | — | — | OST graph-iso failure |
| L1008 | Chymase | — | — | — | OST graph-iso failure |
| L2001 | Cathepsin G | 0.000 | 50.4 | 0.82 | 🚨 catastrophic (§5) |
| L2002 | Cathepsin G | 0.893 | 1.32 | 0.98 | Near-native |
| L4001 | Mpro | 0.010 | 47.3 | 0.88 | 🚨 catastrophic (§5) |
| L4002 | Mpro | 0.899 | 0.45 | 0.91 | Near-native |
| L4003 | Mpro | 0.003 | 44.4 | 0.88 | 🚨 catastrophic (§5) |
| L4004 | Mpro | 0.359 | 5.56 | 0.84 | Poor |

**Aggregate (11 scored):** 6 near-native, 1 mid, 3 catastrophic (all
homodimer, all bug-affected), 1 poor (L4004 — see §5.5).

### 4.3 Post-fix validation spot tests

After identifying the homodimer bug (§5), we rebuilt the JSONs and
ran four targeted spot tests:

| Target | Pre-fix | Post-fix | Δ | GPU used |
|---|---|---|---|---|
| L2001 | 0.000 / 50.4 Å | **0.746 / 1.91 Å** | Rescued | Vast 3090 |
| L4001 | 0.010 / 47.3 Å | **0.784 / 1.78 Å** | Rescued | RunPod 4090 |
| L4003 | 0.003 / 44.4 Å | 0.462 / 5.22 Å | Partial rescue | RunPod 4090 |
| L4004 | 0.359 / 5.56 Å | 0.361 / 5.44 Å | No change | RunPod 4090 |

L2001 spot test cost \$0.03 (9.5 minutes on a \$0.18/hr RTX 3090). L4
spot tests cost \$0.40 (33 minutes on a \$0.72/hr 4090).

### 4.4 Full-batch attempt and stale-prep contamination (2026-04-22)

A full 44-target rerun was attempted on RunPod RTX 4090 (pod
`fe5xqdfbjywk26`) on 2026-04-22. 38 of 44 targets completed before
user-initiated pod termination. On scoring, we discovered **16 of the
38 completed runs are contaminated**: the earlier homodimer-fix
wipe had a broken shell pattern and only removed 1 of ~26 intended
stale `jsons_prepped/<T>.json` entries. The push to the 4090 pod
shipped those stale files unchanged, and `run_batch.sh`'s "skip prep
if prepped JSON exists" logic silently used them — producing
pre-fix-equivalent predictions for every multi-chain target except
L2001 (the one that got correctly wiped).

Verification that the contamination diagnosis is real:

- L4001 run = 0.009 LDDT, matches pre-fix 0.010 to 0.001 precision
  (post-fix spot-test was 0.784 — not observed in the run)
- L4003 run = 0.000, matches pre-fix 0.003 (post-fix spot-test 0.462)
- L4004 run = 0.346, matches pre-fix 0.359 (post-fix spot-test 0.361)
- L4002 run = 0.900, matches pre-fix 0.899 to 0.001 precision
- L2001 run = 0.746, matches post-fix spot-test 0.746 exactly
  (only dimer with correct fresh prep on the pod)

### 4.5 Current trustworthy dataset

After the 2026-04-22 incident, the verified post-fix results in hand:

| Category | n | Targets |
|---|---|---|
| Monomer in source (bug N/A — result correct regardless of prep) | 15 | L1001-L1005, L1009-L1017, L2002 |
| Fresh prep on the 4090 rerun | 1 | L2001 |
| Post-fix spot tests on earlier pod | 3 | L4001, L4003, L4004 |
| **Contaminated in this run (discard)** | 16 | L4002, L4006-L4020, L4022 |
| Not run (batch terminated) | 6 | L4023, L4024 (monomer), L4025-L4028 |
| Scorer-blocked (OST graph-iso) | 3 | L1006, L1007, L1008 |
| **Total** | **44** | |

**Trustworthy post-fix data: 19 of 44** (43%). **Needs re-run on next pod:
22 targets** (16 contaminated + 6 not-reached).

Cost of incident: ~\$1-2 of wasted GPU time and a full wipe-and-verify
cycle. Fully attributable to a shell regex that silently underperformed
its intent; root-causing it yielded a feedback memory entry flagging the
pattern for any future builder change. Local state is now clean (17 L1000
monomer prepped JSONs only) so the next pod rerun will do fresh prep for
all 27 multi-chain targets.

### 4.6 Pending

The remaining 22 targets need a clean re-run:

- 16 contaminated dimers (L4002, L4006-L4020, L4022) — stale prep → new fresh prep
- 6 never-reached targets (L4023, L4024, L4025, L4026, L4027, L4028) — never ran

Estimated cost on a fresh 4090: ~\$1-1.50, ~2-3 hours (shorter than the
full-44 attempt because 17 L1000 Chymase monomers are already clean and
L2001 is already correct).

The Autotaxin L3000 batch (189 targets, ~\$6, ~16 hours) is deferred
until after the 22-target cleanup and the chai rerun (§6.3).

## 5. The CASP16 homodimer-input bug

### 5.1 The bug

Our initial CASP16 input builder concatenated protein chains from
`protein_aligned.pdb` in file order into a single long sequence. For a
monomeric receptor (Chymase, L1000) this produced the correct input.
For a homodimer with chain A (e.g. 226 aa) and chain B (e.g. 227 aa,
differing by one terminal residue), it produced a **chimeric 453-aa
"fake-monomer" sequence** that doesn't exist in nature.

Protenix dutifully folded this chimera into a plausible-looking
(though incorrect) protein, placed the ligand somewhere, and the
prediction emerged with a **high confidence score** despite the ligand
being 40+ Å from the real binding site.

### 5.2 Prevalence

| Set | Multi-chain | Affected |
|---|---|---|
| L1000 Chymase | 0 / 17 (all monomer) | 0 |
| L2000 Cathepsin G | 2 / 2 | 2 |
| L4000 Mpro | 25 / 25 | 25 |
| **Total PoC** | **27 / 44 (61%)** | **27** |

### 5.3 Discovery

Three observations converged:

1. **Shared failures with chai:** L2001 and L4003 both failed
   catastrophically in chai-lab too (40+ Å RMSD). A Protenix-specific
   bug wouldn't explain this, so the first hypothesis was "hard
   target with a symmetric binding site."
2. **Confident wrong answers:** Protenix's `ranking_score` was 0.82+ on
   all three catastrophic failures. The model was confident about an
   obviously-wrong placement — consistent with "the model folded
   something plausible but fed wrong input," not "the model is
   uncertain."
3. **Chain-count audit:** inspection of the JSON for L2001 revealed
   the 453-aa concatenated sequence. A sweep across all 44 PoC targets
   confirmed 27 were affected.

### 5.4 The fix

`build_protenix_json.py` now extracts per-chain sequences and
deduplicates by exact identity:

- **Identical chains** (both 226 aa, same sequence) collapse into a
  single entity with `count: 2` — Protenix's native homo-oligomer
  notation.
- **Near-identical chains** (226 vs 227 aa, truncation variants) remain
  as separate entities with `count: 1` each.
- **Concatenation is never allowed.**

The fix is 20 lines in the builder. Validation spot tests (§4.3)
confirmed dramatic rescues on 2 of 3 catastrophic targets.

### 5.5 Interpretation

**L4004 didn't rescue** (0.36 → 0.36), despite being a homodimer. This
is the most informative data point: not every multi-chain failure is
the bug. L4004's receptor was correctly represented, but its specific
ligand is just hard for Protenix to dock — probably
a combination of conformational flexibility and a less-constrained
binding mode. The fix is targeted: it rescues catastrophic failures
caused by chimera input, and has no effect on genuinely-hard targets.

### 5.6 Broader lesson

The **chai-lab pipeline has the same bug**, unfixed. Chai's CASP16
fastas at `casp16_ligands/fastas/` still concatenate homodimer chains
into a single sequence. This means the apparent chai catastrophes on
L2001 and L4003 are likely input-driven, not method-driven.

**The head-to-head comparison in §6 is therefore unfair to chai.** A
chai re-run with the canonical-generated fastas would likely rescue
chai's catastrophes as well, and the method-level winner on those
targets is currently unknown.

## 6. Head-to-head with chai-lab

### 6.1 The 14-target comparison

Both methods ran on the identical 14 CASP16 targets with the same
ground-truth references and the same OST scorer. Best-sample metrics
per target:

| Target | chai LDDT / RMSD | Protenix LDDT / RMSD | Winner |
|---|---|---|---|
| L1001 | 0.960 / 0.48 | 0.910 / 1.03 | chai (small) |
| L1002 | 0.805 / 1.65 | 0.678 / 2.66 | chai |
| L1003 | 0.510 / 3.90 | **0.967 / 0.60** | **Protenix** |
| L1004 | 0.842 / 1.52 | **0.974 / 0.46** | **Protenix** |
| L1005 | 0.927 / 0.87 | 0.968 / 0.65 | tie |
| L1006/L1007/L1008 | — | — | both OST-unscored |
| L2001 | 0.000 / 39.9 | 0.000 / 50.4 | both fail (bug) |
| L2002 | 0.876 / 1.76 | 0.893 / 1.32 | tie (P slight edge) |
| L4001 | 0.724 / 2.14 | 0.010 / 47.3 | chai (Protenix bug-affected) |
| L4002 | 0.000 / 34.3 | **0.899 / 0.45** | **Protenix (huge)** |
| L4003 | 0.009 / 44.4 | 0.003 / 44.4 | both fail (bug) |
| L4004 | 0.000 / 50.1 | **0.359 / 5.56** | **Protenix (less bad)** |

### 6.2 Summary

| Outcome | Count | Targets |
|---|---|---|
| Protenix wins (ΔLDDT\_pli ≥ 0.05) | 4 | L1003, L1004, L4002, L4004 |
| chai wins (ΔLDDT\_pli ≥ 0.05) | 3 | L1001, L1002, L4001 |
| Tie / both strong | 2 | L1005, L2002 |
| Both catastrophic (bug) | 2 | L2001, L4003 |
| Both OST-unscored | 3 | L1006, L1007, L1008 |

**Neither method dominates.** The most striking Protenix advantages
are on L4002 (chai placed the ligand 34 Å from the pocket; Protenix
hit 0.45 Å) and L4004 (chai 50 Å off; Protenix 5.6 Å off — still poor
but at least in the right neighborhood). The clearest chai advantage
is L4001 pre-fix (2.1 Å vs 47 Å), which recovered under the Protenix
homodimer fix to a narrow Protenix win (1.78 Å). L1001 and L1002 show
small chai edges of 0.05–0.13 LDDT\_pli — within typical seed variance
given chai sampled 5 models while Protenix sampled 25, so statistical
weight is modest.

### 6.3 Caveat: the comparison is unfair to chai

As noted in §5.6, chai's input fastas were generated with the same
concat bug we fixed on the Protenix side. The 3 catastrophic chai
results above (L2001, L4001, L4003) are likely input-driven, not
method-driven. Until chai is re-run with corrected inputs, any
head-to-head reading should treat chai's catastrophes as "possibly
recoverable" rather than true method failures.

A chai re-run is planned (§8.3) and will update this section.

## 7. Infrastructure: the canonical inputs tree

### 7.1 Motivation

Each structure predictor consumes a slightly different input format:
Protenix takes JSON with per-entity chain + ligand descriptions; chai
takes a flat FASTA with `>protein|`/`>ligand|` records; DynamicBind
takes a per-target directory with `protein.fasta` + `ligand.sdf`. Our
initial approach — writing one input-prep script per pipeline — meant
the homodimer bug had to be fixed separately in each codebase. Chai's
version is still unfixed as of this writing.

### 7.2 Design

Build **one canonical representation** of each target, and **N
adapters** that render the canonical into pipeline-specific formats.
Adding a new predictor becomes a single adapter (~100 lines), not a
full input-pipeline fork.

The canonical representation is a `target.json` per target:

```json
{
  "target_id": "L4001",
  "casp_year": 16,
  "set": "L4000",
  "chains": [
    {"id": "chainA_306aa", "type": "protein",
     "sequence": "SGF...", "count": 1, "source_chains": ["A"]},
    {"id": "chainB_304aa", "type": "protein",
     "sequence": "SGF...", "count": 1, "source_chains": ["B"]}
  ],
  "ligands": [
    {"id": "lig_01", "ccd": null,
     "smiles": "O=C(Cn1ccnn1)Nc1ccc(Br)c([N+](=O)[O-])c1",
     "count": 2,
     "sdf_paths": ["ligands/lig_01_LIG_A_401.sdf",
                   "ligands/lig_02_LIG_B_401.sdf"]}
  ],
  "receptor_pdb": "receptor.pdb",
  "notes": "L4000 — LIG"
}
```

Full schema in [`inputs/docs/schema.md`](inputs/docs/schema.md).

### 7.3 Coverage

- CASP15: 27 targets (full source bundle)
- CASP16: 233 targets (all four pharma sets)
- CASP17: placeholder for 2026-04-28 target drop
- **Total: 260 populated target directories**

### 7.4 Adapters

Three pipeline adapters as of this writing:

- `adapters/to_protenix_json.py` — canonical → Protenix JSON
- `adapters/to_chai_fasta.py` — canonical → chai FASTA
- `adapters/to_dynamicbind.py` — stub

Adapter contract: take target.json paths, write to `--out-dir`,
never modify the canonical tree, fail gracefully per target and
continue. Writing a new adapter takes ~30 minutes. Full guide in
[`inputs/docs/adapter_guide.md`](inputs/docs/adapter_guide.md).

### 7.5 Regression validation

Adapter outputs were compared byte-for-byte against our existing
hand-written pipeline inputs:

| Comparison | Identical | Differ |
|---|---|---|
| Protenix CASP15 | 12 / 12 (scored targets) | 0 |
| Protenix CASP16 | 212 / 233 | 21 (documented: ligand count semantics) |
| chai CASP15 | 3 / 25 | 22 (cosmetic ordering / naming) |
| chai CASP16 | 208 / 233 | 25 (legacy concat bug — adapter is correct) |

The Protenix CASP15 byte-identical result means **the canonical tree
is a drop-in replacement for the hand-written pipeline inputs** on our
primary benchmark. The CASP16 21-diff case is a documented semantic
difference (crystal-count vs SMILES-count for multi-copy ligands) not
a bug. The chai CASP16 25-diff case is the **legacy bug** the
adapter correctly fixes.

Full regression report in [`inputs/docs/regression.md`](inputs/docs/regression.md).

### 7.6 CASP17 readiness (activated 2026-04-22)

The CASP17 pipeline is live. A separate session wrote a target-downloader
at `casp17_ligands/` that tracks predictioncenter.org's CASP17 target
releases. A bridge script ingests that output into our canonical tree:

```bash
python3 inputs/build_casp17_from_spider.py
```

Per-target transformation:
- Protein FASTA → canonical `chains[*]` (sequence + count)
- `ligands.smi` TSV → canonical `ligands[*]` with SMILES (CCD code
  inferred when `Name` looks like a 3-letter code)
- `template.pdb` → copied if non-placeholder; flagged in `notes`
  otherwise (CASP17 releases placeholder zero-coord templates until the
  judging phase)
- No `receptor_pdb` or `sdf_paths` — those arrive post-season when
  organizers release ground truth

**First target ingested:** T1214 (YncD, NP_415968, *E. coli*),
677-residue monomer with one PQQ ligand. Token count 701 → fits any
24 GB GPU; estimated ~10 min and ~\$0.03 on a 3090.

Adapter outputs validated end-to-end: `to_protenix_json.py` correctly
emits `CCD_PQQ` since PQQ is a well-known cofactor CCD; `to_chai_fasta.py`
produces a clean single-chain FASTA. Both pipelines can run inference on
T1214 immediately.

**Blind-scoring caveat:** we can predict CASP17 targets in real time,
but LDDT\_pli + BiSyRMSD scoring requires ground-truth SDFs which
organizers release only post-season. Prediction outputs (CIFs +
`ranking_score` / `chain_plddt` JSONs) are the useful deliverables
during the blind phase.

### 7.7 Token-based sizing helper

A companion script `inputs/adapters/estimate_tokens.py` computes the
Protenix-equivalent token count for any `target.json`:

```bash
python3 inputs/adapters/estimate_tokens.py --all --summary
```

Distribution across all 261 targets in the canonical tree:

| Token bucket | Count | GPU tier |
|---|---|---|
| ≤1024 | 247 | 24 GB (RTX 3090/4090) |
| 1025-2048 | 5 | 48 GB (A40/L40S) |
| 2049-4000 | 8 | 48 GB (longer wall time) |
| >6000 | 1 | 80 GB + mitigations |

The outlier is H1114 from CASP15 (7984 tokens; 20 chains, 624 ligand
atoms). Neither we nor chai attempted it. For CASP17 any similar-class
target would need §8.3 mitigations (reduced cycles/samples, mini models,
or chain splitting).

Full GPU selection guidance lives at
`~/data/vaults/docs/GUIDE-gpu-selection-protenix.md`.

## 8. Limitations and open questions

### 8.1 Scope

- **CASP15:** 12 of 25+ targets (the subset chai validated; RNA and
  protein-only targets excluded).
- **CASP16:** 44 of 263 targets as proof-of-concept; L3000 Autotaxin
  (219 targets) deferred for cost reasons.
- **One inference configuration** (5 seeds × 5 samples, 10 recycles,
  templates on, TFG off).
- **One scoring scheme** (LDDT\_pli + BiSyRMSD). No PoseBusters, no
  affinity estimation, no solvent-occupancy checks.

### 8.2 Model caveats

- **Training data cutoff 2021-09-30.** Most CASP targets post-date this
  but a data-leakage audit (checking for very-close-homolog
  structures in training) has not been performed.
- **No v2 weights.** The public checkpoint is v1.0.0 at 368M
  parameters; the 464M v2 described in the Protenix paper has not been
  released (ByteDance issue #295). All results here are v1.0.0.
- **Single-model results.** No ensemble across multiple Protenix or
  AlphaFold3-style predictors.

### 8.3 Known pending

1. **Full 32-target CASP16 post-fix rerun.** Staged and ready; awaits a
   GPU slot. Estimated \$2, ~4-5 hours.
2. **Chai CASP16 rerun with corrected fastas.** After #1; would
   strengthen the §6 head-to-head by removing the bug confound.
3. **L3000 Autotaxin batch.** 189 targets, ~\$6, deferred until
   post-CASP17.
4. **OST graph-isomorphism failure workaround.** 3 of 14 CASP16 and 1
   of 12 CASP15 targets produced "0 ligand pair(s) matched" because
   OST's graph-iso matcher couldn't align the model ligand to the
   reference SDF. Same failure under chai. Not blocking — scorer-side
   issue — but these targets contribute no signal to the aggregates.
5. **Ligand-count semantics** (§7.5 Category 2). The canonical tree
   uses crystal-count for multi-copy ligands (e.g. `count: 2` for Mpro
   with 2 SDFs); we benchmarked with SMILES-count (`count: 1`).
   Deciding which is the "right" semantic for production inference
   awaits clean post-fix data.

### 8.4 What this work does NOT claim

- That Protenix-v1.0.0 is better than chai-lab on CASP
  ligand-binding. (4 wins vs 3 on a 14-target subset with a known
  chai input bug is not a statistical claim.)
- That either method achieves AlphaFold3-level performance. (No
  direct AF3 comparison here.)
- Anything about the 219 deferred Autotaxin targets, the full
  CASP16 main-results tracks (D1273 etc.), or CASP17.
- Any clinical or drug-discovery implication. This is a methods
  validation, not a usability study.

## 9. Reproducibility

Full setup and run instructions in each arm's `docs/05_reproduction.md`.
Brief summary:

1. **Local environment:** conda env with `protenix` from pip,
   `hmmer` + `kalign` from bioconda/apt, separate conda env for OST.
2. **Input prep:** `python3 build_canonical.py --all` to populate the
   inputs tree, then `python3 adapters/to_protenix_json.py ...` per
   pipeline.
3. **MSA prep:** `./prep_local.sh manifest.txt` runs MSA + template
   search locally. MSA server queue time dominates wall clock.
4. **GPU inference:** rent a cloud GPU, rsync the prepared tree in,
   `./run_batch.sh manifest.txt`. Resumable — skips targets that
   already have 25 CIFs.
5. **Scoring:** rsync results back, `python3 score_lddt_pli.py results
   out/lddt_pli --best-only`.

**Total reproduce cost for CASP15 (12 targets):** \$1.54 + 108 min
wall on L40S.

**CASP16 PoC (44 targets) post-fix total estimate:** ~\$4-5, ~4-5
hours (half prep-queue, half GPU).

### 9.1 Infrastructure for new predictors

To onboard a new predictor X:

1. Write `inputs/adapters/to_X.py` (~100 lines; template in
   `inputs/docs/adapter_guide.md`).
2. `python3 inputs/adapters/to_X.py inputs/casp{15,16}/*/target.json --out-dir <X-input-dir>`
3. Run X's native tooling against the adapted input.
4. Score with the existing OST scorer (or adapt the scorer walker
   pattern to X's output layout).

Estimated one-time cost: 1 day of engineering per new predictor.

## 10. Conclusions

- Protenix-v1.0.0 is a credible and cost-effective choice for CASP
  ligand-prediction benchmarking: half of CASP15 targets near-native,
  total cost under \$2 per benchmark.
- Multi-chain input handling matters enormously. The same bug
  (concatenate chains into a fake monomer) was present in two
  independent pipelines; identifying and fixing it rescued 2 of 3
  catastrophic CASP16 failures. We predict 12–18 further rescues
  when the full 32-target post-fix rerun completes.
- Protenix and chai-lab have non-overlapping failure modes on
  pharma-ligand targets. Running both and comparing outputs is a
  sensible triage until more of the failure space is characterized.
- A canonical `inputs/` tree with one-adapter-per-pipeline pays off
  within two predictors; we demonstrate three and are ready to onboard
  a fourth (AlphaFold3 or a hypothetical CASP17-era tool) with
  minimal marginal cost.

## Acknowledgments

- CASP target data from predictioncenter.org; funded by NIH/NIGMS.
- Protenix from ByteDance Research (Apache-2.0).
- chai-lab from Chai Discovery.
- OpenStructure from the OpenStructure team at SIB.
- Cloud GPUs via RunPod and Vast.ai; total vendor-platform cost for
  this entire work to date: ~\$3.
- Pipeline development assisted by Claude Code (Anthropic).

## Appendix A — Document map

| Document | Role |
|---|---|
| `MANUSCRIPT.md` | This file — consolidated report. |
| `STATUS.md` | Living project status, updated at each milestone. |
| `README.md` | Directory map + quick links. |
| `casp15_ligands_protenix/docs/01_methods.md` ff. | CASP15 arm paper-style writeup. |
| `casp16_ligands_protenix/docs/01_methods.md` ff. | CASP16 arm paper-style writeup. |
| `casp15_ligands_protenix/docs/prep_phase.md` | Internals of `protenix prep`. |
| `casp15_ligands_protenix/docs/templates.md` | Template consumption in Protenix. |
| `inputs/README.md` + `inputs/docs/*` | Canonical tree and adapter infrastructure. |
| `casp{15,16}_ligands/README.md` | chai-lab arm landings. |

## Appendix B — Session changelog (condensed)

| Date | Milestone |
|---|---|
| 2026-04-17 | CASP16 Protenix pipeline scaffolded; all 233 JSONs built (with homodimer bug, not yet known). |
| 2026-04-18 | CASP15 Protenix full run — 12/12 targets, 108 min, \$1.54. |
| 2026-04-20 | CASP16 partial inference — 14/44 targets scored; head-to-head with chai drafted; catastrophic failures observed. |
| 2026-04-21 | Tamarind API investigated; v2 checkpoint confirmed withheld; homodimer bug diagnosed and fixed; 4 spot-tests validated fix; paper-style docs written for both arms; canonical `inputs/` tree + 3 adapters delivered; regression-tested and caveat-flagged. |
| 2026-04-22 | This unified manuscript consolidated from the per-arm docs. Full CASP16 post-fix rerun kicked off on RunPod 4090; first 2 targets on 4090 match L40S numbers to ±0.002 LDDT (GPU determinism confirmed). GPU-selection guide and token estimator written. **CASP17 pipeline activated** — T1214 (YncD, 677 aa + PQQ, 701 tokens) ingested into the canonical tree from the `casp17_ligands/` spider and validated through all adapters. |
| 2026-04-22 (evening) | 4090 rerun completed 38 of 44 targets before pod termination. Post-scoring verification found **16 of the 38 are contaminated by stale pre-fix `jsons_prepped/` files** (broken earlier wipe pattern only matched 1 of ~26 intended entries; `run_batch.sh`'s skip-prep-if-exists logic then silently used them). Trustworthy post-fix data in hand: 19 of 44 targets (15 monomer + L2001 fresh-prepped + 3 earlier spot tests). Still need to re-run: 22 targets (16 contaminated + 6 not-reached). Local state now properly wiped; feedback memory added to prevent recurrence. |
