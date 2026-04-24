# 03 — Discussion

## 1. Where Protenix-v1.0.0 wins

**Protein-ligand targets with prominent binding sites, abundant
structural homologs in the PDB, and ligands with moderate-to-large
chemical graphs track the model's strengths.** T1124 (SAH bound to a
methyltransferase scaffold), T1127 / T1127v2 (HEPES in a nucleotide-
binding fold), and T1158v4 (ATP + Mg²⁺ in a kinase-like architecture)
all scored LDDT_pli ≥ 0.93 with RMSD ≤ 0.6 Å on the primary ligand. In
each case:

- The receptor fold has abundant templates available in
  `pdb_seqres_2022-09-28.fasta`.
- The bound ligand is a common cofactor well-represented in the training
  set.
- The binding pocket is the single thermodynamically dominant site on
  the protein.

H1135 (9 ordered K⁺ ions, all LDDT_pli = 1.00) is a special case:
a potassium channel or transporter with a tight cation-binding cage.
Protenix's LDDT_pli on ions is conservatively high because the contact
count is small (lddt_pli_n_contacts ≈ 30–60) — but the RMSDs
(0.08–0.22 Å) are independently rock-solid.

## 2. Where it breaks down

### 2.1 Small or weakly-bound solvent additives

MPD (2-methyl-2,4-pentanediol) on T1127 / T1127v2 is misplaced by ~22 Å
— model places it in the wrong pocket, possibly because MPD in the
training set is often a crystallographic artifact and the model has
learned not to treat it as a meaningful binding event.

Halide ions (Cl⁻) on H1135 show the same pattern: 10 Å errors despite
the same chain being correctly positioned for 9 K⁺ ions. Cl⁻ is weakly
coordinated (no well-defined pocket in most deposited structures),
and Protenix's energy landscape for halides appears shallow.

### 2.2 Heterogeneous ligand sets

T1188 has five ligands of four different chemical classes (DW0 drug
fragment, Cd²⁺ ×2, Co²⁺ ×1, CO). Protenix placed DW0 within 2.5 Å but
two of the four metals ended up >30 Å away — likely in symmetry-related
sites on the homodimer that don't match the ref crystallographic
assignment. OST's symmetry-aware RMSD (BiSyRMSD) should tolerate this
by permutation, but when the model places ions at genuinely non-native
sites, no permutation helps.

### 2.3 Glycan variability

NAG is present on three targets with three different outcomes:

- T1146: LDDT_pli 0.84 — clean, single glycosylation site, well-anchored.
- T1152: LDDT_pli 0.52 — partial; the NAG is in approximately the right
  region but the ring orientation is off by ~90°.
- T1187: LDDT_pli 0.15 — genuinely failed, ~6 Å from the crystal
  position.

This bimodality is consistent with general observations that glycan
placement depends strongly on whether the glycosylation site is
well-conserved and templated vs exploratory.

### 2.4 OST graph-isomorphism failures

Three targets in the broader CASP15+CASP16 PoC set (T1158v3 here;
L1006, L1007, L1008 in CASP16) score 0 ligand pairs not because the
pose is wrong, but because
`openstructure.compare-ligand-structures` cannot match the model
ligand's atomic graph to the reference SDF's atomic graph via graph
isomorphism (even with `--substructure-match` engaged).

Root cause: Protenix emits ligand CIF records with an internal bond-
perception scheme that occasionally differs from the RDKit-canonical
graph of the reference SDF, typically in:

- Tautomer assignment (amide vs imidol form).
- Bond-order disambiguation across aromatic systems.
- Charge placement on zwitterionic or metal-coordinating ligands.

The same three CASP16 targets failed under chai-lab scoring with the
identical error. So this is a **scorer-side tooling issue, not a
model failure**, and is independent of the structure-predictor choice.

Mitigations not yet applied:

1. Force-canonicalize the reference SDF via RDKit MolFromSmiles →
   MolToSmiles → MolFromSmiles roundtrip before scoring.
2. Switch to a scoring tool that ignores bond orders (e.g. PoseBusters,
   LIGPLOT-style contact scoring).
3. Write a narrower coordinate-only RMSD fallback when OST returns 0
   pairs.

None affect the primary conclusion; they'd just convert three `—`
entries in the results table into numbers.

## 3. Confidence calibration

Protenix's `ranking_score` is **target-level, not ligand-level**. Three
distinct failure modes illustrate:

| Case | ranking_score | Actual outcome |
|---|---|---|
| **Trustworthy high-conf** | ≥ 0.95 on T1124, T1127, T1146, T1188 | 5/5 targets had a strong primary ligand result |
| **Blurred confidence** | 0.81 on T1158v4 | Great primary ligands (ATP, Mg²⁺) — low score reflects overall structure uncertainty, not ligand pose |
| **Correctly-low confidence** | 0.40 on T1187 | Poor ligand pose (0.15 LDDT_pli) — model is honest that it doesn't know |
| **Misleading high conf** | 0.94 on T1188 | DW0 OK but Cd²⁺ and Co²⁺ placed ~30 Å away; score diluted by protein correctness |

**Operational recommendation:** For multi-ligand complexes, don't trust
`ranking_score` alone. Use `chain_plddt` from the summary-confidence
JSON to get per-entity confidence, and cross-check the ligand's local
pLDDT against the receptor's.

## 4. Comparison to chai-lab baseline

The CASP15 chai numbers live at `../casp15_ligands/out/...` but were
produced with a slightly different scoring config (chai's
`score_lddt_pli.py` without `--substructure-match`). Rather than
construct a possibly-stale head-to-head here, see the head-to-head
on the CASP16 pharma set in `../casp16_ligands_protenix/` where
both methods were run through this same scorer under identical
conditions.

The headline from that CASP16 comparison (14 shared targets):

- Protenix better on 4 targets
- Chai better on 3 targets
- Tie (both strong or both weak) on 4
- Both fail graph-iso on 3 (scoring issue, not method)

Protenix's wins on CASP16 L4002 (chai's 34 Å → Protenix's 0.45 Å) and
the corresponding chai win on L4001 show that **neither method
dominates**: they have different failure modes on the same hard
targets. For a CASP-like "which single tool to run" decision, current
evidence says run both and compare outputs.

## 5. Pipeline lessons

### 5.1 Homodimer representation matters

Earlier drafts of the CASP16 pipeline concatenated homodimer chains
into a single long sequence. This produced catastrophic failures
(L2001 LDDT_pli 0.00 / 50 Å RMSD) because the model couldn't place the
ligand at the correct inter-chain interface.

Switching to Protenix's native `{"proteinChain": {..., "count": 2}}`
representation (for identical chains) or separate entries for
near-identical chains (truncation variants) fixed L2001 dramatically
(LDDT_pli 0.75, RMSD 1.91 Å).

CASP15's builder was correct from the start — it deduplicates chains
by sequence — so the same targets did not suffer the same regression
here.

### 5.2 MSA queue is the wall-clock bottleneck

On the happy path, GPU inference is ~5 min per target. MSA queue wait
at `protenix-server.com/api/msa` can dwarf that: 2–30 min per uncached
sequence. For the CASP15 run, the prep stage took ~4 hours for 44
CASP16 targets in a cold cache, vs ~100 minutes of GPU time for all
12 CASP15 structures.

Mitigation: run `prep_local.sh` on a cheap local CPU box (or a small
always-on instance) and rsync only the prepped JSONs to GPU pods.
This decouples the free-but-slow MSA server from the paid GPU clock.

### 5.3 Workspace-path consistency

Both our local workstation and the rented GPU pods use the same
absolute path `/workspace/casp15_ligands_protenix/`. This lets the
absolute paths baked into the prepped JSONs
(`/workspace/.../pairing.a3m`) resolve identically on both sides,
eliminating the need for a sed-based path-rewrite step during each
push-pull cycle.
