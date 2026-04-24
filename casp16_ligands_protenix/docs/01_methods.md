# 01 — Methods

## 1. Target selection

### 1.1 CASP16 pharma ligand track

CASP16 introduced a new "pharma ligand" track: four drug-discovery
benchmark sets donated by pharma companies, each with a single crystal-
characterized receptor and many candidate compounds at known binding
poses. The four sets:

| Set | Receptor | # targets | Source PDB | Protein class |
|---|---|---|---|---|
| L1000 | Chymase | 17 | 1T31 | Serine protease |
| L2000 | Cathepsin G | 2 | 1CGH | Cysteine protease |
| L3000 | Autotaxin | 219 | 5M7M | Lipid signaling enzyme |
| L4000 | SARS-CoV-2 Mpro | 25 | 7L11 | Main protease |

Source layout from the CASP16 spider (`../casp16/pharma_ligands/`):

```
exper_struct/<SET>/<SET>_prepared/<TARGET>/
    protein_aligned.pdb         # crystal receptor (pre-aligned across the set)
    ligand_<resid>_<chain>_<n>.pdb   # crystal ligand pose(s)
```

Per-target ligand SMILES are in parallel TSVs at
`../casp16_ligands/<SET>/<TARGET>.tsv` (`ID\tName\tSMILES\tTask`).

### 1.2 PoC subset

Full coverage would be 263 targets. For the PoC reported here we used
`manifest_L1L2L4.txt` (44 targets: L1000 ×17 + L2000 ×2 + L4000 ×25),
excluding L3000's 219 Autotaxin targets. L3000 is deferred because:

- ~15 min/target × 219 = ~55 GPU hours ≈ $47 of L40S time — a big PoC
  without a clear hypothesis to test.
- Confirming pipeline correctness (homodimer fix) and getting preliminary
  head-to-head numbers vs chai was the PoC's primary goal.

## 2. Input construction

### 2.1 JSON construction

`build_protenix_json.py` reads each target's receptor PDB and ligand
TSV, emits a Protenix input JSON:

```json
[{
  "name": "L4001",
  "sequences": [
    {"proteinChain": {"sequence": <chain A seq>, "count": 1}},
    {"proteinChain": {"sequence": <chain B seq>, "count": 1}},
    {"ligand": {"ligand": <SMILES>, "count": 1}}
  ]
}]
```

Per-chain sequence extraction is delegated to `gemmi`. Chains with
identical sequences collapse into `{"proteinChain": {..., "count": N}}`;
chains differing by ≥1 residue (typical in CASP16 — crystal
reconstructions vary by truncation) become separate entities.

### 2.2 The homodimer bug (and fix)

**Earlier versions concatenated all chains into one long sequence**
(e.g. L2001 Cathepsin G: chain A 226 aa + chain B 227 aa =
453-aa chimera). This produced catastrophic ligand misplacement on
multi-chain targets because Protenix saw a nonexistent fake protein
instead of a dimer.

The current (fixed) builder iterates per-chain and dedupes by exact
sequence identity. Affected targets:

| Set | Multi-chain count | Note |
|---|---|---|
| L1000 (Chymase) | 0 / 17 | All monomer — no bug impact |
| L2000 (Cathepsin G) | 2 / 2 | Homodimer, chains differ by 1 aa |
| L4000 (Mpro) | 25 / 25 | All homodimer, chains differ by 1–7 aa |

Post-fix targets produce 2 separate `proteinChain` entries for most
Mpro / Cathepsin complexes; L4009 (perfectly identical chains) becomes
`{"proteinChain": {..., "count": 2}}`.

Validation spot-tests (see [02_results.md §3](02_results.md#3-homodimer-fix-validation))
showed the fix rescues catastrophic failures (RMSD ≥ 40 Å) to
near-native or mid-range outcomes on L2001, L4001, L4003.

### 2.3 Ligand format

Unlike CASP15 (where ligands are HETATM records with canonical CCD
codes), CASP16 pharma targets use **novel pharma SMILES** with no CCD
assignments. Our builder passes the SMILES directly as the `ligand`
field; Protenix's internal CCD cache is bypassed.

This has two practical consequences:

1. No CCD lookup cache hits — every ligand is built from scratch.
2. The OST graph-isomorphism scorer sometimes fails to match model ←→ ref
   (see [03_discussion.md §2.4](03_discussion.md#24-ost-graph-isomorphism-failures)).

## 3. Feature generation

Identical to CASP15 — see [`../../casp15_ligands_protenix/docs/01_methods.md`](../../casp15_ligands_protenix/docs/01_methods.md#2-feature-generation-protenix-prep)
for details. Summary:

- MSA: remote `protenix-server.com/api/msa` (ColabFold-compatible MMseqs2
  API behind a ticket queue). 2 requests per target (`/ticket/msa`,
  `/ticket/pair`).
- Template: local `hmmbuild → hmmsearch → kalign` against
  `pdb_seqres_2022_09_28.fasta` (~300 MB, auto-downloaded).

CASP16 queue patterns: within-set targets (e.g. all Mpro targets) have
17 unique protein sequences for 25 structures due to truncation
variants. `clone_prepped.py` captures partial cache reuse but most
Mpro targets need fresh MSAs.

Total local prep wall time for the full 44 PoC: ~4 hours (dominated
by MSA server queue).

## 4. Inference

### 4.1 Model + sampling

- **Model:** `protenix_base_default_v1.0.0` (368M). The ByteDance
  v2 (464M) checkpoint is withheld — see
  [ByteDance/Protenix#295](https://github.com/bytedance/Protenix/issues/295).
- **Seeds:** 101, 102, 103, 104, 105 (5 seeds).
- **Samples per seed:** 5 (via diffusion, N_step=200).
- **Total structures per target:** 25.
- **Recycles:** 10 (`model.N_cycle = 10`).
- **Templates:** enabled.
- **dtype:** bf16 mixed precision (SampleDiffusion + ConfidenceHead in
  fp32 for numerical stability).

### 4.2 Hardware

All runs on RunPod and Vast.ai cloud GPUs:

| Run | GPU | Role | Wall | Cost |
|---|---|---|---|---|
| L1L2L4 partial (14 of 44 targets, Ctrl-C'd) | L40S 48GB | PoC initial | ~65 min | $0.93 |
| L2001 post-fix spot test | RTX 3090 24GB (Vast) | Bug confirmation | 9.5 min | $0.03 |
| L4001/L4003/L4004 post-fix spot test | RTX 4090 24GB (RunPod) | Bug confirmation | 33 min | $0.40 |

Ran with `--triangle_multiplicative cuequivariance --triangle_attention cuequivariance`.
Fast LayerNorm CUDA kernel JIT-compiles on first target (~1 min,
cached). Checkpoint download (~1 GB) happens on first target then
cached to `$PROTENIX_ROOT_DIR/checkpoint/`.

## 5. Scoring

### 5.1 Ground truth

Reference PDB + SDF per target at `refs/<TARGET>/`, built from
`ligand_<resid>_<chain>_<n>.pdb` files via `../casp16_ligands/prep_references.py`:

- `refs/<T>/receptor.pdb` — protein chains only.
- `refs/<T>/lig_<NN>_ligand_<ID>_<chain>_<copy>.sdf` — each crystal
  ligand as a separate SDF with standardized naming.

Shared with the chai pipeline via symlink `refs/ → ../casp16_ligands/refs/`.

### 5.2 OST command

Identical to CASP15:

```
ost compare-ligand-structures \
    -m <MODEL_CIF> -mf cif \
    -r refs/<T>/receptor.pdb -rf pdb \
    -rl refs/<T>/lig_*.sdf \
    -of json -ft \
    --lddt-pli --rmsd \
    --substructure-match \
    -o <OUT_JSON>
```

Metrics: LDDT_pli (0–1, higher better) and BiSyRMSD (Å, lower better).

### 5.3 Best-sample selection

`score_lddt_pli.py --best-only` picks one sample per target by highest
`ranking_score` (Protenix's confidence head output). Per-target
scoring output in `out/lddt_pli/<T>/seed<S>_sample<N>.json` +
aggregate CSV at `out/lddt_pli/summary.csv`.

## 6. Head-to-head baseline

Chai-lab 0.6.x ran the same 44 targets at `../casp16_ligands/results_A40/`
on an A40 48GB pod. Its outputs live at
`../casp16_ligands/out/lddt_pli_L40S/summary.csv`. The comparison is a
direct cross-reference of the two CSVs — same refs, same scorer, same
24-core OST build, different predictor.

For the head-to-head tables in [02_results.md](02_results.md#2-head-to-head-vs-chai-lab),
we matched rows by `target` and picked each method's best sample per
target. 11 targets were scored by both methods; 3 (L1006-L1008) failed
OST graph-iso matching in both.
