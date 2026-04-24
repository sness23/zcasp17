# 01 — Methods

## 1. Target selection and input preparation

### 1.1 Source data

The CASP15 ligand bundle
`targets_ligand/casp15.targets.ligands.ALL_09.18.2025.tar.gz`
(mirrored from [predictioncenter.org](https://predictioncenter.org/casp15/))
unpacks into `targets_ligand/Targets_ligand/<TARGET>_lig.pdb`. Each file
contains all protein / RNA / DNA chains of the target plus the ligand
HETATM records in a single PDB.

Twenty-five targets are distributed in the bundle; the twelve used here
(Section 1.3) are those validated by the chai-lab reference pipeline at
`../casp15_ligands/`.

### 1.2 JSON construction

`build_protenix_json.py` converts each `<TARGET>_lig.pdb` into the
Protenix input JSON format:

1. Read the PDB with `gemmi`.
2. For each chain, emit a one-letter protein sequence (`protein`), an
   RNA base string (`rna`), or DNA base string (`dna`), skipping waters,
   non-standard residues, and TYR HETATMs (a chai-compatibility
   exclusion).
3. Deduplicate chains by sequence: identical chains collapse into a
   single `{"proteinChain": {"sequence": S, "count": N}}` entity with
   `count` ≥ 2. Slightly-different chains (e.g. a homodimer with one
   truncated subunit) become separate entities.
4. Deduplicate ligand HETATMs by three-letter CCD code: multiple copies
   become `{"ligand": {"ligand": "CCD_X", "count": N}}`.
5. Wrap in a single-element list: `[{"name": TARGET, "sequences": [...]}]`.

Ligands are referenced by CCD code (`CCD_SAH`, `CCD_EPE`, etc.) and
resolved to SMILES at inference time via Protenix's cached RCSB CCD
copy. No local SMILES lookup is performed.

### 1.3 Manifests

| Manifest | Targets | Tokens | Note |
|---|---|---|---|
| `manifest_A40.txt` | T1124, T1127, T1127v2, T1146, T1152, T1187, T1188 | ≤1024 | Fits 24 GB |
| `manifest_A100.txt` | H1135, T1158v1–v4 | 1536–2048 | 48 GB comfortable |
| `manifest_all.txt` | All 12 | mixed | — |

RNA-only targets (R\*) are excluded because the chai baseline excluded
them; targets present in the source bundle but absent from both chai
manifests (T1170, T1181, T1186, H1114, H1171v1-v2, H1172v1-v4) are
excluded for symmetry.

## 2. Feature generation (`protenix prep`)

### 2.1 MSA search

For each `proteinChain` sequence:

1. `POST https://protenix-server.com/api/msa/ticket/msa` with the query
   sequence (no pairing).
2. `POST https://protenix-server.com/api/msa/ticket/pair` with the
   concatenated complex sequence (for taxonomic pairing).
3. Poll `GET /ticket/{ID}` every 60 s until status transitions from
   `PENDING` → `RUNNING` → `COMPLETE`.
4. Download and untar the result; rearrange into `pairing.a3m`
   (taxonomy-tagged) and `non_pairing.a3m` (all hits) under
   `jsons_prepped/<T>/msa/0/`.

Queue wait dominated observed wall time: cached sequences returned in
~10 s; uncached sequences took 2–30 min. The server is a ByteDance-
hosted MMseqs2-App instance; its compute is identical to
`api.colabfold.com` but has its own queue.

### 2.2 Template search

Local `hmmer` + `kalign` against the PDB sequence database
(`pdb_seqres_2022_09_28.fasta`, auto-downloaded to
`$PROTENIX_ROOT_DIR/search_database/`):

1. `hmmbuild` an HMM profile from the two MSAs.
2. `hmmsearch` the profile against `pdb_seqres_2022_09_28.fasta` with
   `-E 100 --incE 100`.
3. Realign hits with `kalign` into `hmmsearch.a3m`.

Typical latency: 0.2 s `hmmbuild`, 7 s `hmmsearch` on 8 CPU cores.

### 2.3 Output

For each target, `jsons_prepped/<T>.json` is the original JSON
augmented with three per-chain paths:

```json
"pairedMsaPath":   ".../jsons_prepped/<T>/msa/0/pairing.a3m",
"unpairedMsaPath": ".../jsons_prepped/<T>/msa/0/non_pairing.a3m",
"templatesPath":   ".../jsons_prepped/<T>/msa/0/hmmsearch.a3m"
```

For full details of what each file contains and how Protenix's model
consumes them, see [prep_phase.md](prep_phase.md).

### 2.4 MSA reuse across targets

Two post-prep optimizations reduce queue cost:

- `clone_prepped.py` — after each successful prep, detects other targets
  in the manifest whose protein sequence is identical to the just-prepped
  target and copies the MSA paths into their prepped JSON, skipping the
  redundant MSA query.
- `prep_local.sh` auto-invokes the cloner after each target.

For CASP15, chain-sequence diversity is high (mostly unique sequences
per target) so the cloner saves little. For CASP16 Mpro (25 targets
sharing the same receptor family with 1–7 residue truncation
differences), reuse is more impactful.

## 3. Inference (`protenix pred`)

### 3.1 Model

`protenix_base_default_v1.0.0` — 368M parameters, AlphaFold3-equivalent
training data cutoff of 2021-09-30. Supports MSA + RNA MSA + template
features. Downloaded automatically from ByteDance TOS to
`$PROTENIX_ROOT_DIR/checkpoint/` on first inference call (~1 GB).

The nominally-current `protenix-v2` (464M params) is referenced in the
Protenix repo but the checkpoint URL returns HTTP 403 (see
[ByteDance/Protenix#295](https://github.com/bytedance/Protenix/issues/295));
weights were never publicly released.

### 3.2 Sampling

- **5 seeds**: 101, 102, 103, 104, 105.
- **5 samples per seed** via diffusion at N_step = 200.
- **10 recycles** (`model.N_cycle = 10`).
- **Templates enabled** (`--use_template true`).
- **MSAs enabled** by default.
- **Mixed precision**: bf16 throughout, with `SampleDiffusion` and
  `ConfidenceHead` pinned to fp32 for numerical stability (Protenix's
  `update_inference_configs` helper lowers these to bf16 only at
  N_token > 2560, which no CASP15 target reaches).

Total: 25 CIF structures per target × 12 targets = **300 structures**.

### 3.3 Kernels

`--triangle_attention cuequivariance --triangle_multiplicative cuequivariance`
(NVIDIA cuequivariance kernels) via Protenix's defaults. The custom
`fast_layer_norm` kernel JIT-compiles on first use per GPU arch (cached
afterward).

### 3.4 Hardware and wall times

Runs on RunPod L40S (48 GB VRAM, PyTorch 2.5.1 + CUDA 12.4 image):

| Manifest | Targets | Wall time | Per-target avg |
|---|---|---|---|
| A40 | 7 | 31 min | 2.5–5 min (first target 7.8 min inc. warm-up) |
| A100 | 5 | 77 min | 12–27 min (H1135 biggest at 27 min) |
| All 12 | 12 | **~108 min** | — |

At L40S's $0.86/hr, total GPU cost ~$1.54.

Spot-tested L2001 (CASP16) on Vast.ai RTX 3090 ($0.18/hr): 578 s
including prep queue. Validates 24 GB VRAM is sufficient for most
CASP15 and CASP16 targets.

## 4. Scoring

### 4.1 Pose selection

For each target, `score_lddt_pli.py --best-only` selects one
"best" sample by highest `ranking_score` from
`<T>_summary_confidence_sample_N.json`, across all 25 samples (5 seeds ×
5 samples). `ranking_score` is Protenix's scalar confidence combining
pLDDT, pTM, ipTM, and clash penalties (see [02_results.md §4](02_results.md#4-confidence-calibration)).

### 4.2 OST comparison

OpenStructure 2.x `compare-ligand-structures` via the `ost` conda env:

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

Two metrics are reported:

- **LDDT_pli** — local Distance Difference Test restricted to
  protein-ligand interaction atoms. Range 0–1; higher is better. Counts
  contacts within a 6 Å radius.
- **BiSyRMSD** — symmetry-aware (permutation-accounting) RMSD of the
  aligned ligand atoms to the reference, in Ångström. Lower is better.

`--substructure-match` tolerates partial ligand resolution and
tautomer-level differences during graph isomorphism matching between
the model CIF's ligand and the reference SDF.

### 4.3 Reference preparation

`refs/<T>/receptor.pdb` and `refs/<T>/lig_*.sdf` are produced from
`<T>_lig.pdb` by `../casp15_ligands/prep_references.py`:

- Protein chains saved as `receptor.pdb`.
- Each HETATM ligand saved as a separate SDF with 3-letter CCD code
  preserved in the filename (`lig_NN_CCD.sdf`).

These refs are shared with the chai baseline via the
`refs -> ../casp15_ligands/refs/` symlink; no Protenix-specific ref
generation is performed.

## 5. Infrastructure

- **Local workstation:** Linux, 64 GB RAM, conda env `protenix`
  (Python 3.11, `protenix` via pip + `bioconda::hmmer` + apt `kalign`).
  Used for JSON construction, local MSA/template prep, and scoring.
- **Remote GPU:** RunPod L40S 48 GB (primary), Vast.ai RTX 3090 24 GB
  (fallback). Ephemeral pods with PyTorch 2.5.1 + CUDA 12.4 base image;
  protenix installed via `pip install protenix` on each fresh pod.
- **Sync convention:** Same absolute path (`/workspace/casp15_ligands_protenix/`)
  on both local and pod, so MSA paths baked into prepped JSONs resolve
  identically. Explicit rsync between local and pod — no live mount.
