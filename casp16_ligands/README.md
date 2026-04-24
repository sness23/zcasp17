# CASP16 Pharma-Ligand Benchmarks (chai-lab arm)

This directory is the **chai-lab** arm of a two-arm benchmark of
AlphaFold3-style protein-ligand structure predictors against the CASP16
pharma-ligand track. The companion **Protenix** arm is at
`../casp16_ligands_protenix/`.

Both arms share ground-truth references (`refs/<T>/receptor.pdb` + SDFs,
produced by `prep_references.py`) so that scoring comparisons are
apples-to-apples.

> ⚠️ **Known issue — homodimer input bug in `fastas/`.**
>
> The FASTA files in `fastas/` were built by `build_chai_fasta.py` with
> a bug: it concatenates homodimer chains from `protein_aligned.pdb`
> into a single long "fake-monomer" sequence. This affects all 27
> targets with multi-chain receptors (all 25 L4000 Mpro targets + both
> L2000 Cathepsin G targets; the 17 L1000 Chymase targets are monomer
> and unaffected).
>
> The Protenix side hit the identical bug and fixed it
> ([see `../casp16_ligands_protenix/docs/03_discussion.md §1`](../casp16_ligands_protenix/docs/03_discussion.md#1-the-homodimer-bug)).
> After the fix, Protenix's L2001 went from LDDT_pli 0.00 to 0.75,
> and L4001 from 0.01 to 0.78 — in other words, the catastrophic
> multi-chain failures were input-driven, not method-driven.
>
> **chai's results in this directory were generated with the buggy
> input** and have not yet been re-run with corrected FASTAs. chai's
> reported catastrophic failures on L2001, L4003 (40+ Å RMSD) may also
> be input-related, not method failures.
>
> **Fix:** regenerate FASTAs from the canonical `inputs/` tree:
> ```bash
> python3 ../inputs/adapters/to_chai_fasta.py \
>     ../inputs/casp16/*/target.json \
>     --out-dir fastas/
> ```
> The canonical adapter correctly emits multi-chain receptors as
> multiple `>protein|...` records, matching the real dimer structure.
>
> Until the chai pipeline is re-run with corrected FASTAs, **the
> head-to-head comparison in
> `../casp16_ligands_protenix/docs/02_results.md §2` understates
> chai's actual performance** on multi-chain targets.

## The pharma-ligand track

CASP16 introduced four drug-discovery benchmark sets, donated by
pharmaceutical partners. Each set has one crystal-characterised receptor
plus many candidate compounds at known binding poses:

| Set | Receptor | Source PDB | # targets | Notes |
|---|---|---|---|---|
| L1000 | Chymase | 1T31 | 17 | Serine protease, monomer |
| L2000 | Cathepsin G | 1CGH | 2 | Cysteine protease, homodimer |
| L3000 | Autotaxin | 5M7M | 219 | Lipid signaling enzyme, largest set |
| L4000 | SARS-CoV-2 Mpro | 7L11 | 25 | Coronavirus main protease, homodimer |

A "PoC" subset (`manifest_L1L2L4.txt`) covers 44 targets across L1/L2/L4
— used for the published head-to-head benchmarks here and in the
Protenix arm. L3000 is deferred (219 targets → ~$47 L40S time).

## Pipeline

Chai-lab takes a FASTA with:

```
>protein|name=<TARGET>-receptor
<single-chain protein sequence>
>ligand|name=<TARGET>-lig-<Name>
<SMILES>
```

And produces 5 model CIF files per target. Our pipeline:

```
  protein_aligned.pdb + <TARGET>.tsv (SMILES)
          │
          ▼  build_chai_fasta.py
      fastas/<T>.fasta
          │
          ▼  run_batch.sh (on GPU pod)
      results_A40/<T>/pred.model_idx_{0..4}.cif
                   + scores.model_idx_*.npz
          │
          └─▶ score_lddt_pli.py → out/lddt_pli/summary.csv
                                  (OST compare-ligand-structures)
```

## Directory contents

```
casp16_ligands/
├── README.md                     # this file — landing for the chai arm
│
│ — pipeline scripts —
├── build_chai_fasta.py           # PDB + SMILES → chai FASTA
├── prep_references.py            # crystal ligand PDBs → SDF + receptor.pdb
├── run_batch.sh                  # chai fold runner (pod-side)
├── score_lddt_pli.py             # OST LDDT_pli + BiSyRMSD scorer
├── summarize_results.py          # per-target aggregate/pTM/ipTM table
├── sanity_check.py               # output completeness + inter-model RMSD
│
│ — manifests —
├── manifest_L1000.txt            # 17 Chymase targets
├── manifest_L2000.txt            #  2 Cathepsin G targets
├── manifest_L3000.txt            # 219 Autotaxin targets (deferred)
├── manifest_L4000.txt            # 25 Mpro targets
├── manifest_L1L2L4.txt           # PoC 44-target subset
│
│ — per-set sources (SMILES TSVs) —
├── L1000/                        # <TARGET>.tsv — one per Chymase target
├── L2000/                        # ditto Cathepsin G
├── L3000/                        # ditto Autotaxin
├── L4000/                        # ditto Mpro
│
│ — built inputs —
├── fastas/                       # generated chai FASTAs (~44 files)
├── refs/                         # prepared ground truth (SHARED via symlink from Protenix arm)
│
│ — outputs —
├── out/lddt_pli_L40S/summary.csv # chai scoring output, the head-to-head reference
└── ...
```

## Results summary (chai-lab on 14 scored of 44 PoC targets)

Full table in [`out/lddt_pli_L40S/summary.csv`](out/lddt_pli_L40S/summary.csv).
Highlights relative to the Protenix arm (from
[`../casp16_ligands_protenix/docs/02_results.md`](../casp16_ligands_protenix/docs/02_results.md#2-head-to-head-vs-chai-lab)):

| Outcome | Count | Targets |
|---|---|---|
| chai wins (LDDT_pli Δ ≥ 0.05) | 3 | L1001, L1002, L4001 |
| Protenix wins | 4 | L1003, L1004, L4002, L4004 |
| Tie / both strong | 2 | L1005, L2002 |
| Both catastrophic (>40 Å) | 2 | L2001, L4003 |
| Both OST-unscored (graph iso) | 3 | L1006, L1007, L1008 |

Neither method dominates. Most striking cases:

- **L4002:** chai placed the ligand 34 Å from the pocket; Protenix hit
  0.45 Å. Large Protenix win.
- **L4001:** chai 2.14 Å; Protenix 47 Å pre-fix, now 1.78 Å post-fix
  (see Protenix arm for the homodimer bug story).

## Running the chai pipeline

From a fresh pod with a 24+ GB GPU:

```bash
cd /workspace/casp16_ligands
pip install chai_lab
apt-get install -y rsync
export CHAI_DOWNLOADS_DIR=/workspace/chai_downloads
mkdir -p $CHAI_DOWNLOADS_DIR

./run_batch.sh manifest_L1L2L4.txt
# → results_A40/<TARGET>/pred.model_idx_*.cif
```

For details see the per-script docstrings and the
[original chai-lab docs](https://github.com/chaidiscovery/chai-lab).

## Related documents

- [`../casp16_ligands_protenix/`](../casp16_ligands_protenix/) — Protenix
  v1.0.0 arm with paper-style docs in `docs/01_methods.md`
  through `docs/05_reproduction.md`.
- [`../casp15_ligands/`](../casp15_ligands/) — CASP15 ligand benchmarks
  (chai arm) and `../casp15_ligands_protenix/` for the Protenix arm.
- [`../casp16/results/ligand_scores.csv`](../casp16/results/ligand_scores.csv)
  — official CASP16 community scores for the main results targets
  (D1273, R1261-4, R1288, T1214 — not the pharma track).
- Curated prose: [`binding_sites.md`](binding_sites.md),
  [`binding_descriptions.md`](binding_descriptions.md).
