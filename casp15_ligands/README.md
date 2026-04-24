# CASP15 Ligand-Binding Benchmarks

This directory is the **chai-lab** arm of a two-arm benchmark of
AlphaFold3-style protein-ligand structure predictors against the CASP15
ligand track. The companion **Protenix** arm is at
`../casp15_ligands_protenix/`.

Both arms share ground-truth references (`refs/<T>/receptor.pdb` + SDFs,
produced by `prep_references.py`) so that scoring comparisons are
apples-to-apples.

## What's here

```
casp15_ligands/
├── README.md                     # this file — the master landing
├── CLAUDE.md                     # pipeline cheat-sheet (chai side)
│
│ — CHAI pipeline scripts —
├── build_chai_fasta.py           # PDB → chai-lab FASTA
├── prefetch_msas.py              # ColabFold MSA prefetch for chai
├── fetch_chai_downloads.sh       # cache chai model weights locally
├── prep_references.py            # build refs/<T>/ from _lig.pdb
├── score_lddt_pli.py             # chai-style OST LDDT_pli scorer
├── summarize_results.py          # per-target score/confidence table
├── sanity_check.py               # output completeness + diversity
├── compare_to_pdb.py             # per-target Kabsch RMSD + footprint overlap
├── run_batch.sh                  # pod-side chai fold runner
│
│ — manifests —
├── manifest_A40.txt              # 7 targets ≤1024 tok (A40/L40S fits)
├── manifest_A100.txt             # 5 targets 1536-2048 tok (A100/H100)
│
│ — source data —
├── targets_ligand/               # CASP15 bundle (unpacked)
├── fastas/                       # chai-style input FASTAs (1 per target)
├── pdb/                          # companion PDB entries
├── refs/                         # prepared refs (SHARED with Protenix arm)
│
│ — outputs —
├── results_A40/                  # 7-target chai folds on A40
├── results_H100/                 # large-target chai folds on H100
├── msas/                         # precomputed MSAs
├── A40_summary.md                # snapshot of latest A40 run
│
│ — analysis writeups (curated) —
├── binding_sites.md              # per-target binding-site notes
├── binding_descriptions.md       # mechanistic descriptions
├── pdb_ligand_analysis.md        # PDB-level ligand inventory
│
│ — metadata —
├── ligand_targets.json           # target → CCD codes
├── ligand_targets_mapping.csv    # target → PDB → CCD + protein metadata
└── ligand_targets_with_dois.{csv,json}  # + DOI crosslinks to literature
```

## Two pipelines, same benchmark

| | chai-lab | Protenix-v1.0.0 |
|---|---|---|
| **Directory** | here | `../casp15_ligands_protenix/` |
| **Model** | chai-1 0.6.x | Protenix v1.0.0 (368M) |
| **Input** | FASTA (`>protein` + `>ligand|...SMILES`) | JSON (`proteinChain` + `ligand CCD_X`) |
| **MSA** | ColabFold (`api.colabfold.com`) | Protenix MSA server (`protenix-server.com/api/msa`) |
| **Templates** | no (in our config) | yes |
| **Samples/target** | 5 models | 25 (5 seeds × 5 samples) |
| **Scoring** | `score_lddt_pli.py` → OST | shared OST scorer |
| **Writeup** | `A40_summary.md`, `binding_sites.md`, `binding_descriptions.md` | `docs/01_methods.md` through `docs/05_reproduction.md` in `_protenix/` |

## Target scope (both arms)

| Manifest | Targets | Token range |
|---|---|---|
| `manifest_A40.txt` | T1124, T1127, T1127v2, T1146, T1152, T1187, T1188 | ≤1024 |
| `manifest_A100.txt` | H1135, T1158v1, T1158v2, T1158v3, T1158v4 | 1536–2048 |
| **Total** | **12** | — |

Excluded from both (and kept consistent between arms):
RNA-only R* targets, T1170, T1181, T1186, H1114, H1171v1–v2,
H1172v1–v4. The source bundle contains these but they were
not folded in either pipeline.

## Quick results snapshot

See `../casp15_ligands_protenix/docs/02_results.md` for the
definitive Protenix-v1.0.0 results. Headline from there:

- 6 / 12 targets near-native (LDDT_pli ≥ 0.9 on primary ligand)
- 3 / 12 mid-range (0.5 – 0.9)
- 2 / 12 catastrophic (LDDT_pli ≤ 0.2)
- 1 / 12 unscored (OST graph-iso failure on XH0)
- Total L40S GPU time for all 12: ~108 min, ~$1.54

For chai results, see `A40_summary.md` (latest snapshot) and the
`out/...` tree (scoring CSV per batch).

## Which pipeline is "better"?

Neither, consistently. On the related CASP16 pharma-ligand head-to-head
(`../casp16_ligands_protenix/` vs `../casp16_ligands/`):

- Protenix wins 4 targets
- chai wins 3 targets
- Tie / both strong or both fail on 7 targets

They have different characteristic failure modes. For pose prediction
today, running both and comparing outputs is a reasonable triage.

## Related documents

- [`../casp15_ligands_protenix/README.md`](../casp15_ligands_protenix/README.md) —
  Protenix-arm landing (paper-style).
- [`../casp15_ligands_protenix/docs/`](../casp15_ligands_protenix/docs/) —
  full methods + results + discussion + limitations + reproduction.
- [`../casp16_ligands/`](../casp16_ligands/) — CASP16 pharma-ligand chai arm.
- [`../casp16_ligands_protenix/`](../casp16_ligands_protenix/) — CASP16 pharma-ligand Protenix arm.
- [`../casp16/results/ligand_scores.csv`](../casp16/results/ligand_scores.csv) —
  official CASP16 community scores (different target subset: D1273, R1261-4, R1288, T1214).

## History

- **2025** — chai-lab runs on CASP15 (pre-existing work, see
  `A40_summary.md` for batch provenance).
- **2026-02** — Protenix-v1.0.0 released; pipeline extended here.
- **2026-04-20/21** — Protenix-v1.0.0 run against CASP15 on RunPod
  L40S; 12/12 folded, 11/12 scored, ~$1.54 total. Writeup in
  `../casp15_ligands_protenix/docs/`.
