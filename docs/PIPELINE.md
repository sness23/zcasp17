# Pipeline architecture

End-to-end design of the zcasp17 structure-prediction and scoring pipeline.

## Design goals

1. **Per-target reproducibility** — every step from FASTA/JSON build → fold → score is committed to disk. Re-running gives the same answer.
2. **Resumable** — partial batches can be stopped and restarted; the runner skips targets with completed outputs.
3. **Apples-to-apples with CASP** — the scoring metric is CASP's own `ost compare-ligand-structures`; the comparison CSV is CASP's public ranking data.
4. **Cheap to iterate** — $1–15 total for a CASP-scale benchmark, assuming a $0.40–$0.90/hr GPU.
5. **Model-agnostic** — swapping Chai-1 for Protenix for Boltz is one adapter script.

## Data flow

```
                              ┌────────────────────┐
   CASP-provided              │  <TARGET>_lig.pdb  │
   ground truth               │  protein_aligned   │
                              │    + ligand_*.pdb  │
                              └─────────┬──────────┘
                                        │
                                        ▼
                              ┌────────────────────┐
                              │  prep_references   │
                              │  (split crystal    │
                              │   → receptor.pdb   │
                              │     + lig_*.sdf    │
                              │   per target)      │
                              └─────────┬──────────┘
                                        ▼
                              ┌────────────────────┐
                              │  refs/<TARGET>/    │
                              │   receptor.pdb     │
                              │   lig_*.sdf        │
                              └────────────────────┘

                              ┌────────────────────┐
                              │  CANONICAL INPUTS  │
                              │  inputs/<casp>/    │
                              │    <T>/target.json │
                              └─────────┬──────────┘
                                        │
                        ┌───────────────┴───────────────┐
                        ▼                               ▼
              ┌────────────────────┐          ┌────────────────────┐
              │ to_protenix_json   │          │ to_chai_fasta      │
              └─────────┬──────────┘          └─────────┬──────────┘
                        │                               │
                        ▼                               ▼
              ┌────────────────────┐          ┌────────────────────┐
              │ Protenix prep      │          │ Chai FASTA files   │
              │ + MSA + templates  │          │                    │
              └─────────┬──────────┘          └─────────┬──────────┘
                        │                               │
                        ▼                               ▼
              ┌────────────────────┐          ┌────────────────────┐
              │ Protenix inference │          │ chai-lab fold      │
              │ GPU pod            │          │ GPU pod            │
              │ 5 seeds × 5 samples│          │ 5 samples          │
              └─────────┬──────────┘          └─────────┬──────────┘
                        └───────────────┬───────────────┘
                                        ▼
                              ┌────────────────────┐
                              │ results/<TARGET>/  │
                              │  pred.*.cif        │
                              │  scores.json       │
                              └─────────┬──────────┘
                                        │
                                        ▼
                              ┌────────────────────┐
                              │  score_lddt_pli    │
                              │  (OST compare-lig  │
                              │   structures)      │
                              └─────────┬──────────┘
                                        ▼
                              ┌────────────────────┐
                              │  out/lddt_pli/     │
                              │   summary.csv      │
                              └─────────┬──────────┘
                                        ▼
                              ┌────────────────────┐
                              │  CASP rankings CSV │
                              │  (ligand_scores)   │
                              │  join + rank       │
                              └────────────────────┘
```

## Components

| Component | Responsibility | Implementation |
|---|---|---|
| **Canonical inputs** | Per-target receptor sequences, ligand SMILES, crystal references | `inputs/<casp>/<T>/target.json` |
| **Adapters** | Canonical → pipeline-specific inputs | `inputs/adapters/to_{protenix,chai,dynamicbind}.py` |
| **SMILES cache** | Idempotent local store of CCD → isomeric SMILES | `.smiles_cache.json` |
| **Manifest** | List of targets to fold in a batch | `manifest_<bucket>.txt` |
| **Batch runner** | Fold every manifest entry, resumable, log-captured | `run_batch.sh` |
| **Weight cache** | Pre-populate predictor weights so every pod starts instantly | `~/data/chai_downloads/`, `~/data/protenix_downloads/` |
| **Reference prep** | Split CASP ground truth → receptor.pdb + ligand SDF | `prep_references.py` |
| **Sanity checker** | Output completeness; pLDDT; inter-model RMSD | `sanity_check.py` |
| **Protein comparator** | Fast CA RMSD vs crystal (per-chain + overall) | `compare_to_pdb.py` |
| **CASP-metric scorer** | `ost compare-ligand-structures` wrapper → per-ligand LDDT_pli + BiSyRMSD | `score_lddt_pli.py` |
| **CASP comparator** | Join our `summary.csv` to CASP's `ligand_scores.csv`; compute rank | pandas notebooks |

## GPU / token budgeting

Chai-1 pads inputs to `AVAILABLE_MODEL_SIZES = [256, 384, 512, 768, 1024, 1536, 2048]`. Tokens = protein residues + ligand heavy atoms.

| Bucket | Fit | Recommended GPU | Hourly |
|---|---|---|---|
| ≤ 1024 | comfortable | A40 48GB | $0.44 |
| 1536 | tight | L40S 48GB | $0.86 |
| 2048 | 80GB required | A100 / H100 | $2–3 |
| > 2048 | can't fold end-to-end | chunk or skip | — |

Batching strategy: one manifest per GPU bucket. A40 for ≤1024, L40S for 1536, H100 NVL for 2048.

## Why the scoring is modular

Three comparison scripts, each with a specific purpose:

1. **`sanity_check.py`** — output completeness + model diversity. Catches broken outputs, OOMs, mis-parsed ligands, and flags targets where the 5 diffusion samples disagree wildly (a red flag even if the best-agg score is fine).
2. **`compare_to_pdb.py`** — fast proxy vs crystal. CA RMSD + ligand centroid distance. Useful for overview; known to flatter the model on ligand quality.
3. **`score_lddt_pli.py`** — CASP-metric, slow (~1–3 s per target). Requires OpenStructure. This is the number to cite.

Separation of concerns: the first two are cheap and need no external tools, so they're the default "did the batch work" check. The third is slower and needs a separate conda env but gives CASP-comparable numbers.

## Weight management

Every AF3-class predictor ships multi-GB weights. Three observations shape the design:

1. **Always** set `$CHAI_DOWNLOADS_DIR` / `$PROTENIX_CACHE` to a persistent path — not the default (inside site-packages).
2. **ESM-2's 5.7 GB file** is the most common source of pod-failure via `IncompleteRead`. Chai-lab's download is atomic (tmp file + rename), so retries work, but pre-populating is faster.
3. A **vault-side cache** rsync'd to every new pod is the winning strategy. One-time 10–20 min local download; instant subsequent pods.

## Resumability

Two principles:

- **Batch runner**: `run_batch.sh` skips a target if its output dir already has N CIFs. Safe to Ctrl-C and restart.
- **Idempotent builders**: regenerating a FASTA or prepped reference produces the same bytes. The SMILES cache dedupes RCSB lookups.

Combined, you can kill pods mid-batch, terminate early, re-fold individual targets, and the pipeline converges.

## Failure modes we've seen and handled

| Failure | Cause | Fix |
|---|---|---|
| `IncompleteRead` on ESM-2 download | Pod network flakiness | Pre-populated weight cache |
| `CUDA out of memory` at 1536 tokens on A40 | Wrong bucket assignment | Manifest partitioning by token count |
| LDDT_pli NaN on single-atom ions | OST graph-iso can't match single atoms | Accept as scoring limitation; flag in summary |
| LDDT_pli NaN on some ligands | OST graph-iso fails on connectivity mismatch | Bond-order assignment from SMILES (fix pending) |
| Catastrophic RMSD on multi-chain targets | Homodimer concat in input builder | Rewrote builder with per-chain dedup |
| Stale prep cache after builder change | `run_batch.sh` skip-if-exists optimization | Wipe prep cache after any builder change |

## Ground truth for CASP17

For CASP17 the architecture doesn't change — only the target loaders. The `inputs/casp17/` directory waits for the target release; once targets drop, `build_canonical.py` populates it from CASP's TSV + PDB bundle, adapters fan out to each predictor, and the batch runner fires.

Time from "CASP targets released" to "predictions on a GPU": minutes, not days. That's the entire point of the canonical inputs tree — it turns CASP17 participation into a drop-and-run operation.
