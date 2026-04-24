# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This directory runs the **chai-lab** protein-ligand folding pipeline against the CASP15 ligand targets. The CASP release bundle `targets_ligand/casp15.targets.ligands.ALL_09.18.2025.tar.gz` unpacks to `targets_ligand/Targets_ligand/<TARGET>_lig.pdb` — these are the ground-truth structures used both as *input specifications* (to derive chain sequences + ligand CCD codes) and later as *references* to score predictions against.

The pipeline is designed to run on RunPod GPU pods: a local machine builds inputs + MSAs, the pod runs chai-lab fold, then results come back here for scoring.

## Pipeline Stages

```
  targets_ligand/*_lig.pdb
          │
          ▼  build_chai_fasta.py   (gemmi + RCSB CCD → SMILES)
      fastas/<TARGET>.fasta          ← protein chains + ligand SMILES entries
          │
          ▼  prefetch_msas.py       (ColabFold MMseqs2 via chai-lab wrapper)
      msas/<TARGET>/*.aligned.pqt    ← optional; falls back to --use-msa-server
          │
          ▼  run_batch.sh           (on GPU pod, one manifest)
      results_A40/<TARGET>/{pred.model_idx_0-4.cif, scores.model_idx_0-4.npz, msa_depth.pdf}
          │
          ├─▶ summarize_results.py  — ranked aggregate / pTM / ipTM table
          ├─▶ sanity_check.py       — file presence, seq integrity, pLDDT, inter-model RMSD, clashes
          └─▶ compare_to_pdb.py     — Kabsch CA RMSD + ligand centroid Δ + 4 Å footprint overlap vs crystal
```

## Common Commands

```bash
# Build all FASTAs from *_lig.pdb (ligand SMILES cached in .smiles_cache.json)
python3 build_chai_fasta.py --all
python3 build_chai_fasta.py T1124                  # one target → stdout

# Pre-fetch ColabFold MSAs for a manifest (CPU; no chai-lab weights needed)
python3 prefetch_msas.py manifest_A40.txt

# Cache all chai-lab model weights locally so RunPod pods can rsync them in
./fetch_chai_downloads.sh                          # → ~/data/chai_downloads (~10 GB)

# On the pod: fold everything in a manifest. Resumable (skips dirs with 5 CIFs)
./run_batch.sh manifest_A40.txt
FASTA_DIR=/workspace/fastas OUT_ROOT=/workspace/results ./run_batch.sh manifest_A100.txt

# Analyze results
python3 summarize_results.py results_A40 --md
python3 sanity_check.py     results_A40 fastas
python3 compare_to_pdb.py   results_A40 -v
```

## Manifests: GPU Sizing

Targets are partitioned by expected chai-lab token count, since chai-lab memory scales super-linearly with tokens:

- `manifest_A40.txt` — ≤1024 tokens, fits A40 48 GB (7 targets: T1124/T1127/T1127v2/T1146/T1152/T1187/T1188)
- `manifest_A100.txt` — 1536–2048 tokens, needs A100 80 GB (T1158v1-v4, H1135)

RNA-only targets (`R*`) are excluded from both.

## Key Gotchas

- **chai-lab requires an empty output dir.** `run_batch.sh` does `rm -rf "$out"` before each fold (unless the skip-if-complete check fires). Don't hand-edit files inside a target results dir expecting them to survive a rerun.
- **MSA strategy is auto-detected.** If `msas/<TARGET>/*.aligned.pqt` exists, `run_batch.sh` passes `--msa-directory`; otherwise falls back to `--use-msa-server` (live ColabFold API).
- **SMILES come from RCSB CCD**, cached in `.smiles_cache.json`. `build_chai_fasta.py` skips waters and (defensively) `TYR` when it appears as HETATM. Unresolved CCD codes raise and skip the target with a warning.
- **Chain matching in `compare_to_pdb.py` is by longest exact sequence prefix** (≥20 residues) — crystal renumbering or mutations will degrade the match. Ligand comparison is centroid-distance + 4 Å footprint overlap, *not* atom-by-atom RMSD (no symmetry-aware substructure matching).
- **pLDDT is read from the CIF B-factor column** in `sanity_check.py` (chai-lab writes it there).

## Data Sources

- `targets_ligand/Targets_ligand/*_lig.pdb` — CASP15 ground-truth structures (protein + ligand)
- `pdb/*.pdb` — companion PDB entries for targets that map to deposited structures
- `ligand_targets_mapping.csv`, `ligand_targets_with_dois.csv/.json` — target → PDB code → protein/ligand metadata
- `binding_sites.md`, `binding_descriptions.md`, `pdb_ligand_analysis.md` — curated notes on binding sites per target
- `A40_summary.md` — checked-in snapshot of the latest A40 batch summary

## Dependencies

- `gemmi` (structure I/O), `numpy` — used by almost every script
- `chai-lab` (only for `prefetch_msas.py` and pod-side folding; CPU-only install is enough for MSA prefetch)
- `requests`, `beautifulsoup4`, `markdownify` — inherited from parent `casp/` spider (not used by the folding pipeline itself)
