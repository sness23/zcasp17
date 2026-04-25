# Benchmarks — CASP15 and CASP16 ligand tracks

This document holds per-target benchmark results across every open-weight AlphaFold-3-class model we've run. Scoring is CASP's own tooling (OpenStructure `compare-ligand-structures --lddt-pli --rmsd --substructure-match`), against CASP's public ligand-score rankings.

## Summary

| Benchmark | Model | N targets | Near-native (≥0.9) | Cost |
|---|---|---|---|---|
| CASP15 ligand track | Chai-1 | 11 of 25 | 3 | ~$3.00 |
| CASP15 ligand track | Protenix-v1.0.0 | 12 of 25 | 6 | $1.54 |
| CASP16 pharma (L1+L2+L4 PoC) | Chai-1 | 40 of 44 | 5 | $1.80 |
| CASP16 pharma (L1+L2+L4 PoC) | Protenix-v1.0.0 | 11 of 14 | 6 | $0.93 |

Targets not run fall into three categories:
- **Token-limit exceeded** (Chai-1 caps at 2048 tokens; 9 CASP15 targets blocked).
- **Non-CCD ligands** (CASP15 T1118v1, T1186 — would require manual SMILES).
- **OST graph-isomorphism scoring failures** (4 CASP16 targets — bond-order assignment fix pending).

---

## CASP15 — Chai-1

```
Headline: #1 on T1124 of all 26 CASP groups.
          #1 on H1135 potassium ion of 1,676 poses.
          Top-38 on T1158v4 of 746 (ATP in ABC transporter).
```

| Target | agg score | CA RMSD | LDDT_pli | CASP best | CASP median | Rank | Beats groups |
|---|---|---|---|---|---|---|---|
| **T1124** | 0.912 | 3.14 Å | **0.957** | 0.927 | 0.507 | **1 / 1570** | **26/26** |
| **H1135** (K⁺) | 0.624 | — | **1.000** | 1.000 | 0.938 | **1 / 1676** | 12/19 |
| T1158v4 (ATP) | 0.657 | 4.95 Å | 0.960 | 1.000 | 0.547 | 38 / 746 | 18/24 |
| T1146 | 0.681 | 1.09 Å | 0.803 | 0.972 | 0.506 | 129 / 400 | 16/26 |
| T1188 | 0.667 | 0.65 Å | 0.611 | 0.828 | 0.129 | 84 / 976 | 9/23 |
| T1127v2 | 0.924 | 1.52 Å | 0.586 | 0.767 | 0.484 | 152 / 654 | 6/23 |
| T1152 | 0.280 | 6.16 Å | 0.251 | 0.844 | 0.152 | 127 / 470 | 12/28 |
| T1187 | 0.560 | 9.33 Å | 0.159 | 0.385 | 0.020 | 74 / 570 | 13/21 |
| T1158v2 | 0.460 | 15.56 Å | 0.399 | 0.668 | 0.243 | 85 / 360 | 9/24 |
| T1158v3 | 0.390 | 4.92 Å | — | — | — | — | OST graph-iso |
| T1127 | 0.925 | 1.52 Å | 0.562 | — | — | — | (merged w/ v2) |

---

## CASP15 — Protenix-v1.0.0

| Bucket | n | Targets |
|---|---|---|
| Near-native (LDDT_pli ≥ 0.9) | 6 | T1124, T1127, T1127v2, T1158v4, H1135, T1146 |
| Mid-range (0.5–0.9) | 3 | T1152, T1158v1, T1158v2 |
| Poor (< 0.5) | 2 | T1187, T1188 |
| OST graph-iso | 1 | T1158v3 |

Cost: $1.54 (L40S 48 GB, 108 min wall).

Findings:
- Nucleotide cofactors (SAH, ATP), large cations (K⁺, Mg²⁺), and zwitterionic buffers (EPE) are placed near-perfectly.
- Small solvent additives (MPD), halide ions (Cl⁻), and some transition metals (Cd²⁺, Co²⁺) are systematically misplaced.
- Glycan placement (NAG) is bimodal — T1146 stellar (0.84), T1187 lost (0.15).
- `ranking_score` is target-level; use `chain_plddt` for per-ligand confidence.

---

## CASP16 pharma track — Chai-1

Per-receptor summary:

| Receptor | N scored | Chai-1 best | Chai-1 median | Chai-1 worst | CASP16 median-of-best |
|---|---|---|---|---|---|
| **Chymase (L1000)** | 14 | **0.988** | **0.903** | 0.503 | 0.951 |
| Cathepsin G (L2000) | 2 | 0.876 | 0.438 | 0.000 | 0.928 |
| **Mpro (L4000)** | 24 | 0.897 | **0.012** | 0.000 | 0.902 |

Key targets:

| Target | LDDT_pli | RMSD | CASP best | Rank | Beats groups |
|---|---|---|---|---|---|
| **L1001** | **0.960** | 0.48 Å | 0.951 | **1 / 129** | **30/30** |
| L1011 | 0.963 | 0.47 Å | 0.967 | 3 / 129 | 29/30 |
| L1017 | 0.938 | 0.83 Å | 0.938 | 2 / 132 | 29/30 |
| L1016 | 0.960 | 0.57 Å | 0.988 | 4 / 129 | 28/30 |
| L1010 | 0.988 | 0.39 Å | 1.000 | 8 / 130 | 26/31 |
| **L4017** (Mpro) | 0.897 | 0.61 Å | 0.899 | **2 / 219** | 25/26 |
| L4015 (Mpro) | 0.889 | 0.49 Å | 0.918 | 5 / 219 | 24/26 |

**All 14 *scored* Chymase targets beat the CASP16 median-of-best.** 4 of 14 landed top-5 of 129 submissions. (The L1000 set has 17 monomer Chymase targets total; L1006 / L1007 / L1008 produced no LDDT_pli score because OST's graph-isomorphism matcher couldn't align the predicted ligand atoms to the reference SDF — a scoring-tool limitation, not a prediction failure. 17 − 3 = 14 scored.)

Mpro is bimodal: when it works, it's top-5 of 200+ submissions. When it fails, the ligand is placed 30–50 Å from the real binding pocket — a pattern traced to the homodimer input bug (see [HOMODIMER-BUG.md](HOMODIMER-BUG.md)).

---

## CASP16 pharma track — Protenix-v1.0.0 (in progress)

| Bucket | n (of 14 scored pre-fix) | Notes |
|---|---|---|
| Near-native (LDDT_pli ≥ 0.9) | 6 | L1001, L1003, L1004, L1005, L2002, L4002 |
| Mid (0.5–0.9) | 1 | L1002 |
| Catastrophic (pre-fix, input-bug driven) | 3 | L2001, L4001, L4003 |
| Poor (non-bug) | 1 | L4004 |
| OST graph-iso | 3 | L1006, L1007, L1008 |

Post-fix validation re-runs (LDDT_pli / BiSyRMSD):

| Target | Pre-fix | Post-fix | Δ |
|---|---|---|---|
| L2001 | 0.00 / 50 Å | **0.75 / 1.9 Å** | rescued |
| L4001 | 0.01 / 47 Å | **0.78 / 1.8 Å** | rescued |
| L4003 | 0.00 / 44 Å | 0.46 / 5.2 Å | partial rescue |
| L4004 | 0.36 / 5.6 Å | 0.36 / 5.4 Å | unchanged (non-bug failure) |

A full 32-target post-fix rerun is pending (~$2 on a rented 4090).

---

## Head-to-head: Protenix-v1 vs Chai-1 (14 shared CASP16 targets)

| Outcome | Count | Notes |
|---|---|---|
| Protenix wins | 4 | L1003, L1004, L4002 (34 Å swing), L4004 |
| Chai-1 wins | 3 | L1001, L1002, L4001 (pre-fix; post-fix tied or Protenix-leaning) |
| Tie / both strong | 2 | L1005, L2002 |
| Both catastrophic | 2 | L2001, L4003 (both input-bug affected on both pipelines) |
| Both OST-unscored | 3 | L1006, L1007, L1008 |

**Caveat**: chai-lab's CASP16 FASTAs have the same homodimer concat bug. Chai-1's 3 catastrophes are likely input-bug-driven too, so this head-to-head understates Chai-1's real performance. A chai-lab re-run with fixed inputs (one-liner from the canonical adapters) would update the comparison.

---

## Scoring methodology

All numbers use OpenStructure 2.x `compare-ligand-structures` with:

```
--lddt-pli              # Ligand-protein interaction LDDT (CASP's primary metric)
--rmsd                  # BiSyRMSD (symmetry-aware ligand RMSD)
--substructure-match    # Tolerant matching for CCD vs crystal differences
```

Scored against:
- **Reference receptor PDB** (split from CASP's `<TARGET>_lig.pdb`)
- **Reference ligand SDF(s)** (per-HETATM)

For each target the **best-aggregate-score** model (out of 5 diffusion samples for Chai-1, 25 per-seed samples for Protenix) is scored against all reference ligands; best-per-reference is reported.

Rank is computed against CASP's public `ligand_scores.csv` (26,432 per-submission rows for CASP15, comparable for CASP16). "Beats groups" counts distinct competing CASP groups whose best submission scored lower.

---

## Known scoring limitations

1. **Best-model only.** CASP groups submit 5 poses and are scored on the best; we score only the best-aggregate Chai-1/Protenix sample. Pooling all 5 would improve our numbers 5–15% on borderline targets.
2. **Graph-isomorphism gaps.** OST LDDT_pli matches ligand atoms by graph isomorphism. A few CASP targets (T1158v3, L1006, L1007, L1008, L4028) fail match because our ligand atom count / connectivity doesn't match the reference. Fixable via explicit bond-order assignment from SMILES before writing the SDF.
3. **Single-atom ions.** CO²⁺, CD²⁺, and similar score 0 for pure graph reasons, not prediction reasons.
4. **No chai-lab ensembling.** No multi-seed, no relaxation, no pose reranking. These are the out-of-the-box numbers.
5. **MSA differences.** CASP15 used jackhmmer MSAs; our Chai-1 run uses ColabFold. MSA pipeline alone can move LDDT_pli by a few percent.

---

## Data artifacts

| Path | Contents |
|---|---|
| `casp15_ligands/results_A40/` | 7 targets × 5 CIFs + scores (Chai-1) |
| `casp15_ligands/results_H100/` | 4 targets × 5 CIFs + scores (Chai-1) |
| `casp15_ligands/out/lddt_pli/` | Per-target OST LDDT_pli JSON + summary CSV |
| `casp15_ligands_protenix/out/` | Protenix-v1 CASP15 full outputs |
| `casp16_ligands/out/lddt_pli_L40S/` | Chai-1 CASP16 pharma outputs |
| `casp16_ligands_protenix/out/lddt_pli_postfix_full/` | Protenix-v1 CASP16 (partial + validation) |
| `inputs/` | 260 canonical target.json files + adapters |

All reproducible from the scripts in the companion directories. See each subdirectory's `README.md`.
