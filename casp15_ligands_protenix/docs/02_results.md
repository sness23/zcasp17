# 02 — Results

All numbers below are from the "best sample per target" selection — the
sample (across 5 seeds × 5 samples each) with the highest
`ranking_score` from Protenix's confidence head.

## 1. Per-target table

Full scoring output is in `out/lddt_pli/summary.csv`. Key columns
reproduced here, grouped by ligand class within each target.

### Near-native (LDDT_pli ≥ 0.9 on the primary ligand)

| Target | Ligand (n) | LDDT_pli | RMSD (Å) | lddt_lp | ranking |
|---|---|---|---|---|---|
| **T1124** | SAH (2) | 0.97, 0.93 | 0.14, 0.33 | 0.92, 0.98 | 0.97 |
| **T1127** | EPE (2) | 0.96, 0.94 | 0.43, 0.58 | 0.95, 0.97 | 0.97 |
| **T1127v2** | EPE (2) | 0.96, 0.94 | 0.42, 0.58 | 0.97, 0.95 | 0.97 |
| **T1158v4** | ATP (2) | 0.95, 0.96 | 0.37, 0.31 | 0.92, 0.91 | 0.81 |
| **T1158v4** | Mg²⁺ (2) | 1.00, 0.98 | 0.27, 0.15 | 1.00, 0.96 | 0.81 |
| **H1135** | K⁺ (9 of 9) | 1.00 each | 0.08–0.22 | 0.97–0.99 | 0.60 |

### Mid-range (LDDT_pli 0.5–0.9)

| Target | Ligand | LDDT_pli | RMSD (Å) | ranking |
|---|---|---|---|---|
| T1146 | NAG | 0.84 | 0.40 | 0.98 |
| T1152 | NAG | 0.52 | 1.67 | 0.62 |
| T1127 / T1127v2 | COA (×2 each) | ~0.57, ~0.52 | ~4.70, ~5.45 | (same target as above) |
| T1158v1 | XPG | 0.65 | 2.09 | 0.90 |
| T1158v2 | P2E | 0.66 | 2.25 | 0.90 |
| T1188 | DW0 (×2) | 0.62, 0.50 | 2.46, 4.96 | 0.94 |
| T1188 | Cd²⁺ (1 of 2) | 0.27 | 5.89 | (same) |

### Failures (LDDT_pli ≤ 0.2 or RMSD > 10 Å)

| Target | Ligand | LDDT_pli | RMSD (Å) | ranking | Why |
|---|---|---|---|---|---|
| T1127 / T1127v2 | MPD | 0.00 | ~21.7 | (high, 0.97) | Small solvent additive misplaced |
| T1158v3 | XH0 | — | — | 0.90 | OST graph-iso match failed |
| T1187 | NAG (×2) | 0.15, 0.14 | 5.9, 5.8 | **0.40** | Genuinely hard; ranking_score agrees |
| T1188 | Cd²⁺ (1 of 2) | 0.00 | 37.2 | 0.94 | Wrong binding pocket |
| T1188 | Co²⁺ | 0.00 | 31.2 | 0.94 | Wrong binding pocket |
| H1135 | Cl⁻ (×3) | 0.08, 0.01, 0.00 | 10.5, 7.3, 10.6 | 0.60 | Halide not captured despite K⁺ success |

## 2. Aggregate statistics

**Scoring coverage**

- Targets with at least one scored ligand: **11/12** (92%).
- Target lost to graph-iso failure: 1 (T1158v3 / XH0).
- Total model-reference ligand pairs compared: 39.

**Primary-ligand accuracy** (the "main" ligand per target, picking the
CCD code with the highest-count HETATM if multiple):

| Bucket | Targets | Share |
|---|---|---|
| Near-native (LDDT_pli ≥ 0.9) | 6 | 50% |
| Mid-range (0.5 – 0.9) | 3 | 25% |
| Poor (< 0.5) or failed | 3 | 25% |

**RMSD distribution (primary ligand):**

| Bucket | Count |
|---|---|
| ≤ 1 Å (atomically correct) | 5 |
| 1 – 2 Å (near-native) | 1 |
| 2 – 5 Å (right pocket, wrong rotamer) | 3 |
| > 5 Å (wrong site) | 2 |

## 3. Results by ligand class

Same underlying numbers, pivoted to see where Protenix is strong vs
weak:

| Class | Representative CCDs | Targets | Typical LDDT_pli | Observation |
|---|---|---|---|---|
| **Nucleotide cofactors** | ATP, SAH | T1124, T1158v4 | 0.93–0.97 | Consistently excellent; high sequence signal + large contact surface |
| **Cations (large)** | K⁺, Mg²⁺ | H1135, T1158v4 | 0.95–1.00 | Perfect when adjacent to well-resolved coordination site; score reflects very few contact atoms (lddt_pli_n_contacts ≈ 30–60) |
| **Zwitterionic buffers** | EPE (HEPES) | T1127, T1127v2 | 0.94–0.96 | Strong — buffer sites are often ordered in CASP targets |
| **Glycans (single)** | NAG | T1146, T1152, T1187 | 0.15–0.84 | Bimodal — T1146 near-native, T1187 lost. Glycosylation site detection is target-specific. |
| **Coenzymes (large)** | COA | T1127, T1127v2 | 0.52–0.57 | Right pocket, wrong fold — CoA's long tail is placed variably |
| **Drug-like novel** | XPG, P2E, DW0, XH0 | T1158v1-3, T1188 | 0.15–0.66 (or graph-iso failure) | Middling; matches CASP expectation that these are the hard targets |
| **Halide ions** | Cl⁻ | H1135 | 0.00–0.08 | Consistently misplaced |
| **Solvent/crystallization additives** | MPD | T1127, T1127v2 | 0.00 | Expected failure — small, weakly bound, often crystallographic artifact |
| **Transition metal ions** | Cd²⁺, Co²⁺ | T1188 | 0.00–0.27 | Misplaced (wrong binding site for 2 of 3 ions) |

## 4. Confidence calibration

Protenix's `ranking_score` is a well-calibrated indicator of
pose quality at the extremes, weaker in the middle:

| ranking_score | n targets | typical LDDT_pli | takeaway |
|---|---|---|---|
| ≥ 0.95 | 5 | ≥ 0.84 (all cases) | high-confidence predictions are reliable |
| 0.80 – 0.95 | 4 | 0.00 – 0.97 (noisy) | **not reliable** in this band |
| ≤ 0.60 | 3 | 0.15 – 1.00 (!) | paradoxical: T1187 low-conf & poor ✓; H1135 low-conf but K⁺ perfect — confidence here is diluted by the misplaced Cl⁻ ions, not a per-ligand score |

**Key observation:** `ranking_score` is a **target-level** confidence,
not per-ligand. For targets with heterogeneous ligand outcomes (H1135:
perfect K⁺ + wrong Cl⁻; T1188: placed DW0 + misplaced CO), a single
low score doesn't tell you which ligands are good. Operationally, use
per-chain pLDDT from `chain_plddt` in the summary-confidence JSON if
you need per-ligand confidence signals.

## 5. Timing and cost

Run on RunPod L40S (48 GB, $0.86/hr) with the paths baked into the
prepped JSONs pointing at `/workspace/casp15_ligands_protenix/`.

### A40 manifest (7 targets, ≤1024 tok)

| Target | Wall (s) |
|---|---|
| T1152 | 467 (first — includes ~100 s checkpoint download + LayerNorm JIT compile) |
| T1146 | 146 |
| T1187 | 168 |
| T1127 | 222 |
| T1127v2 | 271 |
| T1188 | 271 |
| T1124 | 337 |
| **Total** | **1882 s (31 min)** |

### A100 manifest (5 targets, 1536-2048 tok)

| Target | Wall (s) |
|---|---|
| T1158v1 | 756 |
| T1158v2 | 760 |
| T1158v3 | 785 |
| T1158v4 | 742 |
| H1135 | 1608 (biggest target — 2048 tok) |
| **Total** | **4651 s (77 min)** |

### Totals

| Metric | Value |
|---|---|
| Wall time (all 12) | ~108 min |
| GPU cost (L40S) | **~$1.54** |
| Local prep cost | $0 (runs on user workstation; MSA queue time dominates) |
| Total structures generated | 300 (12 × 25 samples) |
| Cost per structure | ~$0.005 |

## 6. What to look at next

- Re-score with `--substructure-match` alternate modes or preprocess
  the XH0 ref SDF to recover T1158v3.
- Run more seeds on T1187 and T1188 to probe whether the catastrophic
  failures are seed-specific or architectural.
- Run the excluded 13 CASP15 targets (T1170, T1181, T1186, H1114,
  H1171v1-v2, H1172v1-v4, R\*) to extend coverage.
- Compare head-to-head against chai baseline on CASP15 targets —
  scoring CSV exists in `../casp15_ligands/out/...` if the chai pipeline
  was run there; not yet directly aligned with this work.

See [03_discussion.md](03_discussion.md) for interpretation.
