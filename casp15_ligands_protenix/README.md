# CASP15 Ligand-Binding Prediction with Protenix-v1.0.0

**Author:** sness
**Last updated:** 2026-04-21
**Status:** 12/12 targets folded, 11/12 scored by OST LDDT_pli.

## Abstract

We reproduced the CASP15 ligand-binding benchmark (12 protein-ligand
complex targets) using the open-source Protenix-v1.0.0 structure
predictor as a drop-in replacement for chai-lab. Protenix produced
near-native poses (LDDT_pli ≥ 0.9, ligand BiSyRMSD < 0.7 Å) on 5 of 12
targets, and competitive poses (LDDT_pli ≥ 0.5) on 3 more. Hard
targets included carbohydrate-binding sites (T1187, partial on T1152)
and a novel drug-like cofactor (T1158v3) where OpenStructure's
ligand-graph matching failed to align the model ligand to the
crystallographic reference. Total compute cost: **~$1.54 on a single
RunPod L40S, ~77 minutes wall time for all 12 targets** (including
25 samples per target).

## Headline results

**Near-native** (top sample LDDT_pli ≥ 0.9 on the primary ligand):

| Target | Ligand | LDDT_pli | RMSD (Å) |
|---|---|---|---|
| T1124 | SAH (×2) | 0.97 / 0.93 | 0.14 / 0.33 |
| T1127 / T1127v2 | EPE (×2) | 0.96 / 0.94 | 0.43 / 0.58 |
| T1158v4 | ATP (×2), Mg²⁺ (×2) | 0.95 / 1.00 | 0.31 / 0.15 |
| H1135 | K⁺ (×9) | 1.00 | 0.09–0.22 |

**Partial success** (LDDT_pli 0.5–0.9 or some ligands hit, others missed):

| Target | Note |
|---|---|
| T1146 | NAG glycan — 0.84 LDDT, 0.40 Å |
| T1152 | NAG — 0.52 LDDT, 1.67 Å |
| T1158v1 | XPG novel ligand — 0.65 LDDT, 2.09 Å |
| T1158v2 | P2E novel ligand — 0.66 LDDT, 2.25 Å |
| T1188 | DW0 placed (0.62); CO²⁺ and one Cd²⁺ misplaced (~30 Å) |
| H1135 (Cl⁻) | 9 of 12 ions perfect (K⁺); 3 Cl⁻ misplaced ~7–10 Å |

**Failures:**

| Target | Issue |
|---|---|
| T1187 | NAG poorly placed — 0.15 LDDT, 5.9 Å. Model's own `ranking_score` of 0.40 correctly flags this as low-confidence. |
| T1158v3 | XH0 ligand — OST could not match via graph isomorphism (see [discussion](docs/03_discussion.md#ost-graph-isomorphism-failures)). |

## Documentation

Scientific write-up split across `docs/`:

- [01_methods.md](docs/01_methods.md) — Pipeline, tools, MSAs, templates, sampling protocol, scoring.
- [02_results.md](docs/02_results.md) — Per-target results, aggregate statistics, confidence calibration.
- [03_discussion.md](docs/03_discussion.md) — Patterns by ligand class, failure-mode analysis, comparison to CASP16 pharma set.
- [04_limitations.md](docs/04_limitations.md) — Scope limits, known biases, what this work doesn't claim.
- [05_reproduction.md](docs/05_reproduction.md) — Step-by-step to reproduce from scratch.

Background docs (carried over from the CASP16 Protenix pipeline):

- [docs/prep_phase.md](docs/prep_phase.md) — What `protenix prep` actually does (MSA server protocol, template search, output paths).
- [docs/templates.md](docs/templates.md) — How Protenix's `TemplateEmbedder` consumes structural templates.

## Target scope

The 25+ CASP15 ligand targets reduced to **12** by restricting to chai's
validated manifests (A40: 7 protein-ligand monomer/dimer targets at ≤1024
tokens; A100: 5 larger or oligomeric targets at 1536–2048 tokens).
RNA-only (R\*) targets and a handful of protein targets not in either
chai manifest (T1170, T1181, T1186, H1114, H1171v1–v2, H1172v1–v4) are
excluded for symmetry with the baseline.

| Manifest | Targets | Token range | GPU class |
|---|---|---|---|
| `manifest_A40.txt` | T1124, T1127, T1127v2, T1146, T1152, T1187, T1188 | ≤1024 | A40 / L40S / 4090 (24 GB works) |
| `manifest_A100.txt` | H1135, T1158v1–v4 | 1536–2048 | A100 / H100 (L40S 48 GB also works) |
| `manifest_all.txt` | all 12, sorted | mixed | any 48 GB class |

## Directory layout

```
casp15_ligands_protenix/
├── README.md                   # this file
├── build_protenix_json.py      # CASP15 PDB → Protenix JSON
├── prep_local.sh               # MSA + template search runner (CPU, local)
├── run_batch.sh                # GPU inference runner (pod)
├── clone_prepped.py            # reuse MSAs across targets w/ identical sequences
├── score_lddt_pli.py           # OST LDDT_pli + BiSyRMSD scorer
├── manifest_A40.txt            # 7 small/medium targets
├── manifest_A100.txt           # 5 large/oligomer targets
├── manifest_all.txt            # all 12
├── jsons/                      # built input JSONs (12)
├── jsons_prepped/              # prepped JSONs + MSA trees (12)
├── results/                    # GPU-generated CIFs (12 × 25 = 300)
├── out/lddt_pli/               # OST scoring output + summary.csv
├── refs -> ../casp15_ligands/refs  # shared ground-truth (symlink)
└── docs/                       # detailed write-up
```

## Data & code provenance

- **Model:** Protenix-v1.0.0 (`protenix_base_default_v1.0.0`, 368M params,
  AlphaFold3-equivalent data cutoff 2021-09-30). The ByteDance v2
  checkpoint is referenced in the Protenix repo but not publicly released
  (see [#295](https://github.com/bytedance/Protenix/issues/295)).
- **Scoring:** OpenStructure `compare-ligand-structures` via conda env `ost`,
  using `--lddt-pli --rmsd --substructure-match`.
- **Ground truth:** CASP15 bundle `Targets_ligand/<T>_lig.pdb`, split into
  `receptor.pdb` + `lig_*.sdf` by `../casp15_ligands/prep_references.py`.
- **Baseline:** chai-lab 0.6.x run on the same targets, same refs — see
  `../casp15_ligands/` for that pipeline (not directly compared in this
  document; the CASP16 comparison in `../casp16_ligands_protenix/` is
  more directly head-to-head).
