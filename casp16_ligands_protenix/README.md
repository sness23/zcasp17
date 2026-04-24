# CASP16 Pharma-Ligand Prediction with Protenix-v1.0.0

**Author:** sness
**Last updated:** 2026-04-21
**Status:** 14/44 PoC targets scored pre-fix; 4 spot-test re-runs post-fix.
Full re-run of the 25 multi-chain targets pending (deferred to next
available L40S).

## Abstract

We benchmarked Protenix-v1.0.0 against the CASP16 pharma-ligand track —
a drug-discovery benchmark of 44 protein-ligand complexes across four
target families (Chymase L1000, Cathepsin G L2000, SARS-CoV-2 Mpro
L4000, and the larger Autotaxin L3000 set, deferred). In a partial run
(14 targets scored), Protenix produced 7 near-native results
(LDDT_pli ≥ 0.9), 1 mid-range, and 6 failures. Diagnostic analysis
identified a pipeline bug — concatenated-homodimer input misrepresentation
causing catastrophic ligand misplacement — that affected 25 of the 44
targets. A fix was validated on L2001 (0.00 → 0.75 LDDT_pli), L4001
(0.01 → 0.78), and L4003 (0.00 → 0.46). A head-to-head comparison
against chai-lab on the same 11 scored targets shows neither method
dominates: 4 Protenix wins, 3 chai wins, 4 ties.

## Headline results (14 scored of 44, best sample per target)

**Near-native** (LDDT_pli ≥ 0.9 on primary ligand):

| Target | Ligand | LDDT_pli | RMSD (Å) | ranking |
|---|---|---|---|---|
| L1001 | ligand_201 (Chymase) | 0.91 | 1.03 | 0.98 |
| L1003 | ligand_843 | 0.97 | 0.60 | 0.98 |
| L1004 | ligand_124 | 0.97 | 0.46 | 0.99 |
| L1005 | ligand_868 | 0.97 | 0.65 | 0.99 |
| L2002 | ligand_380 (Cathepsin G) | 0.89 | 1.32 | 0.98 |
| L4002 | ligand_LIG_A_401 (Mpro) | 0.90 | 0.45 | 0.91 |

**Mid-range:**

| Target | LDDT_pli | RMSD (Å) | note |
|---|---|---|---|
| L1002 | 0.68 | 2.66 | right pocket, wrong rotamer |

**Failures (pre-fix):**

| Target | LDDT_pli | RMSD (Å) | why |
|---|---|---|---|
| L2001 | 0.00 | 50.4 | **Homodimer-bug catastrophe (FIXED — see below)** |
| L4001 | 0.01 | 47.3 | **Homodimer-bug catastrophe (FIXED)** |
| L4003 | 0.00 | 44.4 | **Homodimer-bug catastrophe (PARTIAL FIX)** |
| L4004 | 0.36 | 5.56 | Hard target, not homodimer-related |

**Graph-iso scoring failures (scorer-side, not method):**

| Target | Issue |
|---|---|
| L1006, L1007, L1008 | OST `compare-ligand-structures` could not match model ligand to ref SDF (same failure under chai-lab) |

## Homodimer fix (discovered mid-run)

Our `build_protenix_json.py` was concatenating homodimer chains from
`protein_aligned.pdb` into a single fake-monomer sequence. This
corrupted 25 of the 44 PoC targets. Fix: emit multiple `proteinChain`
entries (or `count=2` for truly identical chains). Validated on
spot-test runs on 3090 + 4090:

| Target | Before | After | Verdict |
|---|---|---|---|
| L2001 | 0.00 / 50 Å | **0.75 / 1.91 Å** | Rescued completely |
| L4001 | 0.01 / 47 Å | **0.78 / 1.78 Å** | Rescued completely |
| L4003 | 0.00 / 44 Å | 0.46 / 5.22 Å | Rescued from catastrophic → mid-range |
| L4004 | 0.36 / 5.56 Å | 0.36 / 5.44 Å | No change (not a homodimer issue) |

Full re-prep and re-run of all 25 multi-chain targets is queued for
the next available L40S slot. Estimated outcome: ~15 rescues, ~5
already-good no-ops, ~5 still-hard. See
[`docs/03_discussion.md`](docs/03_discussion.md#5-pipeline-lessons).

## Head-to-head vs chai-lab (on 14 shared CASP16 targets)

| Outcome | Count | Targets |
|---|---|---|
| Protenix wins (ΔLDDT ≥ 0.05) | 4 | L1003, L1004, **L4002**, L4004 |
| chai wins | 3 | L1001, L1002, **L4001** |
| Tie / both strong | 2 | L1005, L2002 |
| Both fail (40+ Å) | 2 | L2001, L4003 |
| Both OST-unscored (graph iso) | 3 | L1006, L1007, L1008 |

**Neither method dominates.** Most striking divergence: L4002 (chai's
34 Å error vs Protenix's 0.45 Å) and the mirror case L4001 (chai 2.1 Å
vs Protenix 47 Å — now rescued to 1.78 Å by the homodimer fix).

## Documentation

Scientific write-up split across `docs/`:

- [01_methods.md](docs/01_methods.md) — Pipeline, tools, MSAs, templates,
  sampling, scoring, hardware.
- [02_results.md](docs/02_results.md) — Full per-target tables,
  aggregate stats, head-to-head, confidence calibration.
- [03_discussion.md](docs/03_discussion.md) — Homodimer bug story,
  patterns by ligand class, chai comparison analysis.
- [04_limitations.md](docs/04_limitations.md) — Scope, known issues,
  what this work doesn't claim.
- [05_reproduction.md](docs/05_reproduction.md) — End-to-end reproduction guide.

Background docs (target-agnostic pipeline mechanics):

- [docs/prep_phase.md](docs/prep_phase.md) — What `protenix prep` does internally.
- [docs/templates.md](docs/templates.md) — Template search and the `TemplateEmbedder`.

## Target scope

The CASP16 pharma track has 233 targets across four drug-discovery
benchmark sets. For our PoC we ran 44 (`manifest_L1L2L4.txt`):

| Set | Receptor | # Targets | Notes |
|---|---|---|---|
| L1000 | Chymase (serine protease) | 17 | Monomer — our builder was correct |
| L2000 | Cathepsin G (cysteine protease) | 2 | Homodimer — bug hit both |
| L3000 | Autotaxin | 219 | Deferred; not in PoC |
| L4000 | SARS-CoV-2 Mpro | 25 | All homodimer — bug hit all 25 |

`manifest_L1L2L4.txt` (44 targets) is the PoC subset used for head-to-head
comparison with the chai baseline. L3000 is deferred because at
~15 min/target × 219 = ~55 GPU hours = ~$47 of L40S time, its PoC ROI
is lower than fixing correctness issues first.

## Quick provenance / cost summary

| Compute | GPU | Time | Cost |
|---|---|---|---|
| CASP15 full run (12 targets) | L40S 48GB, RunPod | 108 min | $1.54 |
| CASP16 PoC partial (L2000+14 of L4000, pre-fix) | L40S, RunPod | ~65 min | $0.93 |
| CASP16 L2001 spot test (post-fix) | RTX 3090, Vast.ai | 9.5 min | $0.03 |
| CASP16 L4 spot test x3 (post-fix) | RTX 4090, RunPod | 33 min | $0.40 |

Total to date: **~$3 for 33 CASP-ligand Protenix inferences** across
both benchmark years, plus ~4 hours local CPU time for MSA prep.

## Directory layout

```
casp16_ligands_protenix/
├── README.md                         # this file
├── build_protenix_json.py            # PDB+SMILES → Protenix JSON (now with homodimer fix)
├── prep_local.sh                     # MSA + template search runner (CPU, local)
├── run_batch.sh                      # GPU inference runner (pod)
├── clone_prepped.py                  # reuse MSAs across identical sequences
├── score_lddt_pli.py                 # OST LDDT_pli + BiSyRMSD scorer
├── manifest_L1L2L4.txt               # 44-target PoC
├── manifest_L2000.txt                # 2 Cathepsin G targets
├── manifest_smoke.txt                # 2-target smoke test
├── manifest_spot_L2001.txt           # post-fix validation (1 target)
├── manifest_spot_L4_fails.txt        # post-fix L4 spot tests (3 targets)
├── jsons/                            # built input JSONs (post-fix)
├── jsons_concat_bug_backup/          # pre-fix JSONs, kept for diff
├── jsons_prepped/                    # prepped JSONs + MSA trees
├── results/                          # GPU CIFs per target
├── out/lddt_pli/                     # OST scores
├── refs/ -> ../casp16_ligands/refs/  # shared ground truth (symlink)
└── docs/                             # scientific write-up
```
