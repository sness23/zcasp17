# zcasp17 — an independent, fully-open pipeline for CASP17

**Benchmarking the open-weight AlphaFold-3 class of models on the hardest protein–ligand targets in biology — and submitting as an independent team at CASP17 in 2026.**

---

## What this is

`zcasp17` is the public home of a one-person research effort to:

1. **Benchmark every open-weight AF3-class structure predictor** (Chai-1, Protenix, Boltz) against the CASP15 and CASP16 ligand tracks using the exact scoring tools CASP uses — producing the first apples-to-apples, reproducible, cheap-to-rerun comparison.
2. **Submit to CASP17** (running in 2026) as a fully independent, non-institutional team using that pipeline.
3. **Ship the infrastructure** — canonical target definitions, pipeline adapters, scoring harnesses, and a ligand-pose viewer — so any researcher can reproduce, extend, or swap in a new predictor with one adapter script.

The work is already producing results that compete with the organized CASP groups. This repo is the consolidation layer: the plan, the receipts, the code, the CASP17 submission.

---

## Headline results (so far)

From runs against already-published CASP15 and CASP16 targets, scored with the same OpenStructure `compare-ligand-structures --lddt-pli` tool CASP uses:

| Benchmark | Target | Model | LDDT_pli | CASP rank | Beats groups |
|---|---|---|---|---|---|
| CASP15 | **T1124** (SAH / methyltransferase) | Chai-1 | **0.957** | **1 / 1570** poses | **26/26** |
| CASP15 | **H1135** (K⁺ / 12-chain assembly) | Chai-1 | **1.000** | **1 / 1676** | 12/19 |
| CASP15 | T1158v4 (ATP / ABC transporter) | Chai-1 | 0.960 | 38 / 746 | 18/24 |
| CASP16 | **L1001** (Chymase + inhibitor) | Chai-1 | **0.960** | **1 / 129** | **30/30** |
| CASP16 | L4017 (SARS-CoV-2 Mpro) | Chai-1 | 0.897 | 2 / 219 | 25/26 |

Of the 14 CASP16 Chymase targets that scored, **14 of 14** beat the CASP16 median-of-best. Four landed in the top-5 of 129 submissions.

Multi-model coverage (Protenix-v1 + Chai-1) is complete for CASP15 and in flight for CASP16. Near-native predictions (LDDT_pli ≥ 0.9) on **6/12 CASP15 targets** and **6/11 scored CASP16 targets** for Protenix-v1.

**Total compute cost across every benchmark run, including exploration and spot-tests: ~$2.90.**

Full per-target tables, receptor-by-receptor breakdowns, and failure analyses: [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

---

## The homodimer bug — a silent failure mode we found and fixed

Mid-run through the CASP16 pharma benchmark, we noticed that 14 of 24 SARS-CoV-2 Mpro targets were placing the ligand 30–50 Å away from the real binding pocket — catastrophic misses, not near-misses. But the model's `ranking_score` was still high (0.82+).

Investigation found a bug in the CASP16 input builder: homodimer chains from `protein_aligned.pdb` were being concatenated into a single "fake-monomer" sequence. The model was then dutifully folding a chimeric protein that doesn't exist in nature, with the ligand placed in a non-pocket location.

**The same bug exists unfixed in the canonical chai-lab CASP16 pipeline** and likely affected CASP16 submissions from other groups as well.

**Validation of the fix** (LDDT_pli / BiSyRMSD before → after):

| Target | Pre-fix | Post-fix | Δ |
|---|---|---|---|
| L2001 | 0.00 / 50 Å | **0.75 / 1.9 Å** | rescued |
| L4001 | 0.01 / 47 Å | **0.78 / 1.8 Å** | rescued |
| L4003 | 0.00 / 44 Å | 0.46 / 5.2 Å | partial rescue |

Full case study: [`docs/HOMODIMER-BUG.md`](docs/HOMODIMER-BUG.md).

This is the kind of finding that only comes out of running multiple independent models on the same inputs with an instrumented scoring pipeline — a core design choice of this project.

---

## Architecture

```
                     ┌────────────────────────────────────────┐
                     │      CANONICAL INPUTS (260 targets)    │
                     │  inputs/{casp15,casp16,casp17}/<T>/    │
                     │      target.json + receptor.pdb        │
                     │           + ligand_*.sdf               │
                     └───────────┬────────────────────────────┘
                                 │
                     ┌───────────┴────────────────┐
                     │        ADAPTERS            │
                     │   to_protenix_json.py      │
                     │   to_chai_fasta.py         │
                     │   to_boltz.py   (planned)  │
                     │   to_dynamicbind.py (stub) │
                     └───────────┬────────────────┘
                                 │
               ┌─────────────────┴─────────────────┐
               │                                   │
       ┌───────▼────────┐                 ┌────────▼────────┐
       │   Protenix     │                 │    Chai-1       │
       │   v1.0.0       │                 │   chai-lab      │
       └───────┬────────┘                 └────────┬────────┘
               │                                   │
               └─────────────────┬─────────────────┘
                                 │
                                 ▼
                     ┌────────────────────────────┐
                     │  UNIFIED SCORING HARNESS   │
                     │  OpenStructure 2.x         │
                     │  compare-ligand-structures │
                     │  --lddt-pli --rmsd         │
                     │  --substructure-match      │
                     └───────────┬────────────────┘
                                 │
                                 ▼
                     ┌────────────────────────────┐
                     │   CROSS-MODEL COMPARISON   │
                     │   + CASP RANKING JOIN      │
                     │   + casp-viewer (Svelte)   │
                     └────────────────────────────┘
```

Key design properties:

- **Resumable** — batch runners skip targets with completed outputs; pods can be killed and restarted.
- **Idempotent** — regenerating a FASTA or a prepped reference produces the same bytes.
- **Reproducible** — every target's FASTA/JSON, weights cache, and scoring output is committed to disk.
- **GPU-agnostic** — A40, L40S, H100, 3090, 4090 have all been used; the pipeline doesn't care.
- **Cheap** — a CASP-scale benchmark is $1–15, not $X,000.

Full architecture: [`docs/PIPELINE.md`](docs/PIPELINE.md).

---

## Canonical inputs tree — the infrastructure payoff

The main engineering contribution is not any single benchmark — it's the **canonical input tree** (`inputs/manifest.csv` with 260 targets), with one adapter per prediction tool.

```
inputs/
├── manifest.csv                 # master index — 260 targets
├── casp15/<T>/target.json       # 27 CASP15 targets
├── casp16/<T>/target.json       # 233 CASP16 targets
├── casp17/<T>/target.json       # (populates when CASP17 targets drop)
└── adapters/
    ├── to_protenix_json.py
    ├── to_chai_fasta.py
    └── to_dynamicbind.py
```

Adding a new prediction tool = writing one adapter script. Onboarding CASP17 = dropping files into `inputs/casp17/` and running the adapters.

This is what makes CASP17 participation feasible as a solo effort: when the targets drop, the time from "targets released" to "predictions on a GPU" is minutes, not days.

---

## Cost transparency

Every GPU-hour of this project is tracked. Cost is part of the scientific claim: that **a single independent researcher with $15 can reproduce CASP-scale benchmarks** that used to require institutional compute.

| Phase | GPU | Wall time | Cost |
|---|---|---|---|
| CASP15 Chai-1 full (11 targets) | A40 + H100 NVL | 95 min | ~$3.00 |
| CASP15 Protenix full (12 targets) | L40S 48GB | 108 min | $1.54 |
| CASP16 Chai-1 pharma PoC (44) | L40S 48GB | 2 hr | $1.80 |
| CASP16 Protenix partial (14) | L40S 48GB | 65 min | $0.93 |
| Post-fix validation spot-tests | 3090 + 4090 | 43 min | $0.43 |
| **Total to date** | | | **~$7.70** |

Projected for CASP17:

| Phase | Estimate |
|---|---|
| CASP17 live submission on all tractable targets | $10–$30 |
| Retrospective cross-model (Chai + Protenix + Boltz) | $30–$50 |

A full CASP-scale cross-model benchmark costs less than a tank of gas.

---

## Status

- [x] CASP15 Chai-1 benchmark complete (11 targets, scored)
- [x] CASP15 Protenix-v1 benchmark complete (12 targets, scored)
- [x] CASP16 Chai-1 pharma PoC complete (44 targets, 40 scored)
- [x] CASP16 Protenix-v1 PoC in progress (14 scored, 32-target post-fix rerun pending)
- [x] Homodimer-bug discovered, fixed, validated
- [x] Canonical inputs tree live (260 targets)
- [x] Protenix + Chai adapters
- [x] Scoring harness (OpenStructure LDDT_pli + BiSyRMSD)
- [x] CASP17 spider operational, target directory structure ready
- [ ] Boltz-2 adapter
- [ ] Boltz-2 retrospective on CASP15/16
- [ ] CASP17 live participation (targets release mid-2026)
- [ ] `casp-viewer` Svelte app — public beta
- [ ] Preprint (targeting bioRxiv, `q-bio.BM`)

---

## Why this matters

1. **The field needs independent baselines.** Every paper releasing a new structure predictor reports its own numbers on its own slice of data. A single external, reproducible, cross-model benchmark run with public scoring is missing — and it's the fastest way for an academic group to answer "which open model should we use?"

2. **CASP needs participants who can ship fast.** CASP17 will release targets on a tight schedule. A fully-automated pipeline that can go from target release to scored prediction in hours — without institutional MSA queues or committee approvals — is a genuine contribution to the competition even if it doesn't win.

3. **Open tooling compounds.** Each of: the canonical input tree, the OST scoring harness, the Protenix/Chai adapters, and the `casp-viewer` visualization app, is reusable outside CASP. They're the bricks other people's research will be built on.

4. **The cost story is important.** Lowering the marginal cost of a CASP-scale benchmark from "thousands of dollars + months of engineering" to "$15 + an afternoon" changes who can do this work. Graduate students. Independent researchers. Pharmaceutical teams evaluating open models. Anyone who wants to actually check a claim.

---

## Reproducing the benchmarks

The benchmarks that produced the headline numbers live in companion directories (being migrated into this repo):

- `casp15_ligands/` — Chai-1 CASP15
- `casp15_ligands_protenix/` — Protenix-v1 CASP15
- `casp16_ligands/` — Chai-1 CASP16 pharma
- `casp16_ligands_protenix/` — Protenix-v1 CASP16 pharma
- `inputs/` — canonical definitions + adapters
- `casp-viewer/` — Svelte web app for pose visualization

Each has its own `README.md` with full reproduction steps. A fresh CASP15 Chai-1 rerun: ~2 hr on an A40, ~45 min on an H100, total ~$3.

Generic recipe:

```bash
# 1. Clone canonical inputs
cd inputs/ && python3 build_canonical.py

# 2. Emit pipeline-specific inputs
python3 adapters/to_chai_fasta.py \
  inputs/casp17/*/target.json \
  --out-dir casp17_ligands/fastas/

# 3. Rent a GPU (RunPod / Vast.ai)
./run_batch.sh manifest.txt

# 4. Score
python3 score_lddt_pli.py --refs refs/ --results results/
```

See [`docs/PIPELINE.md`](docs/PIPELINE.md) for the full methodology.

---

## CASP17 plan

CASP17 targets release in mid-2026. Plan:

1. **T-minus one month**: finalize Boltz-2 adapter so all three open AF3-class models run through the same pipeline.
2. **Target release day**: spider pulls target list; canonical `target.json` written for each; all three predictor inputs emitted automatically.
3. **Folding**: parallel runs on rented GPUs. Manifest-partitioned by token count — A40 for ≤ 1024, L40S for 1536, H100 for 2048+.
4. **Submission**: top-5 poses per target per predictor; submissions deadline-gated by the CASP schedule.
5. **Retrospective**: cross-model comparison + CASP-official-score join + preprint within 30 days of results release.

Full plan + timeline: [`docs/CASP17-PLAN.md`](docs/CASP17-PLAN.md).

---

## Author

**Steven Ness** — independent researcher, PhD in computer science. Based in Canada.

Previous writing:
- [*A New Kind of Engineer*](https://medium.com/@sness23/a-new-kind-of-engineer-b0d533586a45) — on the emerging role of the "structural biology engineer"
- Medium: [@sness23](https://medium.com/@sness23)
- GitHub: [sness23](https://github.com/sness23)

This is a fully independent research effort. If you are:
- A **CASP17 observer, organizer, or competing team** — I'd love to coordinate on scoring or shared infrastructure.
- A **funder or program officer** (Emergent Ventures, Schmidt Sciences, Wellcome Leap, Astera) — the cross-model benchmark and CASP17 submission need GPU + living costs; cost breakdown above.
- An **academic lab** choosing between open AF3-class predictors — the cross-model comparison in `docs/BENCHMARKS.md` is directly useful. Happy to answer questions.
- A **model developer** (Chai Discovery, ByteDance Protenix, Boltz team) — I am running your tool against independent benchmarks and publishing. Please reach out.

Email: `sness@sness.net`.

---

## License

MIT — see [LICENSE](LICENSE). Predictor and scoring tool licenses differ; check upstream.

---

## Citation

If this work is useful to you, please cite:

```
@software{ness_zcasp17_2026,
  author  = {Ness, Steven},
  title   = {zcasp17: An independent, fully-open pipeline for CASP17},
  year    = {2026},
  url     = {https://github.com/sness23/zcasp17}
}
```

A preprint describing the CASP15/16 cross-model benchmarks + homodimer-bug case study is in preparation for bioRxiv / arXiv `q-bio.BM`.
