# CASP17 participation plan

## The goal

Submit to CASP17 (2026) as an independent, fully-open team using the zcasp17 pipeline. Produce the first cross-model (Chai-1 + Protenix + Boltz-2) CASP submission. Publish a post-competition retrospective comparing per-model performance against the official CASP17 rankings, within 30 days of results release.

## Why this is doable as a solo effort

The infrastructure cost of CASP submission collapses to near-zero when:

- Input construction is automated (canonical tree + adapters).
- MSA generation is outsourced (ColabFold / Protenix MSA server).
- GPU provisioning is on-demand (RunPod / Vast.ai).
- Scoring is reproducible (OST).

All four are already built and tested against CASP15 and CASP16. CASP17 is a drop-in target set, not new infrastructure work.

## Timeline

| Date | Milestone |
|---|---|
| 2026-04 | CASP17 spider operational; registration submitted |
| 2026-05 | Boltz-2 adapter finished; third model live in pipeline |
| 2026-06 | Cross-model retrospective on CASP15/16 published (preprint) |
| ~2026-06 to ~2026-09 | CASP17 target release window (historically ~4 months) |
| T+24 hr per batch | Submit predictions, per predictor, per tractable target |
| 2026-12 | CASP17 results release; full retrospective published |

## Per-target submission plan

For each CASP17 target that drops:

1. **Spider ingests** — target metadata, receptor sequence, ligand SMILES (or CCD). Canonical `target.json` written to `inputs/casp17/<T>/`.
2. **Adapters fan out** — `to_protenix_json.py`, `to_chai_fasta.py`, `to_boltz.py` emit per-predictor inputs.
3. **Manifest partition** — targets bucketed by token count (≤1024 → A40; 1536 → L40S; 2048 → H100).
4. **Parallel folding** — rented GPU pods, one per predictor per bucket. 5 seeds × 5 samples minimum.
5. **Scoring pre-submission** — sanity check (pLDDT, inter-model RMSD) to catch broken outputs; no ground truth available yet so we submit best-aggregate.
6. **Submit** — top-5 poses per target, per predictor. CASP limits typically allow multiple registered groups per lab; we'll register one per predictor so each gets evaluated independently.
7. **Log** — per-target cost, wall time, model confidence. Commit to repo for post-hoc analysis.

## Tractable targets

Chai-1's 2048-token cap rules out ~10% of typical CASP targets. Protenix-v1's practical ceiling is similar. Boltz-2's is higher but not unlimited.

We will:
- Submit for every target that fits in the smallest predictor's context.
- Partially cover larger targets using Protenix/Boltz where Chai-1 can't.
- Skip targets exceeding all three (chain-chunking is out of scope for v1).

Realistic coverage: 60–80% of CASP17 targets. Good submissions on ~50%.

## Submission registration

CASP allows registration as a "group" with a short identifier. Plan:
- Register as `zcasp17` (primary).
- Three per-predictor sub-groups for clean attribution: `zcasp17-chai`, `zcasp17-protenix`, `zcasp17-boltz`.
- Observer status for learning: meta-group `zcasp17-meta` for cross-model aggregates.

## Budget

| Phase | Estimate |
|---|---|
| Boltz-2 adapter + retrospective on CASP15/16 | ~$50 compute + ~40 hr engineering |
| CASP17 live submission (all tractable targets, 3 models) | $10–$30 compute |
| Post-results retrospective (re-run any targets with updated weights, cross-model analysis) | $20–$50 compute |
| **Total GPU-cost floor** | **~$80–$130** |

Engineering hours are the real cost. This is where grant money matters — it's 4–6 months of focused solo work, not compute.

## What would "success" look like

Tier 1 — any single tractable target in the top-10 of submitted groups. Mechanically plausible given CASP15 T1124 (#1 / 26 groups) and CASP16 L1001 (#1 / 30 groups) precedent.

Tier 2 — cross-model comparison ships as the first independent CASP17 retrospective. One researcher, one pipeline, three models, published within 30 days.

Tier 3 — the canonical inputs tree + adapters get reused by at least one academic group for their own CASP17 work. The infrastructure compounds.

Even Tier 1 alone would make this the first fully-independent, non-institutional CASP17 submission of substance. Tier 3 is the real long-term payoff.

## Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| GPU availability during submission window | Low | Three providers (RunPod, Vast, Lambda). Pre-book if needed. |
| CASP17 registration bureaucracy (institutional affiliation required?) | Medium | Verify during registration; fallback to observer status if blocked |
| Token-limit targets dominate | Low | Already rare (<10% in CASP15/16) |
| New pathogenic input bugs found mid-competition | Medium | Cross-model comparison catches most; spot-test validation protocol from the homodimer case is in place |
| Weight-cache desync between predictors | Low | Pinned model versions; commits include weight hashes |

## Open questions

- **Affinity prediction.** CASP16's L1000/L3000 had a "PA" task (pose + affinity). Chai-1/Protenix produce poses, not affinities. A DiffDock or GNINA add-on would slot in but is out of scope for v1.
- **Pose ensembling.** Can we pool multiple samples per predictor or ensemble across predictors? Research question rather than infrastructure question — target after the first retrospective.
- **Ligand bond-order bugs.** The OST graph-iso failures (T1158v3, L1006–L1008) could be fixed with proper RDKit bond-order assignment. Low-lift, high-value; on the fix list before CASP17.

## Registration & contact

Registration in progress at [predictioncenter.org/casp17/](https://predictioncenter.org/casp17/). Observer credentials active; group registration pending.

For collaboration, coordination, or shared infrastructure: `sness@sness.net`.
