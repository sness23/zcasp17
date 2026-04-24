# 03 — Discussion

## 1. The homodimer bug

The most consequential finding of this benchmark was a bug in our
input-preparation pipeline, not a property of Protenix itself.

### 1.1 What happened

Our initial `build_protenix_json.py` for CASP16 treated the receptor
as a single long protein chain by concatenating all protein residues
from `protein_aligned.pdb` in file order. For monomeric receptors
(Chymase / L1000) this is fine — one chain in, one chain out. For
homodimeric receptors (Cathepsin G / L2000 and all Mpro / L4000
targets) this produced a **chimeric fake-monomer** input where chain
A's sequence was immediately followed by chain B's sequence with no
interruption.

Example: L2001 Cathepsin G has chain A (226 aa) and chain B (227 aa,
differing by a single residue at one terminus). Our concatenated
representation gave Protenix a 453-aa protein that doesn't exist in
nature. The model dutifully folded this chimera, produced a passable
protein structure, and placed the ligand somewhere — but that
"somewhere" was not at the inter-chain binding site of the real
dimer. Result: LDDT_pli 0.00, RMSD 50 Å.

### 1.2 How we found it

Three observations converged:

1. **Shared failures with chai.** L2001 and L4003 both failed under
   chai-lab with similarly large RMSDs (~40 Å). A method-specific bug
   wouldn't show up in both pipelines, so our first hypothesis was
   "hard target with symmetric binding site."
2. **High confidence on wrong answer.** Protenix's `ranking_score`
   was 0.82 for L2001 and 0.88 for L4001 — high, despite LDDT_pli = 0.
   This is characteristic of a model that's confident about the
   *protein structure* but the ligand is drawn from a plausible-looking
   but wrong distribution. Inconsistent with "the model is genuinely
   uncertain."
3. **Inventory of multi-chain targets.** Audit showed 25 of 44 PoC
   targets are multi-chain (all L4000 + both L2000). Their affected
   rate tracked perfectly with observed catastrophic failures.

Inspection of `jsons/L2001.json` showed the concatenated 453-aa sequence.

### 1.3 The fix

`build_protenix_json.py` now iterates per-chain, deduplicates by
sequence identity, and emits:

- `{"proteinChain": {..., "count": N}}` if ≥2 chains have identical
  sequences (e.g. L4009).
- Multiple `proteinChain` entries with `count: 1` each if sequences
  differ (common — truncation variants are everywhere in these PDBs).

### 1.4 Spot-test validation

| Target | Bug status | Pre-fix | Post-fix | Outcome |
|---|---|---|---|---|
| L2001 | Dimer, 1-aa-different chains | 0.00 / 50 Å | **0.75 / 1.9 Å** | Full rescue |
| L4001 | Dimer, 2-aa-different chains | 0.01 / 47 Å | **0.78 / 1.8 Å** | Full rescue |
| L4003 | Dimer, 5-aa-different chains | 0.00 / 44 Å | 0.46 / 5.2 Å | Partial rescue |
| L4004 | Dimer, 2-aa-different chains | 0.36 / 5.6 Å | 0.36 / 5.4 Å | No change — bug wasn't the cause |

L4004 confirms the fix isn't a universal improvement — it specifically
rescues *catastrophic* failures caused by the chimera input. Targets
that were merely "poor" for other reasons stay poor.

### 1.5 What this teaches

- **Multi-chain input validation needs explicit testing.** This bug
  slipped past us because L1000 (monomer) results looked great, giving
  false confidence that the builder was correct.
- **Shared failures across models don't mean shared causes.** Our L2001
  catastrophe looked identical to chai's, but they might well have
  different root causes (chai's MSA generation may or may not have
  handled the real dimer correctly).
- **Confident wrong answers are a strong signal.** When `ranking_score`
  is high but the predicted structure violates basic binding-site
  geometry, suspect the input.

## 2. Where Protenix wins vs chai

### 2.1 L4002 (Mpro) — 34 Å swing in Protenix's favor

Both methods ran with comparable MSAs (Mpro is well-represented).
Protenix placed the ligand at 0.45 Å; chai at 34 Å. This target has
one clear binding site and a well-defined ligand; why chai missed is
unclear from this benchmark alone. Could be initialization variance,
could be a confounding crystal artifact, could be chai's attention
over the wrong residues. Worth a closer look when publishing.

### 2.2 Consistent on Chymase (L1003–L1005)

Three Chymase targets (L1003, L1004, L1005) show LDDT_pli ≥ 0.97 with
RMSDs < 0.7 Å — an unambiguous win for Protenix on this class. The
Chymase receptor is a classic serine protease fold with deep PDB
coverage; templates likely contribute meaningfully here. Protenix's
template support (on by default with `--use_template true`) versus
chai's template-off default is plausibly what's driving these wins.

## 3. Where chai wins

### 3.1 L4001 (pre-fix) — chai 2 Å, Protenix 47 Å

This was entirely the homodimer bug. Post-fix Protenix is at 1.78 Å,
recovering to a narrow Protenix win (0.78 vs chai's 0.72 LDDT).

### 3.2 L1001, L1002 — chai slightly ahead

Chai's LDDT_pli leads by 0.05 (L1001: 0.96 vs 0.91) and 0.13 (L1002:
0.81 vs 0.68). Differences in this range could be seed noise — chai
ran 5 seeds total, Protenix 5 seeds × 5 samples = 25 structures.
Best-of-25 should statistically beat best-of-5, so Protenix should
generally win unless chai's single seed happened to draw a particularly
good sample. Not a reliable-enough gap to claim a real method-level
advantage for chai on these targets.

## 4. Ligand-class and target-class observations

### 4.1 Drug-like novel ligands (the whole CASP16 premise)

Unlike CASP15's mix of common cofactors (SAH, ATP, NAG), **every**
CASP16 pharma-ligand target is a novel drug candidate with no CCD
assignment. This has two effects:

- **No cache hit** in Protenix's CCD database; each ligand is built
  from scratch.
- **Higher OST scorer failure rate** on graph-isomorphism matching.
  L1006/1007/1008 failed in both pipelines because of bond-perception
  differences between Protenix's internal representation and the ref
  SDF's RDKit-canonical graph.

### 4.2 Mpro inhibitors

Mpro (L4000, 25 targets, all homodimer) has a well-defined active site
and many X-ray structures in the PDB for context. Post-fix results on
L4001/L4003/L4004 hint at a bimodal outcome: some drugs are placed
well (L4002: 0.90), some land in the right pocket but wrong pose
(L4004 at 5.4 Å), and a small fraction remain truly hard. We should
see more of this pattern once the full 25-target re-run completes.

## 5. Pipeline lessons

### 5.1 Concatenated-homodimer input is a sharp-edge

This bug is embarrassingly easy to make. Any pipeline that "collapses
PDB chains into a sequence" for a multi-chain receptor is vulnerable.
Sanity check: count protein chains in the source PDB; if > 1, emit
multiple `proteinChain` entries.

### 5.2 `ranking_score` on multi-ligand complexes

See [02_results.md §5](02_results.md#5-confidence-calibration). For
multi-ligand targets (CASP16 is always 1 receptor + 1 ligand in the
PoC, but Mpro crystal has the ligand at a clear site in a dimer
interface), `ranking_score` is a protein-dominated metric — it doesn't
reliably track ligand-pose quality. For per-ligand confidence, prefer
`chain_plddt[<ligand chain idx>]` from the summary-confidence JSON.

### 5.3 Spot-test before full re-runs

Three ~10-min spot tests on RTX 3090/4090 (total cost $0.43) confirmed
the homodimer bug fix before committing to a ~$2, 2-hour full rerun.
Always validate a targeted fix on 2-3 affected targets first; the cost
ratio is usually >10× in your favor.

### 5.4 MSA-server queue dominates local prep

For the initial prep of all 44 targets, server queue time was ~3.5 hr
(occasional 10–30 min waits on uncached sequences) vs ~10 min of
actual local hmmer/kalign CPU work. Self-hosting the MSA server
(see [../docs/self_hosted_mmseqs.md](#) — not yet written) would cut
this but at substantial infra overhead (~$800/mo for a 128 GB
r6i.4xlarge, only pays off at hundreds of queries/month).

For now, running prep locally overnight is the right trade-off — the
MSA server's queue doesn't cost us money the way a rented GPU sitting
idle would.
