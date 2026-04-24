# The homodimer bug — case study

Discovered 2026-04-20 during the CASP16 pharma-ligand benchmark. A silent failure mode in our (and chai-lab's) CASP16 input builders that caused catastrophic ligand misplacement — 30–50 Å off, entirely outside the real binding pocket — on 25 of 44 pharma-ligand targets.

This is a writeup of what we found, how we found it, and how it was fixed.

## TL;DR

**Bug:** Our CASP16 input builder concatenated homodimer chains from `protein_aligned.pdb` into a single "fake-monomer" sequence. The model then folded a chimeric protein that doesn't exist in nature and placed the ligand in a non-pocket location.

**Blast radius:** 25 of 44 PoC targets (Cathepsin G + all Mpro). Same bug exists unfixed in the canonical chai-lab CASP16 pipeline.

**Fix:** Rewrite the builder to iterate per-chain with exact-identity deduplication. Near-identical chains (truncation variants) become separate `proteinChain` entries. Only truly identical chains get `count: 2`.

**Validation:** 4 spot-tests. Catastrophic targets jumped from LDDT_pli 0.00 → 0.75/0.78/0.46.

## The symptom

The CASP16 Protenix run started producing obviously-wrong outputs on multi-chain targets:

| Target | `ranking_score` (model's own confidence) | LDDT_pli | BiSyRMSD |
|---|---|---|---|
| L2001 (Cathepsin G) | 0.82 | 0.00 | 50 Å |
| L4001 (Mpro) | 0.87 | 0.01 | 47 Å |
| L4003 (Mpro) | 0.84 | 0.00 | 44 Å |

The model was confident. The pose was catastrophically wrong. **This is the worst kind of failure** — the same signal you'd see on a genuinely hard target (no confident hit) does not appear, so you can't gate on `ranking_score` alone.

## How it was found

Three observations converged:

1. **L2001 and L4003 failed catastrophically in both Protenix and Chai-1.** Two independent predictors with independent training sets producing the same catastrophic failure on the same targets is not a model bug. It's an input bug.

2. **High `ranking_score` on obviously-wrong placements.** If the model were uncertain, you'd expect low confidence. The model was confidently wrong. That pointed at the inputs being well-formed but biologically nonsensical — the classic shape of a data-pipeline bug, not a modeling bug.

3. **Audit: 25 of 44 PoC targets are multi-chain.** Running a quick chain-count script against all 44 inputs revealed the issue structurally: most of the catastrophes were dimers, and our builder was mishandling multi-chain receptors.

The third observation closed the case. A multi-hour inspection of `build_protenix_json.py` revealed the concatenation.

## The bug, in detail

Simplified pre-fix logic:

```python
# WRONG: concatenate all chains
sequence = ""
for chain in pdb.get_chains():
    sequence += extract_sequence(chain)

entity = {
    "proteinChain": {
        "sequence": sequence,
        "count": 1,
    }
}
```

For Chymase (L1000 — single chain), this happened to be correct.

For Cathepsin G (L2000 ×2) and Mpro (L4000 ×25), this concatenated two receptor chains into a single artificial sequence `chain_A + chain_B`. The model was being asked to fold a protein that doesn't exist in nature.

What Protenix-v1 did with that request was reasonable: it folded a structurally-coherent chimera, and placed the ligand somewhere in that chimera's surface — which has no geometric relationship to the actual binding pocket in the real dimer.

## The fix

```python
# RIGHT: per-chain with exact-identity deduplication
chain_seqs = {chain.id: extract_sequence(chain) for chain in pdb.get_chains()}

# Group chains by sequence identity
groups = {}
for cid, seq in chain_seqs.items():
    groups.setdefault(seq, []).append(cid)

# Emit one proteinChain per distinct sequence
entities = []
for seq, members in groups.items():
    entities.append({
        "proteinChain": {
            "sequence": seq,
            "count": len(members),
        }
    })
```

Key design decision: **exact-identity deduplication, not alignment-based similarity**.

Homodimers (two identical chains A/B) → one `proteinChain` with `count: 2`. Near-identical chains (truncation variants where chain B has N-terminal residues 1–5 missing) → two separate `proteinChain` entries with `count: 1` each. The alternative (alignment-based merge with a 99% identity threshold) would lose the truncation information that a few CASP16 targets actually encode.

## Validation

Ran four targets through the post-fix pipeline on rented GPUs (Vast.ai 3090 + RunPod 4090):

| Target | Pre-fix LDDT_pli | Post-fix LDDT_pli | Pre-fix RMSD | Post-fix RMSD | Verdict |
|---|---|---|---|---|---|
| L2001 (Cathepsin G + benzoxazinone) | 0.00 | **0.75** | 50 Å | **1.9 Å** | Rescued |
| L4001 (Mpro + covalent inhibitor) | 0.01 | **0.78** | 47 Å | **1.8 Å** | Rescued |
| L4003 (Mpro + non-covalent inhibitor) | 0.00 | 0.46 | 44 Å | 5.2 Å | Partial rescue |
| L4004 (Mpro + inhibitor) | 0.36 | 0.36 | 5.6 Å | 5.4 Å | No change (not bug-related) |

L4004's failure was unrelated to the homodimer bug — probably a genuine prediction miss on a harder Mpro pose. The fact that it didn't change post-fix is actually good validation that the fix is targeted.

Validation cost: $0.43 on spot GPUs.

## The broader lesson

The same bug exists unfixed in the chai-lab CASP16 FASTA builder. **Any CASP16 submission that used the default chai-lab pipeline on multi-chain receptors probably suffered the same failure mode.** We've flagged this in the chai-lab sibling pipeline with a one-line fix that regenerates FASTAs from the canonical inputs tree.

More broadly: this is an argument for running multiple independent predictors against the same inputs. The cross-model signal — "both pipelines fail catastrophically on the same targets" — was what made the bug discoverable. A single-model benchmark would have let this slide.

This is the thesis of the zcasp17 pipeline: **many models, one evidence framework, cross-model comparison as a primary debugging tool.**

## Process memory

A feedback entry was added to the project memory (`feedback_wipe_stale_prep_after_builder_fix.md`) to flag a related trap: `run_batch.sh` has a skip-prep-if-exists optimization that silently reused pre-fix `jsons_prepped/` files after the builder was updated, yielding pre-fix-equivalent output even on the "fixed" run.

Rule: **after any input-builder change, wipe the prep cache.** The cost of regenerating is minutes; the cost of shipping a contaminated run is the run itself.

## Related

- Pre-fix and post-fix input JSONs are checked into `casp16_ligands_protenix/jsons_prepped/` and `casp16_ligands_protenix/jsons_postfix/`.
- Full audit in `casp16_ligands_protenix/docs/03_discussion.md §1`.
- CASP16 LDDT_pli distribution before and after: `casp16_ligands_protenix/out/lddt_pli_postfix_full/summary.csv`.
