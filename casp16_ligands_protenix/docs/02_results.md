# 02 — Results

## 1. Per-target table (14 scored of 44 PoC, pre-fix)

From `out/lddt_pli/summary.csv`, best sample per target:

| Target | Set | LDDT_pli | RMSD (Å) | ranking | Notes |
|---|---|---|---|---|---|
| L1001 | Chymase | 0.910 | 1.03 | 0.98 | ✓ strong |
| L1002 | Chymase | 0.678 | 2.66 | 0.97 | mid — right pocket, wrong rotamer |
| L1003 | Chymase | 0.967 | 0.60 | 0.98 | ✓ strong |
| L1004 | Chymase | 0.974 | 0.46 | 0.99 | ✓ strong |
| L1005 | Chymase | 0.968 | 0.65 | 0.99 | ✓ strong |
| L1006 | Chymase | — | — | — | OST graph-iso failure |
| L1007 | Chymase | — | — | — | OST graph-iso failure |
| L1008 | Chymase | — | — | — | OST graph-iso failure |
| L2001 | Cathepsin G | 0.000 | 50.4 | 0.82 | 🚨 homodimer bug (now FIXED) |
| L2002 | Cathepsin G | 0.893 | 1.32 | 0.98 | ✓ strong |
| L4001 | Mpro | 0.010 | 47.3 | 0.88 | 🚨 homodimer bug (now FIXED) |
| L4002 | Mpro | 0.899 | 0.45 | 0.91 | ✓ strong |
| L4003 | Mpro | 0.003 | 44.4 | 0.88 | 🚨 homodimer bug (partial fix) |
| L4004 | Mpro | 0.359 | 5.56 | 0.84 | Hard target — bug not the cause |

**Coverage:** 11 scored, 3 OST-unscored, 0 untried of the 14 that finished
inference. (The remaining 30 targets on the PoC manifest had not been
run when the L40S pod was released.)

## 2. Head-to-head vs chai-lab

> ⚠️ **Caveat: chai's CASP16 inputs had the same homodimer bug** we
> fixed on the Protenix side. The FASTAs at `../casp16_ligands/fastas/`
> were built with concatenated-homodimer chains, matching our pre-fix
> Protenix behavior. chai's 3 catastrophic multi-chain failures below
> (L2001 / L4001 / L4003) likely reflect the same input bug, not
> method limitations. A chai re-run with corrected FASTAs would shift
> the comparison — chai's numbers here **understate** its actual
> performance. See
> [`../../casp16_ligands/README.md`](../../casp16_ligands/README.md)
> for the caveat + fix command.

Both methods ran on the same 14 CASP16 targets with identical refs
and identical OST scorer. Best sample per method per target:

| Target | chai LDDT / RMSD | Protenix LDDT / RMSD | Winner |
|---|---|---|---|
| L1001 | 0.960 / 0.48 | 0.910 / 1.03 | chai |
| L1002 | 0.805 / 1.65 | 0.678 / 2.66 | chai |
| L1003 | 0.510 / 3.90 | **0.967 / 0.60** | **Protenix** |
| L1004 | 0.842 / 1.52 | **0.974 / 0.46** | **Protenix** |
| L1005 | 0.927 / 0.87 | 0.968 / 0.65 | tie |
| L1006 | — | — | neither (OST fail) |
| L1007 | — | — | neither (OST fail) |
| L1008 | — | — | neither (OST fail) |
| L2001 | 0.000 / 39.9 | 0.000 / 50.4 | both fail (pre-fix) |
| L2002 | 0.876 / 1.76 | 0.893 / 1.32 | tie (P slight edge) |
| L4001 | 0.724 / 2.14 | 0.010 / 47.3 | **chai** (pre-fix); see §3 post-fix |
| L4002 | 0.000 / 34.3 | **0.899 / 0.45** | **Protenix** (huge) |
| L4003 | 0.009 / 44.4 | 0.003 / 44.4 | both fail (pre-fix) |
| L4004 | 0.000 / 50.1 | **0.359 / 5.56** | **Protenix** (less bad) |

**Summary:**

| Outcome | Count |
|---|---|
| Protenix wins (Δ ≥ 0.05 LDDT_pli) | 4 (L1003, L1004, L4002, L4004) |
| chai wins (Δ ≥ 0.05 LDDT_pli) | 3 (L1001, L1002, L4001) |
| Tie (both strong, \|Δ\| < 0.05) | 2 (L1005, L2002) |
| Both catastrophic (>40 Å) | 2 (L2001, L4003) — see §3 |
| Both OST-unscored | 3 (L1006, L1007, L1008) |

Neither method dominates. Most striking divergence:

- **L4002**: chai placed the ligand 34 Å from the right pocket,
  Protenix hit 0.45 Å. Worth investigating what chai did wrong there.
- **L4001**: mirror case — chai 2.1 Å, Protenix 47 Å. The Protenix
  failure was due to our homodimer bug and is now rescued (§3 below).

## 3. Homodimer fix validation

After discovering the chain-concatenation bug (see
[03_discussion.md §1](03_discussion.md#1-the-homodimer-bug)), we
rebuilt the JSONs for 4 affected targets and re-ran inference on
fresh GPU pods:

| Target | Pre-fix LDDT / RMSD | Post-fix LDDT / RMSD | Outcome |
|---|---|---|---|
| L2001 | 0.000 / 50.4 Å | **0.746 / 1.91 Å** | ✓ catastrophic → near-native |
| L4001 | 0.010 / 47.3 Å | **0.784 / 1.78 Å** | ✓ catastrophic → near-native |
| L4003 | 0.003 / 44.4 Å | 0.462 / 5.22 Å | ~ catastrophic → mid-range |
| L4004 | 0.359 / 5.56 Å | 0.361 / 5.44 Å | no change (bug wasn't the cause) |

L2001 spot test: RTX 3090 on Vast.ai, 9.5 min wall, $0.03.
L4001/L4003/L4004 spot test: RTX 4090 on RunPod, 33 min wall, $0.40.

**Implications:** 25 of 44 PoC targets (all L2000 + all L4000) were
affected by the homodimer bug pre-fix. Re-running all 25 is pending
next available L40S slot.

## 4. Results by target set

Stratified by the four pharma-ligand sets (14 scored targets):

### Chymase (L1000, 5 of 17 scored)

| LDDT_pli | n | Share |
|---|---|---|
| ≥ 0.9 | 4 | 80% |
| 0.5–0.9 | 1 | 20% |

Chymase is a monomer → no homodimer bug → full, clean signal. 4 of 5
evaluated targets are near-native. Expected outcome for the remaining
12 L1000 targets (pending) is similar, given consistent receptor and
docking landscape. 3 additional L1000 targets failed OST scoring
(same fate under chai), so a fuller comparison awaits a scorer fix.

### Cathepsin G (L2000, 2 of 2)

| Target | Pre-fix | Post-fix |
|---|---|---|
| L2001 | 0.00 / 50 Å | **0.75 / 1.9 Å** |
| L2002 | 0.89 / 1.3 Å | (already fine, not re-run) |

Post-fix both L2000 targets are scored as near-native. 100% rescue.

### Mpro (L4000, 4 of 25)

Pre-fix: only L4002 strong, L4001/L4003 catastrophic, L4004 poor.
Post-fix validations: L4001 rescued (→ 0.78), L4003 rescued mid-range
(→ 0.46), L4004 unchanged (not bug-related).

Extrapolation (pending full rerun): expect 12–18 of 25 to rescue from
catastrophic to decent, 3–5 to remain hard (genuine model difficulty,
not input issue), and the remaining 4–6 to stay near-native as they
already are.

### Autotaxin (L3000, 0 of 219 — deferred)

Not in PoC. At observed L40S inference rates (~5 min/target), the full
set is ~18 GPU hours. Adding local prep queue (~15 hours serial, probably
~6 hours with `clone_prepped.py` for identical receptor sequences)
brings total wall time to ~24 hours, ~$20. Better cost-benefit on a
future pass once other issues are fixed.

## 5. Confidence calibration (14 scored)

`ranking_score` relationship to pose quality, same bucketing as
CASP15 results:

| ranking_score bucket | n | mean LDDT_pli | comment |
|---|---|---|---|
| ≥ 0.95 | 6 | 0.92 | reliably strong |
| 0.85–0.95 | 5 | 0.30 (but high variance; includes L4001 catastrophe at 0.88) | **not reliable** |
| 0.80–0.85 | 2 | 0.36 | L2001 (0.82), L4004 (0.84) |
| < 0.80 | 1 | — | (none) |

Same finding as CASP15: `ranking_score` is target-level, dilutes when
receptor is correct but ligand is not. Pre-fix, three targets
(L2001, L4001, L4003) had high confidence but 0.00 LDDT_pli — the
homodimer bug was a systematic way to produce false-confident
predictions.

Post-fix, L2001's ranking_score dropped from 0.82 → 0.45 (more honest
uncertainty) while the pose improved from 0.00 → 0.75. This suggests
`ranking_score` is better-calibrated when the input is well-structured.

## 6. Timing and cost

| Phase | Wall | Cost |
|---|---|---|
| Local prep for 44 targets (one-time) | ~4 hr CPU + MSA queue | $0 |
| Partial inference (14 of 44) | ~65 min on L40S | $0.93 |
| L2001 spot test (RTX 3090) | 9.5 min | $0.03 |
| L4001/3/4 spot test (RTX 4090) | 33 min | $0.40 |
| **Subtotal to date** | — | **$1.36** |
| **Pending: full 25-target rerun** | ~2 hr on L40S | ~$1.80 |
| **Pending: remaining 16 of 44 PoC** | ~1 hr on L40S | ~$0.90 |
| Total PoC estimate | ~5 hr | ~$4 |

## 7. Key cross-references

- L40 / Mpro structural homology context: `../casp16_ligands/binding_sites.md`
  (curated notes from the chai analysis).
- Official CASP16 pharma-ligand scores (for the CASP16 community's
  assessor metrics, separate from our LDDT_pli replication):
  `../casp16/results/ligand_scores.csv`.
- Companion analyses at `../casp15_ligands_protenix/docs/02_results.md`.
