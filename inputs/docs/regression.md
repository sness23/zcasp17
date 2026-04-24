# Regression report: adapter outputs vs hand-written pipeline inputs

Last run: 2026-04-21.

A diff of every target's adapter-produced input against the hand-written
input in the legacy pipelines (`../casp{15,16}_ligands*/`). Verifies
the `inputs/` canonical tree + adapters reproduce what we've been
actually running, or documents any intentional divergence.

## Results

| Comparison | identical | differ | no-ref-file |
|---|---|---|---|
| Protenix CASP15 | **12 / 12** run targets | 0 | 15 (canonical-only: R*, H1114, H1171/H1172, T1170, T1181, T1186) |
| Protenix CASP16 | 212 / 233 | 21 (Category 2 below) | 0 |
| chai CASP15 | 3 / 25 | 22 (Categories A + B below) | 2 (canonical-only) |
| chai CASP16 | 208 / 233 | 25 (Category 3 below) | 0 |

## Known diff categories

### Category 1 — cosmetic (fixed)

**What:** ligand entity ordering in Protenix CASP15 JSONs; chain ID
string style; ligand name template.

**Was:** adapter emitted `Counter.most_common()` order (K⁺×9 before
Cl⁻×3 on H1135); chain IDs like `chainA_379aa`; ligand names like
`lig-SAH-1`.

**Now:** `build_canonical.py` uses `Counter.items()` (first-seen
source order); chai adapter uses `source_chains` original IDs; ligand
name template `lig{i}-{LABEL}`.

**Status:** fixed. Protenix CASP15 is byte-identical on all 12 benchmark
targets.

### Category 2 — Protenix CASP16 ligand `count` semantics

**What:** for CASP16 pharma dimers (L4000 set, L2001), adapter emits
`{"ligand": {..., "count": 2}}` (two crystal SDFs → 2 copies),
whereas the hand-written JSONs we actually ran had `count: 1`
(one SMILES in the TSV, regardless of crystallographic multiplicity).

**Impact:** At inference time, `count: N` tells Protenix to predict N
ligand copies. `count: 1` predicts one, which OST then matches against
one of the crystal references (leaving the other unmatched). Neither
is strictly "wrong":
- `count: 1` is closer to the original "one SMILES per TSV row" semantics.
- `count: 2` is closer to the crystal's 2-copy stoichiometry.

**Status:** **not fixed** per decision to stick with `count: 1`
semantics for the existing benchmark numbers. Adapter output still
produces `count: N` from the canonical (accurate to the crystal),
so if we want count=1 behavior for a rerun, we need either:

- A flag on the adapter to override `count: 1` for all ligands (e.g.
  `--predict-one-per-ligand`).
- A `predict_count` field in the canonical schema separate from crystal `count`.

Revisit after the full CASP16 post-fix rerun when we have clean numbers.

### Category 3 — chai CASP16 fastas have the homodimer concat bug

**What:** the hand-written chai CASP16 FASTAs at
`../casp16_ligands/fastas/` concatenate homodimer chains into a single
long receptor sequence (e.g. L4001 = 610-aa chimera of chains A+B).
This is the same bug we fixed on the Protenix side — but the chai
pipeline's `build_chai_fasta.py` was never updated.

**Impact:** chai's CASP16 results (at `../casp16_ligands/results_A40/`
and the scoring in `../casp16_ligands/out/lddt_pli_L40S/`) were
generated from buggy input. The catastrophic failures we attributed to
chai on L2001, L4003 (both methods at 40+ Å) may not be chai's fault
at all — chai may well recover if re-run with the adapter's correctly-
split fastas.

**Status:** flagged in `../casp16_ligands/README.md`. Rebuilding chai
CASP16 fastas from the canonical tree is a one-liner:

```bash
python3 adapters/to_chai_fasta.py casp16/*/target.json \
    --out-dir ../casp16_ligands/fastas/
```

Running chai again with the new fastas would update the
head-to-head comparison in
`../casp16_ligands_protenix/docs/02_results.md §2` — currently that
comparison is Protenix-on-clean-input vs chai-on-buggy-input, which
understates chai's actual capability.

### Category A — chai multi-chain interleave order

**What:** some CASP15 chai FASTAs (H1135, H1171, H1172 variants) have
chains emitted in source PDB order (ABCDEFG), interleaving chains of
different sequences. The adapter groups chains by unique sequence
(all of sequence-1, then all of sequence-2).

**Impact:** chai doesn't parse FASTA record order semantically, so
this is purely cosmetic. Semantically equivalent.

**Status:** unfixed. Low-value cosmetic match would require the chai
adapter to re-emit in source-chain order, which contradicts the
"group by unique sequence" design. Accepting the difference.

### Category B — chai RNA chain naming

**What:** CASP15 RNA targets (R\*) use chain name `0` in the source
PDB. Legacy chai builder emits `>rna|name=<T>-chain0`. Adapter's
single-chain branch emits `>rna|name=<T>-receptor`.

**Impact:** cosmetic. Chai ignores the label text.

**Status:** unfixed. One-line change possible if needed — add
`source_chains[0]` to the single-chain header template in
`to_chai_fasta.py`.

## Untested coverage

The regression runs on **all 260 canonical targets**, but only ~55 of
them (12 CASP15 + 44 CASP16 PoC) have corresponding hand-written
inputs to compare against. The remaining 205 (mostly L3000 Autotaxin
and CASP15 excluded targets) are generated-only; no ground truth for
comparison.

When we extend benchmark runs to those targets (L3000 batch, CASP15
full set), the adapter output is the authoritative source — there's
no hand-written version to regress against.

## Running the regression check

```bash
cd ~/data/vaults/casp/inputs

# regenerate adapter outputs to scratch
rm -rf /tmp/_adapter_check && mkdir -p /tmp/_adapter_check/{protenix15,protenix16,chai15,chai16}
python3 adapters/to_protenix_json.py casp15/*/target.json --out-dir /tmp/_adapter_check/protenix15
python3 adapters/to_protenix_json.py casp16/*/target.json --out-dir /tmp/_adapter_check/protenix16
python3 adapters/to_chai_fasta.py    casp15/*/target.json --out-dir /tmp/_adapter_check/chai15
python3 adapters/to_chai_fasta.py    casp16/*/target.json --out-dir /tmp/_adapter_check/chai16

# diff against hand-written
for y in 15 16; do
    for tool in protenix chai; do
        ref_tool_dir=$([ "$tool" = "protenix" ] && echo "casp${y}_ligands_protenix/jsons" || echo "casp${y}_ligands/fastas")
        ref="$HOME/data/vaults/casp/$ref_tool_dir"
        gen="/tmp/_adapter_check/${tool}${y}"
        n_same=$(diff -rq "$ref" "$gen" 2>/dev/null | grep -c "Only in $gen" || true)
        diff -rq "$ref" "$gen" 2>&1 | head
    done
done
```
