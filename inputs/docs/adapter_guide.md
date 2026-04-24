# Writing a new pipeline adapter

An adapter converts `target.json` (canonical) into whatever input
format a specific structure-prediction pipeline expects.

## Existing adapters

| File | Target pipeline | Output |
|---|---|---|
| `adapters/to_protenix_json.py` | Protenix v1.0.0+ | `<T>.json` (Protenix native JSON) |
| `adapters/to_chai_fasta.py` | chai-lab 0.6.x+ | `<T>.fasta` (chai FASTA with `>protein`, `>ligand` blocks) |
| `adapters/to_dynamicbind.py` | DynamicBind | per-target dir with `protein.fasta`, `ligand.sdf`, etc. |

## Adapter contract

Every adapter must:

1. **Read one or more `target.json` paths** as positional arguments.
2. **Write outputs to an `--out-dir`** (creating it if absent).
3. **Never modify the canonical tree** — adapters are pure functions.
4. **Log `wrote <path>` per target** so the user can verify.
5. **Handle failures gracefully** — print `[error] <src>: <msg>` and
   continue with the remaining targets.

The boilerplate looks like this:

```python
#!/usr/bin/env python3
"""Adapter: canonical target.json → <pipeline X> input format."""
import argparse, json, sys
from pathlib import Path


def convert(target_json_path: Path, out_root: Path) -> None:
    d = json.loads(target_json_path.read_text())
    # ... transform d into pipeline X's format ...
    # ... write to out_root / <something> ...
    pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target_jsons", nargs="+", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for src in args.target_jsons:
        try:
            convert(src, args.out_dir)
            ok += 1
        except Exception as e:
            print(f"[error] {src}: {e}", file=sys.stderr)
    print(f"\n{ok}/{len(args.target_jsons)} converted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

## Mapping canonical → pipeline concepts

The canonical schema exposes these primitives:

| Canonical concept | Common downstream mapping |
|---|---|
| `chains[*]` with `type: "protein"` | protein-chain input, FASTA record, input entity |
| `chains[*]` with `type: "rna"` / `"dna"` | NA input (if pipeline supports; otherwise skip or error) |
| `chains[*].count > 1` | repeat the chain N times, or use a native "count" / "copies" field |
| `ligands[*]` with `ccd` populated | use pipeline's CCD-by-code input if available (`CCD_SAH`, `ccd:SAH`, etc.) |
| `ligands[*]` with `smiles` only | SMILES input |
| `ligands[*].count > 1` | repeat the ligand N times, or use native copy-count |
| `receptor_pdb` | usually NOT fed to the predictor (it's a ground-truth file), but some docking tools consume it |

## Common gotchas

### 1. Chain counts

Different pipelines handle `count` differently:

- **Protenix:** native `{"count": N}` field in the JSON. Adapter just copies.
- **chai-lab:** no count; emit N separate `>protein|...` FASTA entries with
  numeric suffix IDs.
- **DynamicBind:** single-chain assumption; adapter needs to decide
  whether to concatenate (discouraged) or error out on multi-chain input.

### 2. Ligand count

Same pattern as chains. Critically, **don't swap count for duplicating
the SMILES in one field** — if the pipeline wants "this ligand, twice,"
encode that correctly so the predictor samples 2 copies.

### 3. CCD vs SMILES

Many pipelines don't accept CCD codes directly. If the canonical has
`ccd: "SAH"` but no `smiles`, the adapter has two options:

a. **Look up SMILES at adapter time** — extra work but keeps the
   pipeline input self-contained.
b. **Fail gracefully** with a clear error — the build_canonical should
   have populated SMILES; if it didn't, the adapter should refuse.

Our adapters prefer (b): they treat missing SMILES as a
canonical-tree bug, not an adapter-time problem.

### 4. NA chains

Some pipelines don't support nucleic acids. Adapter should:

- Error out with a clear message if ANY `chains[*].type != "protein"`,
  or
- Skip NA entities and log a warning, or
- Silently drop them (discouraged — silent data loss).

Pick one and document in the adapter's docstring.

## Running adapters at scale

Batch-convert all CASP15+16 targets in one go:

```bash
# Protenix — CASP15 and CASP16 into their respective pipeline dirs
python3 adapters/to_protenix_json.py \
    casp15/*/target.json \
    --out-dir ../casp15_ligands_protenix/jsons/

python3 adapters/to_protenix_json.py \
    casp16/*/target.json \
    --out-dir ../casp16_ligands_protenix/jsons/

# chai — same target → different pipeline dir
python3 adapters/to_chai_fasta.py \
    casp15/*/target.json \
    --out-dir ../casp15_ligands/fastas/

python3 adapters/to_chai_fasta.py \
    casp16/*/target.json \
    --out-dir ../casp16_ligands/fastas/
```

## Adding a new adapter

Template checklist:

- [ ] Create `adapters/to_<pipeline>.py`.
- [ ] Read input from `target.json`, write output to `--out-dir`.
- [ ] Handle the count / CCD / multi-chain gotchas above.
- [ ] Add an entry to the table in this doc and in `inputs/README.md`.
- [ ] Smoke-test on one CASP15 + one CASP16 target before batch use.

## When the canonical schema isn't enough

If a pipeline needs something the canonical schema doesn't express
(e.g. specific residue modifications, covalent bond annotations, or
pocket-restraint hints), there are three options:

1. **Extend the schema** — add an optional field, update
   `schema.md`, update `build_canonical.py` to populate it, update
   all adapters that might use it.
2. **Add a pipeline-specific sidecar file** — e.g.
   `inputs/casp15/T1124/protenix_restraints.json` alongside
   `target.json`. Adapters for that pipeline read both.
3. **Encode the extra data in the `notes` field** — okay for
   one-off annotations, discouraged for structured data.

Option 1 is preferred when the field is useful to multiple pipelines;
option 2 when the data is truly pipeline-specific.
