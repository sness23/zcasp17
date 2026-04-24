#!/usr/bin/env python3
"""Adapter: canonical target.json → Protenix input JSON.

Protenix expects:
    [{
      "name": "<target>",
      "sequences": [
        {"proteinChain": {"sequence": "...", "count": N}},
        {"rnaSequence":   {"sequence": "...", "count": N}},
        {"dnaSequence":   {"sequence": "...", "count": N}},
        {"ligand":        {"ligand": "CCD_<code>" or "<SMILES>", "count": N}},
        {"ion":           {"ion": "MG", "count": N}}
      ]
    }]

Rules applied by this adapter:
- Every `chains[*]` in target.json becomes a matching entity in the output.
- Ligands: use `CCD_<code>` if `ccd` is set, otherwise the raw SMILES string.
- Ions (Mg, K, etc.) are currently NOT auto-promoted — they appear as
  regular ligands. If we need true ion entities, flag in conventions.md.

Usage:
    python3 to_protenix_json.py <target.json>... --out-dir DIR
    python3 to_protenix_json.py ../casp15/*/target.json --out-dir /workspace/casp15_ligands_protenix/jsons
"""
import argparse
import json
import sys
from pathlib import Path


def convert(target_json_path: Path) -> list:
    d = json.loads(target_json_path.read_text())
    sequences = []

    for chain in d["chains"]:
        kind = chain["type"]
        entry_key = {
            "protein": "proteinChain",
            "rna": "rnaSequence",
            "dna": "dnaSequence",
        }[kind]
        sequences.append({
            entry_key: {"sequence": chain["sequence"], "count": chain["count"]}
        })

    for lig in d.get("ligands", []):
        if lig.get("ccd"):
            lig_str = f"CCD_{lig['ccd']}"
        elif lig.get("smiles"):
            lig_str = lig["smiles"]
        else:
            print(f"  [warn] {d['target_id']} ligand {lig.get('id','?')}: no ccd or smiles — skipped",
                  file=sys.stderr)
            continue
        sequences.append({"ligand": {"ligand": lig_str, "count": lig["count"]}})

    return [{"name": d["target_id"], "sequences": sequences}]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target_jsons", nargs="+", type=Path,
                    help="one or more target.json paths")
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="output directory (will be created)")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for src in args.target_jsons:
        try:
            out = convert(src)
            tid = out[0]["name"]
            dst = args.out_dir / f"{tid}.json"
            dst.write_text(json.dumps(out, indent=2) + "\n")
            print(f"wrote {dst}")
            ok += 1
        except Exception as e:
            print(f"[error] {src}: {e}", file=sys.stderr)
    print(f"\n{ok}/{len(args.target_jsons)} converted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
