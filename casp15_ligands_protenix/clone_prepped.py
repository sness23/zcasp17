#!/usr/bin/env python3
"""Clone a prepped JSON's MSA/template paths to any sibling target with the
same protein sequence.

Protenix's MSA + template search is keyed to the protein sequence alone —
ligand SMILES doesn't affect it.  So if target A has been prepped and
target B has the same protein sequence, we can produce a valid prepped
JSON for B by copying A's sequences[0] (which carries the MSA paths) and
substituting B's ligand.

This is a cheap sidecar: safe to run while prep_local.sh is mid-batch;
just drops new files into jsons_prepped/ which the main loop will then
SKIP.

Usage:
    python3 clone_prepped.py                        # clone everything possible
    python3 clone_prepped.py --manifest manifest_L1L2L4.txt
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JSONS = ROOT / "jsons"
PREPPED = ROOT / "jsons_prepped"


def load_target(path: Path) -> dict:
    return json.load(path.open())[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", help="restrict to targets in this manifest")
    args = ap.parse_args()

    if args.manifest:
        with open(args.manifest) as f:
            wanted = {
                line.strip() for line in f
                if line.strip() and not line.startswith("#")
            }
    else:
        wanted = None

    # Index already-prepped targets by their protein sequence
    seq_to_prepped: dict[str, Path] = {}
    for p in sorted(PREPPED.glob("*.json")):
        try:
            d = load_target(p)
            seq = d["sequences"][0]["proteinChain"]["sequence"]
            seq_to_prepped.setdefault(seq, p)
        except (KeyError, json.JSONDecodeError, IndexError):
            continue

    if not seq_to_prepped:
        print("no prepped JSONs yet — nothing to clone from")
        return 0

    print(f"found {len(seq_to_prepped)} unique prepped sequence(s)")

    cloned = skipped = unmatched = 0
    for src in sorted(JSONS.glob("*.json")):
        name = src.stem
        if wanted and name not in wanted:
            continue
        dest = PREPPED / f"{name}.json"
        if dest.exists():
            skipped += 1
            continue
        try:
            d = load_target(src)
            seq = d["sequences"][0]["proteinChain"]["sequence"]
        except (KeyError, json.JSONDecodeError, IndexError):
            continue
        donor_path = seq_to_prepped.get(seq)
        if donor_path is None:
            unmatched += 1
            continue

        donor = load_target(donor_path)
        out = [{
            "name": name,
            "sequences": [
                # clone protein chain (keeps MSA + template path fields)
                json.loads(json.dumps(donor["sequences"][0])),
                # use this target's own ligand (and any other non-protein entries)
                *d["sequences"][1:],
            ],
        }]
        if "covalent_bonds" in d:
            out[0]["covalent_bonds"] = d["covalent_bonds"]

        dest.write_text(json.dumps(out, indent=2) + "\n")
        cloned += 1
        print(f"  cloned {name}.json <- {donor_path.stem} "
              f"(ligand: {d['sequences'][1].get('ligand', {}).get('ligand', '?')[:40]}...)")

    print(f"\ncloned={cloned} skipped_exists={skipped} "
          f"unmatched_seq={unmatched}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
