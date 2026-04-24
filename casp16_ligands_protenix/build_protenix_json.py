#!/usr/bin/env python3
"""Build Protenix input JSONs for CASP16 pharma-track targets.

Drop-in parallel to build_chai_fasta.py. Same inputs
(protein_aligned.pdb + <TARGET>.tsv), different output shape — Protenix
takes JSON not FASTA.

Per target we emit:

    [
      {
        "name": "<TARGET>",
        "sequences": [
          {"proteinChain": {"sequence": "<SEQ>", "count": 1}},
          {"ligand":       {"ligand":   "<SMILES>", "count": 1}}
        ]
      }
    ]

MSA and template paths are NOT filled in here — `protenix prep` on the
pod rewrites the JSON to add pairedMsaPath / unpairedMsaPath / templatesPath.

Usage:
    python3 build_protenix_json.py --all            # all four sets
    python3 build_protenix_json.py --set L2000      # one set
    python3 build_protenix_json.py --set L1000 --target L1001
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent
CASP16_LIGANDS = ROOT.parent / "casp16_ligands"
EXPER = ROOT.parent / "casp16" / "pharma_ligands" / "exper_struct"
JSON_DIR = ROOT / "jsons"
SETS = ("L1000", "L2000", "L3000", "L4000")


def receptor_path(set_id: str, target: str) -> Path:
    return EXPER / set_id / f"{set_id}_prepared" / target / "protein_aligned.pdb"


def tsv_path(set_id: str, target: str) -> Path:
    return CASP16_LIGANDS / set_id / f"{target}.tsv"


def chain_sequences_from_pdb(pdb_path: Path) -> list[str]:
    """Extract per-chain protein sequences (empty chains dropped)."""
    s = gemmi.read_structure(str(pdb_path))
    chain_seqs = []
    for chain in s[0]:
        seq = []
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                seq.append(info.one_letter_code.upper())
        if seq:
            chain_seqs.append("".join(seq))
    return chain_seqs


def parse_tsv(tsv: Path) -> tuple[str, str] | None:
    with tsv.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return None
    return rows[0]["Name"], rows[0]["SMILES"]


def build_json(set_id: str, target: str) -> list | None:
    rec = receptor_path(set_id, target)
    tsv = tsv_path(set_id, target)
    if not rec.exists() or not tsv.exists():
        print(
            f"[warn] missing inputs for {target} (rec={rec.exists()} tsv={tsv.exists()})",
            file=sys.stderr,
        )
        return None
    parsed = parse_tsv(tsv)
    if parsed is None:
        print(f"[warn] empty TSV for {target}", file=sys.stderr)
        return None
    _name, smiles = parsed

    # Dedup chains: identical sequences → one entity with count=N.
    # Distinct sequences (e.g. truncation variants) → separate entities.
    from collections import OrderedDict
    seq_counts: "OrderedDict[str, int]" = OrderedDict()
    for seq in chain_sequences_from_pdb(rec):
        seq_counts[seq] = seq_counts.get(seq, 0) + 1

    sequences = [
        {"proteinChain": {"sequence": seq, "count": n}}
        for seq, n in seq_counts.items()
    ]
    sequences.append({"ligand": {"ligand": smiles, "count": 1}})

    return [{"name": target, "sequences": sequences}]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="all four sets")
    ap.add_argument("--set", dest="set_id", choices=SETS, help="single set")
    ap.add_argument("--target", help="single target (e.g. L1001)")
    args = ap.parse_args()

    if args.all:
        sets = SETS
    elif args.set_id:
        sets = (args.set_id,)
    else:
        ap.error("specify --all or --set <L1000|L2000|L3000|L4000>")

    JSON_DIR.mkdir(exist_ok=True)
    total = 0
    for sid in sets:
        tsvs = sorted((CASP16_LIGANDS / sid).glob("*.tsv"))
        if args.target:
            tsvs = [p for p in tsvs if p.stem == args.target]
        for tsv in tsvs:
            target = tsv.stem
            data = build_json(sid, target)
            if data is None:
                continue
            out = JSON_DIR / f"{target}.json"
            out.write_text(json.dumps(data, indent=2) + "\n")
            total += 1
            print(f"wrote {out}")
    print(f"\n{total} JSONs in {JSON_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
