#!/usr/bin/env python3
"""Build Protenix input JSONs from CASP15 *_lig.pdb ground-truth files.

CASP15 input shape is different from CASP16:
  - One combined PDB per target: targets_ligand/Targets_ligand/<TARGET>_lig.pdb
    contains all protein chains, (optional) RNA/DNA chains, and ligand HETATMs.
  - Ligand identity is a 3-letter PDB CCD code — we reference it as CCD_<id>
    in the Protenix JSON (Protenix will look it up in its own CCD cache at
    inference time). No SMILES conversion needed.

Per target we emit jsons/<TARGET>.json:

    [{
      "name": "<TARGET>",
      "sequences": [
        {"proteinChain": {"sequence": "<SEQ1>", "count": 2}},    # homodimer dedup
        {"ligand":       {"ligand": "CCD_SAH", "count": 2}},
        ...
      ]
    }]

Chains with the same sequence are merged and `count` is set accordingly —
this is the Protenix-native way to model homo-oligomers. Multiple copies
of the same CCD ligand likewise merge.

Usage:
    python3 build_protenix_json.py T1124             # one target -> jsons/T1124.json
    python3 build_protenix_json.py --all             # every *_lig.pdb in the CASP15 source
    python3 build_protenix_json.py --manifest manifest_A40.txt   # only targets in a manifest
    python3 build_protenix_json.py --skip-na         # drop RNA/DNA chains (e.g. for prot-lig-only batches)
"""
import argparse
import json
import sys
from collections import OrderedDict
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent
CASP15 = ROOT.parent / "casp15_ligands"
PDBDIR = CASP15 / "targets_ligand" / "Targets_ligand"
JSON_DIR = ROOT / "jsons"

# HETATM residue names that should NOT be emitted as a ligand.
# (Mirrors the chai build_chai_fasta.py exclusion list.)
SKIP_AS_LIGAND = {"HOH", "WAT", "DOD", "TYR"}
RNA_RES = {"A", "U", "G", "C"}
DNA_RES = {"DA", "DT", "DG", "DC"}


def extract_entities(pdb_path: Path, skip_na: bool = False):
    """Parse a CASP15 _lig.pdb into (chains, ligand_ccds).

    chains: list of (chain_id, type, sequence) — type is 'protein'/'rna'/'dna'
    ligand_ccds: list of CCD codes (one per *instance* of a ligand in the file)
    """
    s = gemmi.read_structure(str(pdb_path))
    model = s[0]
    chains: list[tuple[str, str, str]] = []
    ligand_ccds: list[str] = []
    seen_lig_per_chain: dict[str, set] = {}

    for chain in model:
        prot_seq, na_seq, na_kind = [], [], None
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                prot_seq.append(info.one_letter_code.upper())
            elif res.name in RNA_RES:
                na_seq.append(res.name)
                na_kind = na_kind or "rna"
            elif res.name in DNA_RES:
                na_seq.append(res.name[1])  # strip leading 'D' → A/T/G/C
                na_kind = na_kind or "dna"
            elif res.het_flag == "H" and res.name not in SKIP_AS_LIGAND:
                key = (chain.name, res.seqid.num, res.name)
                seen = seen_lig_per_chain.setdefault(chain.name, set())
                if key not in seen:
                    seen.add(key)
                    ligand_ccds.append(res.name)
        if prot_seq:
            chains.append((chain.name, "protein", "".join(prot_seq)))
        if na_seq and not skip_na:
            chains.append((chain.name, na_kind, "".join(na_seq)))
    return chains, ligand_ccds


def build_json(target: str, pdb_path: Path, skip_na: bool) -> list:
    chains, ligand_ccds = extract_entities(pdb_path, skip_na=skip_na)

    # Dedup chains by (type, sequence) → Protenix count
    seq_counts: "OrderedDict[tuple[str,str], int]" = OrderedDict()
    for _cid, kind, seq in chains:
        seq_counts[(kind, seq)] = seq_counts.get((kind, seq), 0) + 1

    sequences = []
    for (kind, seq), count in seq_counts.items():
        if kind == "protein":
            entity = {"proteinChain": {"sequence": seq, "count": count}}
        elif kind == "rna":
            entity = {"rnaSequence": {"sequence": seq, "count": count}}
        elif kind == "dna":
            entity = {"dnaSequence": {"sequence": seq, "count": count}}
        else:
            raise ValueError(f"unknown chain kind: {kind}")
        sequences.append(entity)

    # Dedup ligands by CCD code → count
    from collections import Counter
    for ccd, count in Counter(ligand_ccds).items():
        sequences.append({"ligand": {"ligand": f"CCD_{ccd}", "count": count}})

    return [{"name": target, "sequences": sequences}]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", help="single target id (e.g. T1124)")
    ap.add_argument("--all", action="store_true",
                    help="emit JSON for every *_lig.pdb in the CASP15 source dir")
    ap.add_argument("--manifest", help="only emit targets listed in this manifest")
    ap.add_argument("--skip-na", action="store_true",
                    help="omit RNA/DNA chains (useful for prot-lig only batches)")
    args = ap.parse_args()

    if args.all or args.manifest:
        targets = []
        if args.manifest:
            with open(args.manifest) as f:
                wanted = {line.strip() for line in f
                          if line.strip() and not line.startswith("#")}
        else:
            wanted = None
        for pdb in sorted(PDBDIR.glob("*_lig.pdb")):
            t = pdb.stem.replace("_lig", "")
            if wanted is None or t in wanted:
                targets.append((t, pdb))
    elif args.target:
        pdb = PDBDIR / f"{args.target}_lig.pdb"
        if not pdb.exists():
            ap.error(f"no such target file: {pdb}")
        targets = [(args.target, pdb)]
    else:
        ap.error("specify a target, --all, or --manifest")

    JSON_DIR.mkdir(exist_ok=True)
    ok = 0
    for name, pdb in targets:
        try:
            data = build_json(name, pdb, skip_na=args.skip_na)
        except Exception as e:
            print(f"[WARN] {name}: {e}", file=sys.stderr)
            continue
        out = JSON_DIR / f"{name}.json"
        out.write_text(json.dumps(data, indent=2) + "\n")
        print(f"wrote {out}")
        ok += 1
    print(f"\n{ok} JSONs in {JSON_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
