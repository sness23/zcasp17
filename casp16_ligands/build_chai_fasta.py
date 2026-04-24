#!/usr/bin/env python3
"""Build chai-lab FASTAs for CASP16 pharma-track targets (L1001, L2002, …).

CASP16 pharma format is different from CASP15:
  - Receptor: one common PDB per set
    casp16/pharma_ligands/exper_struct/<SET>/<SET>_prepared/<TARGET>/protein_aligned.pdb
  - Ligand SMILES: one TSV per target (single data row)
    casp16_ligands/<SET>/<TARGET>.tsv
  - Crystal ligand: ligand_*.pdb in the same target folder

For each target we emit one chai-lab FASTA:
    >protein|name=<TARGET>-receptor
    <sequence from protein_aligned.pdb>
    >ligand|name=<TARGET>-lig-<Name>
    <SMILES>

Usage:
    python3 build_chai_fasta.py L2000           # write fastas/L2001.fasta, fastas/L2002.fasta
    python3 build_chai_fasta.py --all           # all four sets
    python3 build_chai_fasta.py L1000 --target L1001
"""
import argparse
import csv
import sys
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent
CASP16 = ROOT.parent / "casp16"
EXPER = CASP16 / "pharma_ligands" / "exper_struct"
FASTA_DIR = ROOT / "fastas"
SETS = ("L1000", "L2000", "L3000", "L4000")


def receptor_path(set_id: str, target: str) -> Path:
    return EXPER / set_id / f"{set_id}_prepared" / target / "protein_aligned.pdb"


def ligand_xtal_path(set_id: str, target: str) -> Path | None:
    """Return the first ligand_*.pdb in the target folder (CASP uses ligand_<resid>_<chain>_<copy>.pdb)."""
    d = EXPER / set_id / f"{set_id}_prepared" / target
    hits = sorted(d.glob("ligand_*.pdb"))
    return hits[0] if hits else None


def protein_sequence_from_pdb(pdb_path: Path) -> str:
    """One-letter protein sequence from a prepared PDB. Concatenate chains in file order."""
    s = gemmi.read_structure(str(pdb_path))
    out = []
    for chain in s[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                out.append(info.one_letter_code.upper())
    return "".join(out)


def parse_tsv(tsv_path: Path) -> tuple[str, str] | None:
    """Return (ligand_name, SMILES) from a CASP16 pharma target TSV."""
    with tsv_path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        return None
    row = rows[0]
    return row["Name"], row["SMILES"]


def build_fasta(set_id: str, target: str) -> str | None:
    rec = receptor_path(set_id, target)
    tsv = ROOT / set_id / f"{target}.tsv"
    if not rec.exists() or not tsv.exists():
        print(f"[warn] missing inputs for {target} (rec={rec.exists()} tsv={tsv.exists()})",
              file=sys.stderr)
        return None
    seq = protein_sequence_from_pdb(rec)
    parsed = parse_tsv(tsv)
    if parsed is None:
        print(f"[warn] empty TSV for {target}", file=sys.stderr)
        return None
    name, smiles = parsed
    lines = [
        f">protein|name={target}-receptor",
        seq,
        f">ligand|name={target}-lig-{name}",
        smiles,
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("set_id", nargs="?", help="L1000 / L2000 / L3000 / L4000")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--target", help="single target within the set (e.g. L1001)")
    args = ap.parse_args()

    if args.all:
        selected_sets = SETS
    elif args.set_id:
        selected_sets = (args.set_id,)
    else:
        ap.error("specify a set or --all")

    FASTA_DIR.mkdir(exist_ok=True)
    total = 0
    for sid in selected_sets:
        tsvs = sorted((ROOT / sid).glob("*.tsv"))
        if args.target:
            tsvs = [p for p in tsvs if p.stem == args.target]
        for tsv in tsvs:
            target = tsv.stem
            fasta = build_fasta(sid, target)
            if fasta is None:
                continue
            out = FASTA_DIR / f"{target}.fasta"
            out.write_text(fasta)
            total += 1
            print(f"wrote {out}")
    print(f"\n{total} FASTAs in {FASTA_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
