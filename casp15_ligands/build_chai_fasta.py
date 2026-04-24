#!/usr/bin/env python3
"""Build a chai-lab FASTA from a CASP15 ligand target.

Reads `targets_ligand/Targets_ligand/<TARGET>_lig.pdb` and writes a FASTA
with one entry per protein chain and one entry per unique ligand
(SMILES fetched from RCSB Chemical Component Dictionary).

Usage:
    python build_chai_fasta.py T1124
    python build_chai_fasta.py T1124 --out fastas/T1124.fasta
    python build_chai_fasta.py --all          # write fastas/<TARGET>.fasta for every target
"""
import argparse
import json
import os
import sys
import urllib.request
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent
PDBDIR = ROOT / "targets_ligand" / "Targets_ligand"
SMILES_CACHE = ROOT / ".smiles_cache.json"


def load_cache() -> dict:
    if SMILES_CACHE.exists():
        return json.loads(SMILES_CACHE.read_text())
    return {}


def save_cache(cache: dict) -> None:
    SMILES_CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))


def fetch_smiles(ccd_id: str, cache: dict) -> str:
    """Fetch isomeric SMILES for a 3-letter PDB CCD code."""
    if ccd_id in cache:
        return cache[ccd_id]
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_id}"
    with urllib.request.urlopen(url, timeout=10) as r:
        data = json.load(r)
    desc = data.get("rcsb_chem_comp_descriptor", {})
    smiles = desc.get("SMILES_stereo") or desc.get("SMILES")
    if not smiles:
        raise RuntimeError(f"No SMILES for {ccd_id}")
    cache[ccd_id] = smiles
    save_cache(cache)
    return smiles


# Standard amino acids and common modified residues that should be treated
# as part of the protein chain (chai-lab supports modified residues via
# the AAA(SEP)AAA notation, but we keep it simple here and skip them).
SKIP_AS_LIGAND = {"HOH", "WAT", "DOD", "TYR"}  # TYR sometimes appears as HETATM


def extract_chains_and_ligands(pdb_path: Path):
    """Return (chains: list[(chain_id, type, sequence)], ligands: list[ccd_id]).

    Chain type is one of: 'protein', 'rna', 'dna'.
    """
    s = gemmi.read_structure(str(pdb_path))
    model = s[0]
    chains: list[tuple[str, str, str]] = []
    ligands: list[str] = []
    seen_lig_per_chain: dict[str, set] = {}

    rna_residues = {"A", "U", "G", "C"}
    dna_residues = {"DA", "DT", "DG", "DC"}

    for chain in model:
        prot_seq, na_seq, na_kind = [], [], None
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                prot_seq.append(info.one_letter_code.upper())
            elif res.name in rna_residues:
                na_seq.append(res.name)
                na_kind = na_kind or "rna"
            elif res.name in dna_residues:
                na_seq.append(res.name[1])  # strip leading D for chai FASTA
                na_kind = na_kind or "dna"
            elif res.het_flag == "H" and res.name not in SKIP_AS_LIGAND:
                key = (chain.name, res.seqid.num, res.name)
                seen = seen_lig_per_chain.setdefault(chain.name, set())
                if key not in seen:
                    seen.add(key)
                    ligands.append(res.name)
        if prot_seq:
            chains.append((chain.name, "protein", "".join(prot_seq)))
        if na_seq:
            chains.append((chain.name, na_kind, "".join(na_seq)))
    return chains, ligands


def build_fasta(target: str, pdb_path: Path, cache: dict) -> str:
    chains, ligands = extract_chains_and_ligands(pdb_path)
    out: list[str] = []
    for chain_id, kind, seq in chains:
        out.append(f">{kind}|name={target}-chain{chain_id}")
        out.append(seq)
    for i, ccd in enumerate(ligands, start=1):
        smiles = fetch_smiles(ccd, cache)
        out.append(f">ligand|name={target}-lig{i}-{ccd}")
        out.append(smiles)
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", help="CASP15 target id (e.g. T1124)")
    ap.add_argument("--out", help="output FASTA path (default: print to stdout)")
    ap.add_argument("--all", action="store_true",
                    help="write fastas/<TARGET>.fasta for every *_lig.pdb")
    args = ap.parse_args()

    cache = load_cache()

    if args.all:
        outdir = ROOT / "fastas"
        outdir.mkdir(exist_ok=True)
        for pdb in sorted(PDBDIR.glob("*_lig.pdb")):
            target = pdb.stem.replace("_lig", "")
            try:
                fasta = build_fasta(target, pdb, cache)
            except Exception as e:
                print(f"[WARN] {target}: {e}", file=sys.stderr)
                continue
            (outdir / f"{target}.fasta").write_text(fasta)
            print(f"wrote {outdir / (target + '.fasta')}")
        return 0

    if not args.target:
        ap.error("target required (or use --all)")
    pdb_path = PDBDIR / f"{args.target}_lig.pdb"
    if not pdb_path.exists():
        ap.error(f"no such target file: {pdb_path}")
    fasta = build_fasta(args.target, pdb_path, cache)
    if args.out:
        Path(args.out).write_text(fasta)
        print(f"wrote {args.out}")
    else:
        sys.stdout.write(fasta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
