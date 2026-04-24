#!/usr/bin/env python3
"""Prep CASP16 pharma reference structures into receptor.pdb + per-ligand SDFs
for `ost compare-ligand-structures`.

CASP16 layout (per target):
    casp16/pharma_ligands/exper_struct/<SET>/<SET>_prepared/<TARGET>/
        protein_aligned.pdb
        ligand_<resid>_<chain>_<copy>.pdb        # one crystal ligand pose

We write:
    casp16_ligands/refs/<TARGET>/
        receptor.pdb                              # filtered to ATOM records
        lig_01_<NAME>.sdf                         # crystal ligand pose
        ligand_map.json                           # metadata

Usage:
    python3 prep_references.py --targets L2001,L2002
    python3 prep_references.py --set L2000
    python3 prep_references.py --all
"""
import argparse
import json
import sys
from pathlib import Path

import gemmi
from rdkit import Chem

ROOT = Path(__file__).resolve().parent
CASP16 = ROOT.parent / "casp16"
EXPER = CASP16 / "pharma_ligands" / "exper_struct"
OUT_DIR = ROOT / "refs"


def set_for(target: str) -> str:
    return {"L1": "L1000", "L2": "L2000", "L3": "L3000", "L4": "L4000"}[target[:2]]


def write_receptor(src_pdb: Path, out_pdb: Path) -> None:
    lines = []
    with src_pdb.open() as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("TER"):
                lines.append(line)
    lines.append("END\n")
    out_pdb.write_text("".join(lines))


def residue_to_pdb_block(residue: gemmi.Residue, chain_name: str) -> str:
    s = gemmi.Structure()
    m = gemmi.Model("1")
    c = gemmi.Chain(chain_name[0])
    c.add_residue(residue.clone())
    m.add_chain(c)
    s.add_model(m)
    return s.make_pdb_string()


def ligand_pdb_to_sdf(src_pdb: Path, out_sdf: Path) -> bool:
    """Convert a ligand PDB file to SDF. Bonds inferred from proximity; OST
    fault-tolerant mode reconstructs bond orders anyway."""
    struct = gemmi.read_structure(str(src_pdb))
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    # The prepared ligand PDBs contain a single HETATM residue.
    for model in struct:
        for chain in model:
            for res in chain:
                block = residue_to_pdb_block(res, chain.name)
                raw = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False,
                                           proximityBonding=True)
                if raw is None:
                    return False
                out_sdf.write_text(Chem.MolToMolBlock(raw) + "$$$$\n")
                return True
    return False


def prep_target(target: str) -> dict:
    set_id = set_for(target)
    src_dir = EXPER / set_id / f"{set_id}_prepared" / target
    if not src_dir.is_dir():
        return {"target": target, "status": "missing_source"}
    receptor = src_dir / "protein_aligned.pdb"
    ligs = sorted(src_dir.glob("ligand_*.pdb"))
    if not receptor.exists() or not ligs:
        return {"target": target, "status": "missing_files"}

    out = OUT_DIR / target
    out.mkdir(parents=True, exist_ok=True)
    write_receptor(receptor, out / "receptor.pdb")

    map_entries = []
    for i, lig_pdb in enumerate(ligs, start=1):
        name = lig_pdb.stem  # e.g. ligand_201_C_1
        sdf = out / f"lig_{i:02d}_{name}.sdf"
        ok = ligand_pdb_to_sdf(lig_pdb, sdf)
        map_entries.append({"src": lig_pdb.name, "sdf": sdf.name, "ok": ok})
    (out / "ligand_map.json").write_text(json.dumps(map_entries, indent=2))
    return {"target": target, "status": "ok", "n_ligands": len(ligs)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", help="comma-separated list (e.g. L2001,L2002)")
    ap.add_argument("--set", dest="set_id", help="one of L1000/L2000/L3000/L4000")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    targets: list[str] = []
    if args.all:
        for sid in ("L1000", "L2000", "L3000", "L4000"):
            targets.extend(sorted(p.stem for p in (ROOT / sid).glob("*.tsv")))
    elif args.set_id:
        targets = sorted(p.stem for p in (ROOT / args.set_id).glob("*.tsv"))
    elif args.targets:
        targets = [t.strip() for t in args.targets.split(",")]
    else:
        ap.error("specify --all, --set, or --targets")

    OUT_DIR.mkdir(exist_ok=True)
    for t in targets:
        r = prep_target(t)
        print(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
