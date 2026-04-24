#!/usr/bin/env python3
"""Split each CASP15 ligand-target reference into receptor.pdb + ligand SDFs.

OpenStructure's ligand-scoring pipeline requires ligands as SDF when the
receptor is PDB, and gemmi's automatic PDB->CIF conversion doesn't always
classify carbohydrates/lipids as non-polymer entities. Splitting explicitly
sidesteps both problems and gives us per-ligand files keyed by residue.

For each <TARGET>_lig.pdb we write:
    refs/<TARGET>/receptor.pdb          # all ATOM records (protein)
    refs/<TARGET>/lig_<i>_<CCD>.sdf     # each unique HETATM ligand residue
    refs/<TARGET>/ligand_map.json       # metadata: chain, resid, ccd, sdf path

Bond orders come from the cached SMILES dictionary (.smiles_cache.json),
pulled in earlier by build_chai_fasta.py.

Usage:
    python3 prep_references.py                          # all targets
    python3 prep_references.py --targets T1124,T1146    # subset
"""
import argparse
import json
import sys
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent
REF_DIR = ROOT / "targets_ligand" / "Targets_ligand"
OUT_DIR = ROOT / "refs"
SMILES_CACHE = ROOT / ".smiles_cache.json"

SKIP = {"HOH", "WAT", "DOD"}


def load_smiles_cache() -> dict:
    if SMILES_CACHE.exists():
        return json.loads(SMILES_CACHE.read_text())
    return {}


def is_aa_or_na(residue_name: str) -> bool:
    info = gemmi.find_tabulated_residue(residue_name)
    return bool(info and (info.is_amino_acid() or info.is_nucleic_acid()))


def write_receptor(source_pdb: Path, out_pdb: Path) -> None:
    """Copy only ATOM/TER/END records from source_pdb, dropping HETATMs."""
    kept_prefixes = ("ATOM  ", "TER   ", "TER\n", "END\n", "END   ")
    out_lines = []
    with source_pdb.open() as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("TER"):
                out_lines.append(line)
    out_lines.append("END\n")
    out_pdb.write_text("".join(out_lines))


def residue_to_pdb_block(residue: gemmi.Residue, chain_name: str) -> str:
    """Build a minimal valid PDB fragment with only this residue's atoms."""
    s = gemmi.Structure()
    m = gemmi.Model("1")
    c = gemmi.Chain(chain_name[0])
    c.add_residue(residue.clone())
    m.add_chain(c)
    s.add_model(m)
    return s.make_pdb_string()


def ligand_to_sdf(pdb_block: str, template_smiles: str | None, out_sdf: Path) -> bool:
    """Write a single ligand residue to SDF.

    Skip bond-order template assignment — it's segfaulting on some ligands and
    OST's fault-tolerant mode reconstructs bonds from distances anyway. Raw
    atom coordinates are sufficient for OST's symmetry-aware matching.
    """
    raw = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False,
                               proximityBonding=True)
    if raw is None:
        return False
    out_sdf.write_text(Chem.MolToMolBlock(raw) + "$$$$\n")
    return True


def prep_target(target: str, smiles_cache: dict) -> dict:
    pdb_path = REF_DIR / f"{target}_lig.pdb"
    if not pdb_path.exists():
        return {"target": target, "status": "missing_reference"}
    out_dir = OUT_DIR / target
    out_dir.mkdir(parents=True, exist_ok=True)

    structure = gemmi.read_structure(str(pdb_path))
    write_receptor(pdb_path, out_dir / "receptor.pdb")

    ligand_entries = []
    i = 0
    for model in structure:
        for chain in model:
            for res in chain:
                if res.name in SKIP:
                    continue
                if is_aa_or_na(res.name):
                    continue
                i += 1
                ccd = res.name
                block = residue_to_pdb_block(res, chain.name)
                out_sdf = out_dir / f"lig_{i:02d}_{ccd}.sdf"
                smiles = smiles_cache.get(ccd)
                ok = ligand_to_sdf(block, smiles, out_sdf)
                ligand_entries.append({
                    "ccd": ccd,
                    "chain": chain.name,
                    "seqid": res.seqid.num,
                    "sdf": str(out_sdf.relative_to(ROOT)),
                    "with_template": smiles is not None,
                    "ok": ok,
                })

    (out_dir / "ligand_map.json").write_text(json.dumps(ligand_entries, indent=2))
    return {"target": target, "status": "ok", "n_ligands": len(ligand_entries)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", help="comma-separated target whitelist")
    args = ap.parse_args()

    cache = load_smiles_cache()
    OUT_DIR.mkdir(exist_ok=True)

    targets = sorted(p.stem.replace("_lig", "") for p in REF_DIR.glob("*_lig.pdb"))
    if args.targets:
        wanted = {t.strip() for t in args.targets.split(",")}
        targets = [t for t in targets if t in wanted]

    for t in targets:
        r = prep_target(t, cache)
        print(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
