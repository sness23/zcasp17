#!/usr/bin/env python3
"""Populate the inputs/ canonical tree from the existing CASP15 and CASP16
source data under ../casp15_ligands/, ../casp16_ligands/, and ../casp16/.

Produces, per target:
    inputs/casp{15,16}/<TARGET>/
        target.json
        receptor.fasta
        receptor.pdb
        ligands.tsv
        ligands/lig_NN_<CCD>.sdf   (copied/symlinked from ../refs/)

Usage:
    python3 build_canonical.py --casp15
    python3 build_canonical.py --casp16
    python3 build_canonical.py --all
    python3 build_canonical.py --manifest          # regenerate manifest.csv
"""
import argparse
import csv
import json
import shutil
import sys
import urllib.request
from collections import Counter, OrderedDict
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parent                           # inputs/
CASP_ROOT = ROOT.parent                                          # ~/data/vaults/casp

# CASP15 source locations
CASP15_BUNDLE = CASP_ROOT / "casp15_ligands" / "targets_ligand" / "Targets_ligand"
CASP15_REFS = CASP_ROOT / "casp15_ligands" / "refs"

# CASP16 source locations
CASP16_EXPER = CASP_ROOT / "casp16" / "pharma_ligands" / "exper_struct"
CASP16_LIGANDS_ROOT = CASP_ROOT / "casp16_ligands"
CASP16_REFS = CASP_ROOT / "casp16_ligands" / "refs"
CASP16_SETS = ("L1000", "L2000", "L3000", "L4000")

# Residue sets
RNA_RES = {"A", "U", "G", "C"}
DNA_RES = {"DA", "DT", "DG", "DC"}
SKIP_AS_LIGAND = {"HOH", "WAT", "DOD", "TYR"}

SMILES_CACHE_PATH = ROOT / ".smiles_cache.json"


def load_smiles_cache() -> dict:
    if SMILES_CACHE_PATH.exists():
        return json.loads(SMILES_CACHE_PATH.read_text())
    return {}


def save_smiles_cache(cache: dict) -> None:
    SMILES_CACHE_PATH.write_text(json.dumps(cache, indent=2, sort_keys=True))


def fetch_ccd_smiles(ccd: str, cache: dict) -> str | None:
    """Fetch isomeric SMILES for a 3-letter PDB CCD code via RCSB. Cached."""
    if ccd in cache:
        return cache[ccd]
    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}"
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.load(r)
        desc = data.get("rcsb_chem_comp_descriptor", {})
        smiles = desc.get("SMILES_stereo") or desc.get("SMILES")
    except Exception as e:
        print(f"  [warn] CCD lookup failed for {ccd}: {e}", file=sys.stderr)
        smiles = None
    cache[ccd] = smiles
    save_smiles_cache(cache)
    return smiles


def extract_entities_from_pdb(pdb_path: Path, skip_na: bool = False):
    """Return (chains, ligand_ccds) from a CASP-style PDB.

    chains: list of (source_chain_id, type, sequence)
    ligand_ccds: list of CCD codes (one per HETATM instance in the file)
    """
    s = gemmi.read_structure(str(pdb_path))
    chains = []
    lig_ccds = []
    seen_lig_per_chain: dict[str, set] = {}

    for chain in s[0]:
        prot_seq, na_seq, na_kind = [], [], None
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                prot_seq.append(info.one_letter_code.upper())
            elif res.name in RNA_RES:
                na_seq.append(res.name)
                na_kind = na_kind or "rna"
            elif res.name in DNA_RES:
                na_seq.append(res.name[1])
                na_kind = na_kind or "dna"
            elif res.het_flag == "H" and res.name not in SKIP_AS_LIGAND:
                key = (chain.name, res.seqid.num, res.name)
                seen = seen_lig_per_chain.setdefault(chain.name, set())
                if key not in seen:
                    seen.add(key)
                    lig_ccds.append(res.name)
        if prot_seq:
            chains.append((chain.name, "protein", "".join(prot_seq)))
        if na_seq and not skip_na:
            chains.append((chain.name, na_kind, "".join(na_seq)))
    return chains, lig_ccds


def dedup_chains(raw_chains):
    """Merge identical sequences into count>1 entities, keeping source_chains."""
    merged: "OrderedDict[tuple[str,str], dict]" = OrderedDict()
    for src_id, kind, seq in raw_chains:
        key = (kind, seq)
        if key in merged:
            merged[key]["count"] += 1
            merged[key]["source_chains"].append(src_id)
        else:
            synthetic_id = f"chain{src_id}_{len(seq)}aa" if len(raw_chains) > 1 else src_id
            merged[key] = {
                "id": synthetic_id,
                "type": kind,
                "sequence": seq,
                "count": 1,
                "source_chains": [src_id],
            }
    return list(merged.values())


def write_receptor_pdb(source_pdb: Path, dst_pdb: Path) -> None:
    """Write source PDB with HETATMs stripped."""
    s = gemmi.read_structure(str(source_pdb))
    for model in s:
        for chain in model:
            to_remove = [i for i, r in enumerate(chain) if r.het_flag == "H"]
            for i in reversed(to_remove):
                del chain[i]
    s.write_pdb(str(dst_pdb))


def write_receptor_fasta(target_id: str, chains: list, dst_fasta: Path) -> None:
    with dst_fasta.open("w") as f:
        for c in chains:
            hdr = f">{target_id}|{c['id']}|type={c['type']}|count={c['count']}"
            f.write(hdr + "\n")
            f.write(c["sequence"] + "\n")


def write_ligands_tsv(ligands: list, dst_tsv: Path) -> None:
    with dst_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["lig_id", "ccd", "smiles", "count", "sdf_paths"])
        for L in ligands:
            w.writerow([
                L["id"],
                L.get("ccd") or "",
                L.get("smiles") or "",
                L["count"],
                ";".join(L["sdf_paths"]),
            ])


def build_casp15_target(target: str, lig_pdb: Path, smiles_cache: dict) -> tuple[dict, list[Path]] | None:
    """Build target.json for a CASP15 target. Returns (target_dict, sdf_files_to_copy)."""
    raw_chains, lig_ccds = extract_entities_from_pdb(lig_pdb)
    if not raw_chains:
        print(f"  [skip] {target}: no protein/NA chains")
        return None
    chains = dedup_chains(raw_chains)

    # Ligands: group by CCD, resolve SMILES via RCSB
    ligand_entities = []
    ccd_counts = Counter(lig_ccds)

    # Look up available SDFs in refs/<TARGET>/
    refs_dir = CASP15_REFS / target
    sdf_by_ccd: dict[str, list[Path]] = {}
    if refs_dir.exists():
        for sdf in sorted(refs_dir.glob("lig_*.sdf")):
            # filenames are like lig_01_SAH.sdf
            parts = sdf.stem.split("_")
            if len(parts) >= 3:
                ccd = parts[2]
                sdf_by_ccd.setdefault(ccd, []).append(sdf)

    # Iterate in first-seen (source) order, not count-descending, so the
    # canonical layout matches the order legacy builders used.
    for i, (ccd, n) in enumerate(ccd_counts.items(), start=1):
        smiles = fetch_ccd_smiles(ccd, smiles_cache)
        sdfs = sdf_by_ccd.get(ccd, [])
        ligand_entities.append({
            "id": f"lig_{i:02d}",
            "ccd": ccd,
            "smiles": smiles,
            "count": n,
            "sdf_paths": [f"ligands/{p.name}" for p in sdfs[:n]],
        })

    sdf_files = []
    for entity in ligand_entities:
        sdfs = sdf_by_ccd.get(entity["ccd"], [])
        sdf_files.extend(sdfs[:entity["count"]])

    target_dict = {
        "target_id": target,
        "casp_year": 15,
        "set": None,
        "chains": chains,
        "ligands": ligand_entities,
        "receptor_pdb": "receptor.pdb",
        "notes": "",
    }
    return target_dict, sdf_files


def build_casp16_target(set_id: str, target: str) -> tuple[dict, list[Path]] | None:
    """Build target.json for a CASP16 target."""
    receptor_pdb_src = CASP16_EXPER / set_id / f"{set_id}_prepared" / target / "protein_aligned.pdb"
    tsv_src = CASP16_LIGANDS_ROOT / set_id / f"{target}.tsv"

    if not receptor_pdb_src.exists() or not tsv_src.exists():
        print(f"  [skip] {target}: missing source files")
        return None

    raw_chains, _ = extract_entities_from_pdb(receptor_pdb_src)
    chains = dedup_chains(raw_chains)

    # Read SMILES + ligand name from the CASP16 TSV
    rows = list(csv.DictReader(tsv_src.open(), delimiter="\t"))
    if not rows:
        return None
    lig_name, smiles = rows[0].get("Name", ""), rows[0].get("SMILES", "")

    # Find crystal ligand SDFs in refs/<TARGET>/
    refs_dir = CASP16_REFS / target
    sdfs = sorted(refs_dir.glob("lig_*.sdf")) if refs_dir.exists() else []

    ligand_entity = {
        "id": "lig_01",
        "ccd": None,                       # pharma targets lack CCDs
        "smiles": smiles,
        "count": len(sdfs) if sdfs else 1,
        "sdf_paths": [f"ligands/{p.name}" for p in sdfs],
    }
    if lig_name:
        ligand_entity["source_name"] = lig_name

    target_dict = {
        "target_id": target,
        "casp_year": 16,
        "set": set_id,
        "chains": chains,
        "ligands": [ligand_entity],
        "receptor_pdb": "receptor.pdb",
        "notes": f"{set_id} — {lig_name}" if lig_name else set_id,
    }
    return target_dict, list(sdfs)


def write_target_dir(root: Path, casp_year: int, target: str,
                      target_dict: dict, sdf_files: list[Path],
                      source_receptor_pdb: Path) -> None:
    dst = root / f"casp{casp_year}" / target
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "ligands").mkdir(exist_ok=True)

    # target.json
    (dst / "target.json").write_text(json.dumps(target_dict, indent=2) + "\n")

    # receptor.fasta + receptor.pdb
    write_receptor_fasta(target, target_dict["chains"], dst / "receptor.fasta")
    write_receptor_pdb(source_receptor_pdb, dst / "receptor.pdb")

    # ligands/*.sdf (copy so the canonical tree is self-contained)
    for sdf in sdf_files:
        shutil.copy2(sdf, dst / "ligands" / sdf.name)

    # ligands.tsv redundant view
    write_ligands_tsv(target_dict["ligands"], dst / "ligands.tsv")


def build_casp15(dst_root: Path, smiles_cache: dict) -> int:
    pdbs = sorted(CASP15_BUNDLE.glob("*_lig.pdb"))
    ok = 0
    for pdb in pdbs:
        target = pdb.stem.replace("_lig", "")
        print(f"[casp15] {target}")
        try:
            result = build_casp15_target(target, pdb, smiles_cache)
            if result is None:
                continue
            target_dict, sdfs = result
            write_target_dir(dst_root, 15, target, target_dict, sdfs, pdb)
            ok += 1
        except Exception as e:
            print(f"  [error] {target}: {e}", file=sys.stderr)
    return ok


def build_casp16(dst_root: Path) -> int:
    ok = 0
    for set_id in CASP16_SETS:
        tsvs = sorted((CASP16_LIGANDS_ROOT / set_id).glob("*.tsv"))
        for tsv in tsvs:
            target = tsv.stem
            print(f"[casp16/{set_id}] {target}")
            try:
                result = build_casp16_target(set_id, target)
                if result is None:
                    continue
                target_dict, sdfs = result
                source_pdb = CASP16_EXPER / set_id / f"{set_id}_prepared" / target / "protein_aligned.pdb"
                write_target_dir(dst_root, 16, target, target_dict, sdfs, source_pdb)
                ok += 1
            except Exception as e:
                print(f"  [error] {target}: {e}", file=sys.stderr)
    return ok


def build_manifest(dst_root: Path) -> int:
    """Walk inputs/casp*/ and emit manifest.csv."""
    rows = []
    for year_dir in sorted(dst_root.glob("casp*")):
        year = int(year_dir.name[4:])
        for tgt_dir in sorted(year_dir.iterdir()):
            if not (tgt_dir / "target.json").exists():
                continue
            d = json.loads((tgt_dir / "target.json").read_text())
            n_chains = sum(c["count"] for c in d["chains"])
            n_unique = len(d["chains"])
            n_lig_ent = len(d["ligands"])
            n_lig_copies = sum(L["count"] for L in d["ligands"])
            primary_ccd = (d["ligands"][0].get("ccd") if d["ligands"] else "") or ""
            rows.append({
                "casp_year": year,
                "set": d.get("set") or "",
                "target_id": d["target_id"],
                "n_chains": n_chains,
                "n_unique_seqs": n_unique,
                "n_ligand_entities": n_lig_ent,
                "n_ligand_copies": n_lig_copies,
                "primary_ccd": primary_ccd,
            })
    if not rows:
        print("no targets found in canonical tree", file=sys.stderr)
        return 0
    manifest = dst_root / "manifest.csv"
    with manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {manifest} ({len(rows)} rows)")
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--casp15", action="store_true")
    ap.add_argument("--casp16", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--manifest", action="store_true")
    ap.add_argument("--targets", help="comma-separated whitelist (scopes to year-dir already chosen)")
    args = ap.parse_args()

    if not any([args.casp15, args.casp16, args.all, args.manifest]):
        ap.error("specify --casp15 / --casp16 / --all / --manifest")

    cache = load_smiles_cache()
    total = 0

    if args.all or args.casp15:
        print("=== CASP15 ===")
        total += build_casp15(ROOT, cache)
    if args.all or args.casp16:
        print("=== CASP16 ===")
        total += build_casp16(ROOT)

    # always regenerate manifest at the end (unless user asked only for manifest)
    if args.manifest or total > 0:
        print("\n=== manifest ===")
        build_manifest(ROOT)

    print(f"\ndone: {total} target directories written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
