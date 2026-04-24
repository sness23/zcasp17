#!/usr/bin/env python3
"""Estimate Protenix `N_token` per target from canonical target.json files.

Protenix counts one token per:
  - protein residue
  - RNA / DNA nucleotide
  - ligand heavy atom

Multi-copy entities (chain `count > 1`, ligand `count > 1`) contribute
proportionally. This script walks the canonical `inputs/casp{15,16,17}/`
tree (or an explicit list of target.json paths), computes the estimate,
and recommends a GPU tier based on the thresholds in
`~/data/vaults/docs/GUIDE-gpu-selection-protenix.md §4`.

Use this RIGHT AFTER a new target drops to decide which GPU to rent
before committing to a pod.

Usage:
    python3 estimate_tokens.py casp15/*/target.json casp16/*/target.json
    python3 estimate_tokens.py casp17/T1234/target.json --json
    python3 estimate_tokens.py --manifest manifest.csv  # all targets in manifest.csv
    python3 estimate_tokens.py --summary                # aggregate histogram
"""
import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent    # inputs/

# RDKit is the accurate way to count heavy atoms from SMILES.
# Fallback regex counts atomic symbols crudely if RDKit isn't available.
try:
    from rdkit import Chem
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False


def heavy_atoms_from_smiles(smiles: str) -> int:
    """Return count of non-hydrogen atoms in a SMILES string."""
    if not smiles:
        return 0
    if RDKIT_OK:
        try:
            m = Chem.MolFromSmiles(smiles)
            if m is not None:
                return m.GetNumHeavyAtoms()
        except Exception:
            pass
    # fallback: count uppercase-start element tokens; approximate
    # (recognizes Cl, Br, Si as two-char; misses rare elements)
    import re
    two_char = {"Cl", "Br", "Si", "Na", "Mg", "Al", "Se", "Zn", "Fe",
                "Ca", "Mn", "Co", "Ni", "Cu", "As", "Hg", "Pb"}
    s = smiles
    # strip bracket contents but count the element inside
    tokens = re.findall(r"Cl|Br|Si|[A-Z][a-z]?", s.replace("[", "").replace("]", ""))
    return sum(1 for t in tokens if t.upper() not in ("H",))


def ccd_to_heavy_atoms(ccd: str) -> int:
    """Rough heavy-atom count for common CCD codes. Approximation when
    SMILES is missing."""
    common = {
        "SAH": 26, "ATP": 31, "ADP": 27, "GTP": 32, "GDP": 28,
        "NAG": 14, "BMA": 11, "MAN": 11, "GAL": 11, "GLC": 11,
        "EPE": 15, "MPD": 8, "COA": 48,
        "HEM": 43, "NAD": 44, "NDP": 44, "FAD": 53,
        "HOH": 1, "K": 1, "NA": 1, "MG": 1, "CA": 1, "ZN": 1, "FE": 1,
        "CL": 1, "BR": 1, "CO": 1,
        "DW0": 30, "CD": 1,
        # pharma CCDs from CASP16 and beyond typically have SMILES in
        # the target.json so we rarely fall back to this.
    }
    return common.get(ccd.upper(), 25)  # 25 = default drug-sized guess


def estimate(target_json_path: Path) -> dict:
    d = json.loads(target_json_path.read_text())
    if isinstance(d, list):
        d = d[0]

    # Protein / NA residue tokens
    poly_tokens = 0
    n_chains = 0
    for chain in d.get("chains", []):
        n_chains += chain.get("count", 1)
        poly_tokens += len(chain["sequence"]) * chain.get("count", 1)

    # Ligand heavy-atom tokens
    lig_tokens = 0
    n_lig_copies = 0
    for lig in d.get("ligands", []):
        count = lig.get("count", 1)
        n_lig_copies += count
        # Prefer SMILES (accurate). Fall back to CCD lookup.
        if lig.get("smiles"):
            atoms = heavy_atoms_from_smiles(lig["smiles"])
        elif lig.get("ccd"):
            atoms = ccd_to_heavy_atoms(lig["ccd"])
        else:
            atoms = 25
        lig_tokens += atoms * count

    total = poly_tokens + lig_tokens
    return {
        "target_id": d.get("target_id", target_json_path.parent.name),
        "casp_year": d.get("casp_year"),
        "set": d.get("set") or "",
        "n_chains": n_chains,
        "poly_tokens": poly_tokens,
        "lig_tokens": lig_tokens,
        "total_tokens": total,
        "n_lig_copies": n_lig_copies,
        "recommended_gpu": recommend_gpu(total),
    }


def recommend_gpu(n_token: int) -> str:
    """Map token count to GPU tier per GUIDE-gpu-selection-protenix.md §4."""
    if n_token <= 1024:
        return "24GB (RTX 3090/4090)"
    if n_token <= 2048:
        return "48GB (A40/L40S/6000Ada)"
    if n_token <= 4000:
        return "48GB (L40S preferred for wall time)"
    if n_token <= 6000:
        return "80GB (A100/H100)"
    return "80GB + mitigations (§5 of guide)"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target_jsons", nargs="*", type=Path,
                    help="paths to target.json files (or use --manifest/--all)")
    ap.add_argument("--all", action="store_true",
                    help="scan every target.json under inputs/casp*/")
    ap.add_argument("--manifest", type=Path,
                    help="CSV with column 'target_id' — restrict to these targets")
    ap.add_argument("--json", action="store_true",
                    help="emit JSON instead of a table")
    ap.add_argument("--summary", action="store_true",
                    help="append an aggregate summary (distribution + counts per bucket)")
    ap.add_argument("--sort", default="total",
                    choices=["total", "name", "set"],
                    help="sort order for the table")
    args = ap.parse_args()

    if not RDKIT_OK:
        print("[warn] RDKit not importable — falling back to regex heavy-atom count "
              "(less accurate). Install rdkit in the active env for precise results.",
              file=sys.stderr)

    paths: list[Path] = []
    if args.all:
        for year_dir in sorted(ROOT.glob("casp*")):
            paths.extend(sorted(year_dir.glob("*/target.json")))
    elif args.target_jsons:
        paths = args.target_jsons
    elif args.manifest:
        with args.manifest.open() as f:
            wanted = {r["target_id"] for r in csv.DictReader(f)}
        for year_dir in sorted(ROOT.glob("casp*")):
            for tj in year_dir.glob("*/target.json"):
                if tj.parent.name in wanted:
                    paths.append(tj)
    else:
        ap.error("provide target_jsons, or --all, or --manifest")

    rows = []
    for p in paths:
        try:
            rows.append(estimate(p))
        except Exception as e:
            print(f"[error] {p}: {e}", file=sys.stderr)

    # sort
    if args.sort == "total":
        rows.sort(key=lambda r: r["total_tokens"])
    elif args.sort == "name":
        rows.sort(key=lambda r: r["target_id"])
    elif args.sort == "set":
        rows.sort(key=lambda r: (r["casp_year"], r["set"], r["target_id"]))

    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        print(f"{'target':12} {'year':>4} {'set':>7} {'chains':>6} "
              f"{'poly':>5} {'lig':>5} {'total':>6}  GPU")
        print("-" * 78)
        for r in rows:
            print(f"{r['target_id']:12} {r['casp_year']:>4} {r['set']:>7} "
                  f"{r['n_chains']:>6} {r['poly_tokens']:>5} "
                  f"{r['lig_tokens']:>5} {r['total_tokens']:>6}  "
                  f"{r['recommended_gpu']}")

    if args.summary:
        print("\n=== aggregate ===")
        print(f"total targets: {len(rows)}")
        print(f"total-token stats: min={min(r['total_tokens'] for r in rows)}  "
              f"median={sorted(r['total_tokens'] for r in rows)[len(rows)//2]}  "
              f"max={max(r['total_tokens'] for r in rows)}")
        buckets = {
            "≤1024 (24GB)":         0,
            "1025-2048 (48GB)":     0,
            "2049-4000 (48GB hot)": 0,
            "4001-6000 (80GB)":     0,
            ">6000 (mitigations)":  0,
        }
        for r in rows:
            t = r["total_tokens"]
            if   t <= 1024: buckets["≤1024 (24GB)"] += 1
            elif t <= 2048: buckets["1025-2048 (48GB)"] += 1
            elif t <= 4000: buckets["2049-4000 (48GB hot)"] += 1
            elif t <= 6000: buckets["4001-6000 (80GB)"] += 1
            else:           buckets[">6000 (mitigations)"] += 1
        for k, v in buckets.items():
            bar = "#" * min(60, v)
            print(f"  {k:26} {v:>5}  {bar}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
