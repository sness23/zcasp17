#!/usr/bin/env python3
"""Score chai-lab predictions with OpenStructure's LDDT_pli + BiSyRMSD.

Wraps `ost compare-ligand-structures` — the same scoring tool CASP used.
For each target + model, runs the comparison and parses the JSON output.

Requires OpenStructure (we installed it in a conda env named `ost`) and
gemmi for converting the reference PDBs to CIF on the fly.

Usage:
    python3 score_lddt_pli.py results_A40 out/lddt_pli
    python3 score_lddt_pli.py results_A40 out/lddt_pli --best-only
    python3 score_lddt_pli.py results_A40 out/lddt_pli --targets T1124,T1146

The raw per-target JSON is preserved under <out_dir>/<target>/model_idx_<N>.json.
A combined CSV at <out_dir>/summary.csv is overwritten each run.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent
REF_DIR = ROOT / "targets_ligand" / "Targets_ligand"
PREP_REFS = ROOT / "refs"  # populated by prep_references.py
OST_BIN = Path(os.environ.get(
    "OST_BIN", os.path.expanduser("~/anaconda3/envs/ost/bin/ost")))


def pdb_to_cif(pdb_path: Path, out_cif: Path) -> None:
    """Convert a PDB reference to mmCIF (required by ost compare-ligand-structures)."""
    s = gemmi.read_structure(str(pdb_path))
    s.setup_entities()
    s.make_mmcif_document().write_file(str(out_cif))


def best_agg_model_idx(target_dir: Path) -> int | None:
    npzs = sorted(target_dir.glob("scores.model_idx_*.npz"))
    if not npzs:
        return None
    best, best_idx = -1.0, None
    for npz in npzs:
        d = np.load(npz)
        agg = float(d["aggregate_score"][0])
        if agg > best:
            best = agg
            best_idx = int(npz.stem.split("_")[-1])
    return best_idx


def run_ost(model_cif: Path, ref_pdb: Path, ref_ligand_sdfs: list[Path],
            out_json: Path) -> bool:
    """Run `ost compare-ligand-structures`. Return True on success."""
    cmd = [
        str(OST_BIN), "compare-ligand-structures",
        "-m", str(model_cif), "-mf", "cif",
        "-r", str(ref_pdb),   "-rf", "pdb",
        "-rl", *[str(p) for p in ref_ligand_sdfs],
        "-of", "json", "-ft",
        "--lddt-pli", "--rmsd",
        "-o", str(out_json),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        return False
    if not out_json.exists():
        return False
    try:
        d = json.loads(out_json.read_text())
        if d.get("status") == "FAILURE":
            return False
    except json.JSONDecodeError:
        return False
    return True


def parse_result(json_path: Path) -> list[dict]:
    """One row per model ligand. OST assigns LDDT_pli and RMSD independently
    (they can map the same model ligand to different reference ligands), so we
    index by model ligand and carry both assignments if present."""
    d = json.loads(json_path.read_text())
    rows: dict[str, dict] = {}

    for a in d.get("lddt_pli", {}).get("assigned_scores", []):
        m = a["model_ligand"]
        r = rows.setdefault(m, {"model_ligand": m})
        r["lddt_pli"] = a.get("score")
        r["lddt_pli_ref"] = Path(a["reference_ligand"]).name \
            if a.get("reference_ligand", "").endswith(".sdf") else a.get("reference_ligand")
        r["lddt_pli_n_contacts"] = a.get("lddt_pli_n_contacts")
        r["lddt_pli_coverage"] = a.get("coverage")

    for a in d.get("rmsd", {}).get("assigned_scores", []):
        m = a["model_ligand"]
        r = rows.setdefault(m, {"model_ligand": m})
        r["bisy_rmsd"] = a.get("score")
        r["bisy_rmsd_ref"] = Path(a["reference_ligand"]).name \
            if a.get("reference_ligand", "").endswith(".sdf") else a.get("reference_ligand")
        r["lddt_lp"] = a.get("lddt_lp")
        r["bs_bb_rmsd"] = a.get("bb_rmsd")
    return list(rows.values())


def score_target(target: str, target_results_dir: Path, out_dir: Path,
                 best_only: bool) -> list[dict]:
    ref_dir = PREP_REFS / target
    ref_pdb = ref_dir / "receptor.pdb"
    ref_ligs = sorted(ref_dir.glob("lig_*.sdf"))
    if not ref_pdb.exists():
        print(f"[{target}] skip: no prepped receptor at {ref_pdb} "
              f"(run prep_references.py first)", file=sys.stderr)
        return []
    if not ref_ligs:
        print(f"[{target}] skip: no prepped ligand SDFs in {ref_dir}",
              file=sys.stderr)
        return []

    if best_only:
        idx = best_agg_model_idx(target_results_dir)
        model_indices = [idx] if idx is not None else []
    else:
        model_indices = sorted(
            int(p.stem.split("_")[-1])
            for p in target_results_dir.glob("pred.model_idx_*.cif")
        )

    rows: list[dict] = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for idx in model_indices:
        model_cif = target_results_dir / f"pred.model_idx_{idx}.cif"
        out_json = out_dir / f"model_idx_{idx}.json"
        ok = run_ost(model_cif, ref_pdb, ref_ligs, out_json)
        if not ok:
            print(f"[{target}] model {idx}: ost failed", file=sys.stderr)
            continue
        parsed = parse_result(out_json)
        for r in parsed:
            r["target"] = target
            r["model_idx"] = idx
            rows.append(r)
        print(f"[{target}] model {idx}: scored {len(parsed)} ligand pair(s)")
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    if not rows:
        path.write_text("")
        return
    keys = ["target", "model_idx", "model_ligand",
            "lddt_pli", "lddt_pli_ref", "lddt_pli_n_contacts", "lddt_pli_coverage",
            "bisy_rmsd", "bisy_rmsd_ref", "lddt_lp", "bs_bb_rmsd"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in keys})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--best-only", action="store_true",
                    help="score only the best-aggregate model per target")
    ap.add_argument("--targets", help="comma-separated target whitelist")
    args = ap.parse_args()

    targets = sorted(
        p.name for p in args.results_dir.iterdir()
        if p.is_dir() and not p.name.startswith("_")
    )
    if args.targets:
        wanted = {t.strip() for t in args.targets.split(",")}
        targets = [t for t in targets if t in wanted]

    all_rows: list[dict] = []
    for t in targets:
        rows = score_target(
            t, args.results_dir / t, args.out_dir / t, args.best_only
        )
        all_rows.extend(rows)

    write_csv(all_rows, args.out_dir / "summary.csv")
    print(f"\nwrote {args.out_dir / 'summary.csv'} ({len(all_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
