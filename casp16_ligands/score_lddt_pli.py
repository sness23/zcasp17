#!/usr/bin/env python3
"""Score chai-lab CASP16 predictions with OpenStructure's LDDT_pli + BiSyRMSD.

Mirrors casp15_ligands/score_lddt_pli.py but reads prepped refs from
casp16_ligands/refs/<TARGET>/ (populated by prep_references.py).

Usage:
    python3 score_lddt_pli.py results_A40 out/lddt_pli
    python3 score_lddt_pli.py results_A40 out/lddt_pli --best-only
    python3 score_lddt_pli.py results_A40 out/lddt_pli --targets L2001,L2002
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PREP_REFS = ROOT / "refs"
OST_BIN = Path(os.environ.get(
    "OST_BIN", os.path.expanduser("~/anaconda3/envs/ost/bin/ost")))


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
    if res.returncode != 0 or not out_json.exists():
        return False
    try:
        d = json.loads(out_json.read_text())
        return d.get("status") != "FAILURE"
    except json.JSONDecodeError:
        return False


def parse_result(json_path: Path) -> list[dict]:
    d = json.loads(json_path.read_text())
    rows: dict[str, dict] = {}
    for a in d.get("lddt_pli", {}).get("assigned_scores", []):
        m = a["model_ligand"]
        r = rows.setdefault(m, {"model_ligand": m})
        r["lddt_pli"] = a.get("score")
        r["lddt_pli_ref"] = Path(a["reference_ligand"]).name \
            if a.get("reference_ligand", "").endswith(".sdf") else a.get("reference_ligand")
        r["lddt_pli_n_contacts"] = a.get("lddt_pli_n_contacts")
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
    if not ref_pdb.exists() or not ref_ligs:
        print(f"[{target}] skip: no prepped refs in {ref_dir}", file=sys.stderr)
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
        if not run_ost(model_cif, ref_pdb, ref_ligs, out_json):
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
            "lddt_pli", "lddt_pli_ref", "lddt_pli_n_contacts",
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
    ap.add_argument("--best-only", action="store_true")
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
        rows = score_target(t, args.results_dir / t, args.out_dir / t, args.best_only)
        all_rows.extend(rows)

    write_csv(all_rows, args.out_dir / "summary.csv")
    print(f"\nwrote {args.out_dir / 'summary.csv'} ({len(all_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
