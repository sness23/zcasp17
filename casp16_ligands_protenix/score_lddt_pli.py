#!/usr/bin/env python3
"""Score Protenix CASP16 predictions with OpenStructure's LDDT_pli + BiSyRMSD.

Adapter for casp16_ligands/score_lddt_pli.py — same OST command, same refs
(symlinked at ./refs -> ../casp16_ligands/refs), but walks Protenix's nested
output layout:

    <results_dir>/<TARGET>/<SEED>/<TARGET>_<SEED>_sample_<N>.cif
    <results_dir>/<TARGET>/<SEED>/<TARGET>_<SEED>_summary_confidence_sample_<N>.json

"Best" sample selection uses `ranking_score` from the summary-confidence JSON
(Protenix's equivalent of chai's `aggregate_score`).

Usage:
    python3 score_lddt_pli.py results_L40S out/lddt_pli
    python3 score_lddt_pli.py results_L40S out/lddt_pli --best-only
    python3 score_lddt_pli.py results_L40S out/lddt_pli --targets L2001,L2002
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PREP_REFS = ROOT / "refs"
OST_BIN = Path(os.environ.get(
    "OST_BIN", os.path.expanduser("~/anaconda3/envs/ost/bin/ost")))

# <TARGET>_sample_<N>.cif   (seed comes from the seed_<N>/ parent dir)
CIF_RE = re.compile(r"^(?P<target>.+?)_sample_(?P<n>\d+)\.cif$")


def iter_samples(target_dir: Path):
    """Yield (seed:int, sample:int, cif_path, summary_json_path) for a target.

    Actual Protenix output layout:
        <target_dir>/<name>/seed_<N>/predictions/<name>_sample_<k>.cif
    We rglob for CIFs and derive seed from the seed_<N>/ parent dir.
    """
    for cif in sorted(target_dir.rglob("*_sample_*.cif")):
        m = CIF_RE.match(cif.name)
        if not m:
            continue
        n = int(m.group("n"))
        seed_dir = cif.parent.parent
        if not seed_dir.name.startswith("seed_"):
            continue
        try:
            seed = int(seed_dir.name.split("_", 1)[1])
        except ValueError:
            continue
        summary = cif.parent / cif.name.replace(
            f"_sample_{n}.cif", f"_summary_confidence_sample_{n}.json"
        )
        yield seed, n, cif, summary


def best_sample(target_dir: Path):
    """Return (seed, n, cif, summary) with the highest ranking_score."""
    best, best_tuple = -float("inf"), None
    for seed, n, cif, summary in iter_samples(target_dir):
        if not summary.exists():
            continue
        try:
            score = float(json.loads(summary.read_text()).get("ranking_score", -1))
        except (json.JSONDecodeError, TypeError):
            continue
        if score > best:
            best, best_tuple = score, (seed, n, cif, summary)
    return best_tuple


def run_ost(model_cif: Path, ref_pdb: Path, ref_ligand_sdfs: list[Path],
            out_json: Path) -> bool:
    cmd = [
        str(OST_BIN), "compare-ligand-structures",
        "-m", str(model_cif), "-mf", "cif",
        "-r", str(ref_pdb),   "-rf", "pdb",
        "-rl", *[str(p) for p in ref_ligand_sdfs],
        "-of", "json", "-ft",
        "--lddt-pli", "--rmsd",
        "--substructure-match",      # tolerate partial-resolution / tautomer / bond-order
        "-o", str(out_json),          # mismatches when matching model ligand to SDF ref
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_json.exists():
        return False
    try:
        return json.loads(out_json.read_text()).get("status") != "FAILURE"
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
        picked = best_sample(target_results_dir)
        samples = [picked] if picked is not None else []
    else:
        samples = list(iter_samples(target_results_dir))

    rows: list[dict] = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for seed, n, cif, summary in samples:
        out_json = out_dir / f"seed{seed}_sample{n}.json"
        if not run_ost(cif, ref_pdb, ref_ligs, out_json):
            print(f"[{target}] seed {seed} sample {n}: ost failed", file=sys.stderr)
            continue
        parsed = parse_result(out_json)
        # Pull ranking_score for cross-sample ordering in the CSV
        rscore = None
        if summary.exists():
            try:
                rscore = float(json.loads(summary.read_text()).get("ranking_score", None))
            except (json.JSONDecodeError, TypeError):
                pass
        for r in parsed:
            r["target"] = target
            r["seed"] = seed
            r["sample"] = n
            r["ranking_score"] = rscore
            rows.append(r)
        print(f"[{target}] seed {seed} sample {n}: scored {len(parsed)} ligand pair(s)")
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    if not rows:
        path.write_text("")
        return
    keys = ["target", "seed", "sample", "ranking_score", "model_ligand",
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
    ap.add_argument("--best-only", action="store_true",
                    help="only score the top-ranking_score sample per target")
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
