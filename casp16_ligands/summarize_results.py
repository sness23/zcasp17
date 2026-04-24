#!/usr/bin/env python3
"""Summarize chai-lab batch results.

Walks a results directory of `<target>/scores.model_idx_*.npz` files and
prints a per-target ranked summary plus an overall table.

Usage:
    python summarize_results.py /path/to/results
    python summarize_results.py /path/to/results --md > summary.md
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np


def per_target(target_dir: Path) -> dict | None:
    score_files = sorted(target_dir.glob("scores.model_idx_*.npz"))
    if not score_files:
        return None
    rows = []
    for sf in score_files:
        d = np.load(sf)
        idx = int(sf.stem.split("_")[-1])
        rows.append({
            "model": idx,
            "agg":  float(d["aggregate_score"][0]),
            "ptm":  float(d["ptm"][0]),
            "iptm": float(d["iptm"][0]),
            "clash": bool(d["has_inter_chain_clashes"][0]),
        })
    rows.sort(key=lambda r: -r["agg"])  # best first
    return {"target": target_dir.name, "models": rows}


def fmt_text(targets: list[dict]) -> str:
    out = []
    out.append(f"{'Target':<12}{'Best agg':>10}{'pTM':>7}{'ipTM':>7}{'Clash':>7}  Best CIF")
    out.append("-" * 72)
    for t in targets:
        best = t["models"][0]
        cif = f"{t['target']}/pred.model_idx_{best['model']}.cif"
        clash = "YES" if best["clash"] else "no"
        out.append(f"{t['target']:<12}{best['agg']:>10.3f}{best['ptm']:>7.3f}{best['iptm']:>7.3f}{clash:>7}  {cif}")
    return "\n".join(out) + "\n"


def fmt_md(targets: list[dict]) -> str:
    out = ["| Target | Best agg | pTM | ipTM | Clash | Best CIF |",
           "|---|---:|---:|---:|---|---|"]
    for t in targets:
        best = t["models"][0]
        cif = f"{t['target']}/pred.model_idx_{best['model']}.cif"
        clash = "⚠️ YES" if best["clash"] else "no"
        out.append(f"| {t['target']} | {best['agg']:.3f} | {best['ptm']:.3f} | {best['iptm']:.3f} | {clash} | `{cif}` |")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--md", action="store_true", help="output markdown table")
    ap.add_argument("--json", action="store_true", help="output JSON")
    args = ap.parse_args()

    if not args.results_dir.is_dir():
        print(f"not a directory: {args.results_dir}", file=sys.stderr)
        return 2

    targets = []
    for sub in sorted(args.results_dir.iterdir()):
        if sub.name.startswith("_") or not sub.is_dir():
            continue
        t = per_target(sub)
        if t:
            targets.append(t)

    if not targets:
        print("no completed targets found", file=sys.stderr)
        return 1

    targets.sort(key=lambda t: -t["models"][0]["agg"])

    if args.json:
        print(json.dumps(targets, indent=2))
    elif args.md:
        print(fmt_md(targets))
    else:
        print(fmt_text(targets))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
