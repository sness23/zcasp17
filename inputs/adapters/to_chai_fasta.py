#!/usr/bin/env python3
"""Adapter: canonical target.json → chai-lab input FASTA.

chai-lab expects a FASTA with header types:
    >protein|name=<name>
    <sequence>
    >rna|name=<name>
    <sequence>
    >dna|name=<name>
    <sequence>
    >ligand|name=<name>
    <SMILES>

Chain `count > 1` is expanded into N separate entries (chai doesn't
have a native "count" field). Ligands get one entry per copy, all
with the same SMILES.

Usage:
    python3 to_chai_fasta.py <target.json>... --out-dir DIR
"""
import argparse
import json
import sys
from pathlib import Path


def convert(target_json_path: Path) -> str:
    d = json.loads(target_json_path.read_text())
    name = d["target_id"]
    lines: list[str] = []

    # Single-receptor CASP16-style targets get a single ">protein|name=<T>-receptor"
    # header matching the legacy chai fasta convention. Multi-chain CASP15-style
    # targets emit per-chain headers using the ORIGINAL source chain IDs (not the
    # canonical-synthetic `chainA_379aa`).
    casp_year = d.get("casp_year")
    multi_chain = len(d["chains"]) > 1 or any(c["count"] > 1 for c in d["chains"])

    if not multi_chain and d["chains"]:
        # CASP16 pharma style: one receptor, one header
        c = d["chains"][0]
        tag = c["type"]
        lines.append(f">{tag}|name={name}-receptor")
        lines.append(c["sequence"])
    else:
        for chain in d["chains"]:
            tag = chain["type"]
            source_ids = chain.get("source_chains") or [chain.get("id", "chain")]
            for i, src_id in enumerate(source_ids):
                lines.append(f">{tag}|name={name}-chain{src_id}")
                lines.append(chain["sequence"])

    # Ligand headers: legacy template is `{target}-lig{i}-{LABEL}` where LABEL is
    # the CCD (CASP15) or the source name from the TSV (CASP16: "LIG", "201", etc.).
    global_i = 0
    for lig in d.get("ligands", []):
        smiles = lig.get("smiles")
        if not smiles:
            if lig.get("ccd"):
                smiles = f"ccd:{lig['ccd']}"
            else:
                print(f"  [warn] {name} ligand {lig.get('id','?')}: no smiles or ccd — skipped",
                      file=sys.stderr)
                continue
        label = lig.get("ccd") or lig.get("source_name") or lig["id"]
        for i in range(lig["count"]):
            global_i += 1
            # CASP15 legacy: one lig entry per crystal copy, numbered 1..N
            # CASP16 legacy: single lig entry per target (count=1 always)
            if casp_year == 16 and lig["count"] == 1:
                lines.append(f">ligand|name={name}-lig-{label}")
            else:
                lines.append(f">ligand|name={name}-lig{global_i}-{label}")
            lines.append(smiles)

    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target_jsons", nargs="+", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for src in args.target_jsons:
        try:
            fasta = convert(src)
            tid = json.loads(src.read_text())["target_id"]
            dst = args.out_dir / f"{tid}.fasta"
            dst.write_text(fasta)
            print(f"wrote {dst}")
            ok += 1
        except Exception as e:
            print(f"[error] {src}: {e}", file=sys.stderr)
    print(f"\n{ok}/{len(args.target_jsons)} converted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
