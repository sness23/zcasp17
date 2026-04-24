#!/usr/bin/env python3
"""Adapter: canonical target.json → DynamicBind input layout (stub).

DynamicBind expects per-target:
    <target>/
        protein.fasta            # or protein.pdb if docking against a known structure
        ligand.sdf               # or a SMILES .smi file
        config.yaml              # DynamicBind config

This adapter produces that layout from the canonical target.json.
Extend as needed for specific DynamicBind versions.

Usage:
    python3 to_dynamicbind.py <target.json>... --out-dir DIR
"""
import argparse
import json
import shutil
import sys
from pathlib import Path


def convert(target_json_path: Path, out_root: Path) -> None:
    d = json.loads(target_json_path.read_text())
    name = d["target_id"]
    dst = out_root / name
    dst.mkdir(parents=True, exist_ok=True)

    # protein.fasta — concatenate all chain sequences in a standard FASTA
    # (DynamicBind doesn't distinguish chains the way Protenix does).
    lines = []
    for chain in d["chains"]:
        for i in range(chain["count"]):
            suffix = f"-{i+1}" if chain["count"] > 1 else ""
            lines.append(f">{name}_{chain['id']}{suffix}")
            lines.append(chain["sequence"])
    (dst / "protein.fasta").write_text("\n".join(lines) + "\n")

    # protein.pdb — copy the canonical receptor
    src_pdb = target_json_path.parent / d["receptor_pdb"]
    if src_pdb.exists():
        shutil.copy2(src_pdb, dst / "protein.pdb")

    # ligand — use first ligand's first SDF (DynamicBind does one ligand at a time)
    ligands = d.get("ligands", [])
    if ligands:
        L = ligands[0]
        if L.get("sdf_paths"):
            src_sdf = target_json_path.parent / L["sdf_paths"][0]
            if src_sdf.exists():
                shutil.copy2(src_sdf, dst / "ligand.sdf")
        # also dump SMILES for fallback
        if L.get("smiles"):
            (dst / "ligand.smi").write_text(f"{L['smiles']} {name}\n")

    # manifest of this target's inputs
    manifest = {
        "target_id": name,
        "n_chains": sum(c["count"] for c in d["chains"]),
        "ligand_smiles": ligands[0]["smiles"] if ligands else None,
        "source_target_json": str(target_json_path),
    }
    (dst / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target_jsons", nargs="+", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for src in args.target_jsons:
        try:
            convert(src, args.out_dir)
            tid = json.loads(src.read_text())["target_id"]
            print(f"wrote {args.out_dir / tid}/")
            ok += 1
        except Exception as e:
            print(f"[error] {src}: {e}", file=sys.stderr)
    print(f"\n{ok}/{len(args.target_jsons)} converted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
