#!/usr/bin/env python3
"""Populate inputs/casp17/<T>/ from ../casp17_ligands/ (the CASP17 live spider).

CASP17 is a blind prediction track — we have:
  - protein.fasta  (sequence, from predictioncenter.org)
  - ligands.smi    (TSV with ID/Name/SMILES/Task, 1+ rows)
  - template.pdb   (placeholder or real structural template)

We do NOT have:
  - crystal ligand SDF (organizers release post-season)
  - ground-truth receptor PDB

So target.json has `sdf_paths: []` and `receptor_pdb: ""` — adapters
and the predictor itself will work; the scorer needs a separate code
path for blind predictions.

Usage:
    python3 build_casp17_from_spider.py                # all targets in spider
    python3 build_casp17_from_spider.py --target T1214
"""
import argparse
import csv
import json
import shutil
import sys
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parent                    # inputs/
SPIDER = ROOT.parent / "casp17_ligands"
SPIDER_TARGETS = SPIDER / "targets"


def parse_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Return list of (header, sequence) pairs."""
    entries = []
    header = None
    buf = []
    with fasta_path.open() as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(buf)))
                header = line[1:]
                buf = []
            else:
                buf.append(line)
        if header is not None:
            entries.append((header, "".join(buf)))
    return entries


def parse_ligands_smi(smi_path: Path) -> list[dict]:
    """Parse the CASP17 ligands.smi TSV (ID\\tName\\tSMILES\\tTask)."""
    ligs = []
    with smi_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        # strip trailing whitespace-filled columns
        reader.fieldnames = [fn.strip() for fn in reader.fieldnames if fn.strip()]
        for row in reader:
            smi = (row.get("SMILES") or "").strip()
            if not smi:
                continue
            ligs.append({
                "source_id": (row.get("ID") or "").strip(),
                "source_name": (row.get("Name") or "").strip(),
                "smiles": smi,
                "task": (row.get("Task") or "").strip(),
            })
    return ligs


def looks_like_placeholder_pdb(pdb_path: Path) -> bool:
    """CASP17 templates are sometimes all-zero placeholders before release."""
    if not pdb_path.exists():
        return True
    with pdb_path.open() as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # coords in cols 30-54 per PDB spec
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if (x, y, z) != (0.0, 0.0, 0.0):
                        return False
                except ValueError:
                    continue
    return True


def build_target(target: str) -> dict | None:
    src = SPIDER_TARGETS / target
    fasta = src / "protein.fasta"
    smi = src / "ligands.smi"
    template = src / "template.pdb"

    if not fasta.exists():
        print(f"  [skip] {target}: no protein.fasta", file=sys.stderr)
        return None

    # Chains from FASTA
    entries = parse_fasta(fasta)
    # Dedup by sequence, preserve first-seen order
    seq_counts: "OrderedDict[str, dict]" = OrderedDict()
    for i, (hdr, seq) in enumerate(entries):
        if seq not in seq_counts:
            src_id = hdr.split()[0] if hdr else f"chain{i+1}"
            seq_counts[seq] = {
                "id": f"chain{src_id}_{len(seq)}aa" if len(entries) > 1 else "A",
                "type": "protein",
                "sequence": seq,
                "count": 1,
                "source_chains": [src_id],
            }
        else:
            seq_counts[seq]["count"] += 1

    chains = list(seq_counts.values())

    # Ligands
    ligand_entities = []
    if smi.exists():
        parsed = parse_ligands_smi(smi)
        for i, L in enumerate(parsed, start=1):
            # Name is commonly the CCD code (e.g. "PQQ")
            ccd = L["source_name"] if (len(L["source_name"]) in (2, 3)
                                        and L["source_name"].isalnum()
                                        and L["source_name"].isupper()) else None
            ligand_entities.append({
                "id": f"lig_{i:02d}",
                "ccd": ccd,
                "smiles": L["smiles"],
                "count": 1,
                "sdf_paths": [],
                "source_name": L["source_name"],
                "source_task": L["task"],
            })

    placeholder = looks_like_placeholder_pdb(template)

    target_dict = {
        "target_id": target,
        "casp_year": 17,
        "set": None,
        "chains": chains,
        "ligands": ligand_entities,
        "receptor_pdb": "",                      # no crystal ref yet
        "template_pdb": ("template.pdb" if template.exists() and not placeholder
                         else ""),
        "notes": f"CASP17 blind target. "
                 f"{'Template placeholder (all-zero coords).' if placeholder else 'Template available.'} "
                 f"No crystal ref yet — scorer cannot run until organizers release truth.",
    }
    return target_dict


def write_target_dir(target: str, d: dict) -> None:
    dst = ROOT / "casp17" / target
    dst.mkdir(parents=True, exist_ok=True)

    (dst / "target.json").write_text(json.dumps(d, indent=2) + "\n")

    # receptor.fasta — per-chain entries, canonical-style
    with (dst / "receptor.fasta").open("w") as f:
        for c in d["chains"]:
            for i in range(c["count"]):
                suf = f"-{i+1}" if c["count"] > 1 else ""
                f.write(f">{target}|{c['id']}{suf}|type={c['type']}\n")
                f.write(c["sequence"] + "\n")

    # ligands.tsv — stable view
    with (dst / "ligands.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["lig_id", "ccd", "smiles", "count", "source_name", "source_task"])
        for L in d["ligands"]:
            w.writerow([L["id"], L.get("ccd") or "", L["smiles"], L["count"],
                        L.get("source_name", ""), L.get("source_task", "")])

    # Copy template.pdb if non-placeholder
    src_template = SPIDER_TARGETS / target / "template.pdb"
    if d.get("template_pdb") and src_template.exists():
        shutil.copy2(src_template, dst / "template.pdb")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", help="single target id (e.g. T1214)")
    args = ap.parse_args()

    if not SPIDER_TARGETS.exists():
        print(f"error: no CASP17 spider data at {SPIDER_TARGETS}", file=sys.stderr)
        return 2

    if args.target:
        targets = [args.target]
    else:
        targets = sorted(p.name for p in SPIDER_TARGETS.iterdir() if p.is_dir())

    ok = 0
    for t in targets:
        print(f"[casp17] {t}")
        d = build_target(t)
        if d is None:
            continue
        write_target_dir(t, d)
        n_chains = sum(c["count"] for c in d["chains"])
        n_ligs = len(d["ligands"])
        placeholder = "template_pdb" in d and not d["template_pdb"]
        print(f"  wrote casp17/{t}/   chains={n_chains}  ligands={n_ligs}  "
              f"{'(template=placeholder)' if placeholder else '(template=real)'}")
        ok += 1

    print(f"\n{ok}/{len(targets)} targets populated in {ROOT / 'casp17'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
