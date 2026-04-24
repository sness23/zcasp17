#!/usr/bin/env python3
"""Pre-fetch ColabFold MSAs locally for a chai-lab batch.

Reads protein chains from each FASTA in a manifest, calls the ColabFold
MMseqs2 server (via chai-lab's own wrapper), and writes one .aligned.pqt
file per unique sequence into <out_dir>/<target>/.

The resulting per-target directory is what you pass to:
    chai-lab fold --msa-directory <out_dir>/<target> ...

Usage:
    python3 prefetch_msas.py manifest_A100.txt
    python3 prefetch_msas.py manifest_A100.txt --out msas
    python3 prefetch_msas.py manifest_A100.txt --fasta-dir fastas

Requires chai-lab installed locally (CPU is fine — no model weights touched).
"""
import argparse
import sys
from pathlib import Path

# chai-lab does the heavy lifting; we just orchestrate per-target.
from chai_lab.data.dataset.msas.colabfold import generate_colabfold_msas

ROOT = Path(__file__).resolve().parent
COLABFOLD_URL = "https://api.colabfold.com"


def parse_protein_chains(fasta_path: Path) -> list[str]:
    """Return the list of protein sequences in a chai-lab FASTA."""
    seqs, current_is_protein = [], False
    buf: list[str] = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if current_is_protein and buf:
                seqs.append("".join(buf))
            current_is_protein = line.startswith(">protein")
            buf = []
        elif current_is_protein:
            buf.append(line.strip())
    if current_is_protein and buf:
        seqs.append("".join(buf))
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest", type=Path, help="manifest with one target per line")
    ap.add_argument("--fasta-dir", type=Path, default=ROOT / "fastas")
    ap.add_argument("--out", type=Path, default=ROOT / "msas",
                    help="where to write per-target msa subdirs")
    args = ap.parse_args()

    targets = [
        ln.strip() for ln in args.manifest.read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    print(f"manifest: {args.manifest} ({len(targets)} targets)")

    args.out.mkdir(exist_ok=True)
    for i, target in enumerate(targets, start=1):
        fasta = args.fasta_dir / f"{target}.fasta"
        if not fasta.exists():
            print(f"[{i}/{len(targets)}] {target}: MISSING fasta ({fasta})")
            continue
        target_dir = args.out / target
        if target_dir.exists() and any(target_dir.glob("*.aligned.pqt")):
            print(f"[{i}/{len(targets)}] {target}: SKIP (msas exist)")
            continue
        target_dir.mkdir(exist_ok=True)
        seqs = parse_protein_chains(fasta)
        if not seqs:
            print(f"[{i}/{len(targets)}] {target}: no protein chains, skipping")
            continue
        print(f"[{i}/{len(targets)}] {target}: fetching MSAs for {len(seqs)} chain(s)...")
        try:
            generate_colabfold_msas(
                protein_seqs=seqs,
                msa_dir=target_dir,
                msa_server_url=COLABFOLD_URL,
                search_templates=False,
            )
            n_pqt = len(list(target_dir.glob("*.aligned.pqt")))
            print(f"[{i}/{len(targets)}] {target}: OK ({n_pqt} .aligned.pqt files)")
        except Exception as e:
            print(f"[{i}/{len(targets)}] {target}: FAIL — {e}", file=sys.stderr)

    print(f"\nDone. Per-target MSA dirs in {args.out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
