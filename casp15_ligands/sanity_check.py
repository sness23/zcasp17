#!/usr/bin/env python3
"""Sanity-check chai-lab fold results.

Per target, verifies:
  - 5 CIFs + 5 score npzs + msa_depth.pdf present
  - All 5 CIFs parse cleanly
  - Protein chain sequences in the predicted CIF match the input FASTA
  - Each input ligand is present in the output (by atom count)
  - CA-CA bond lengths are physical (3.0–4.5 Å for adjacent residues)
  - pLDDT distribution per model (chai-lab writes pLDDT to B-factor)
  - Inter-model CA RMSD (sample diversity)
  - Score consistency across 5 models
  - No inter-chain clashes flag

Usage:
    python3 sanity_check.py results_A40 fastas
"""
import argparse
import json
import sys
from pathlib import Path
from statistics import mean, stdev

import gemmi
import numpy as np


# ---------- helpers ---------- #

def parse_input_fasta(fasta_path: Path):
    """Return ({chain_name -> protein_seq}, [ligand_smiles])."""
    proteins = {}
    ligands = []
    name, kind, buf = None, None, []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if name and kind == "protein":
                proteins[name] = "".join(buf)
            elif name and kind == "ligand":
                ligands.append("".join(buf))
            kind, hdr = line[1:].split("|", 1)
            name = hdr.split("name=")[-1].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name and kind == "protein":
        proteins[name] = "".join(buf)
    elif name and kind == "ligand":
        ligands.append("".join(buf))
    return proteins, ligands


def chain_summary(cif_path: Path):
    """Return per-chain info from a predicted CIF: {chain_id: {seq, n_atoms, kind}}."""
    s = gemmi.read_structure(str(cif_path))
    out = {}
    for chain in s[0]:
        seq = []
        n_atoms = 0
        kind = "ligand"
        for res in chain:
            n_atoms += len(res)
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid() and info.one_letter_code.isalpha():
                seq.append(info.one_letter_code.upper())
                kind = "protein"
            elif info and info.is_nucleic_acid():
                kind = "rna" if len(res.name) == 1 else "dna"
                seq.append(res.name[-1])
        out[chain.name] = {"seq": "".join(seq), "n_atoms": n_atoms, "kind": kind}
    return out


def ca_coords(cif_path: Path):
    """Return dict chain_id -> Nx3 array of CA coordinates (protein only)."""
    s = gemmi.read_structure(str(cif_path))
    out = {}
    for chain in s[0]:
        cas = []
        for res in chain:
            ca = next((a for a in res if a.name == "CA"), None)
            if ca is not None:
                cas.append([ca.pos.x, ca.pos.y, ca.pos.z])
        if cas:
            out[chain.name] = np.asarray(cas)
    return out


def ca_bond_lengths(coords: np.ndarray) -> np.ndarray:
    if len(coords) < 2:
        return np.array([])
    return np.linalg.norm(np.diff(coords, axis=0), axis=1)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD after optimal alignment (Kabsch). Equal lengths required."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Pc_aligned = Pc @ R.T
    return float(np.sqrt(((Pc_aligned - Qc) ** 2).sum(axis=1).mean()))


def plddt_stats(cif_path: Path):
    """Return (mean, min, frac>=70, frac>=90) for protein CA pLDDT (B-factor column)."""
    s = gemmi.read_structure(str(cif_path))
    vals = []
    for chain in s[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if not (info and info.is_amino_acid()):
                continue
            ca = next((a for a in res if a.name == "CA"), None)
            if ca is not None:
                vals.append(ca.b_iso)
    if not vals:
        return None
    arr = np.asarray(vals)
    return {
        "mean": float(arr.mean()),
        "min": float(arr.min()),
        "frac_ge_70": float((arr >= 70).mean()),
        "frac_ge_90": float((arr >= 90).mean()),
        "n_residues": int(len(arr)),
    }


# ---------- per-target check ---------- #

def check_target(target_dir: Path, fasta_dir: Path) -> dict:
    name = target_dir.name
    rep = {"target": name, "issues": [], "ok": True}

    fasta_path = fasta_dir / f"{name}.fasta"
    if not fasta_path.exists():
        rep["issues"].append(f"no input FASTA at {fasta_path}")
        rep["ok"] = False
        return rep

    in_proteins, in_ligands = parse_input_fasta(fasta_path)
    rep["n_input_protein_chains"] = len(in_proteins)
    rep["n_input_ligands"] = len(in_ligands)
    rep["input_total_protein_residues"] = sum(len(s) for s in in_proteins.values())

    # File presence
    cifs = sorted(target_dir.glob("pred.model_idx_*.cif"))
    npzs = sorted(target_dir.glob("scores.model_idx_*.npz"))
    msa_pdf = target_dir / "msa_depth.pdf"
    rep["n_cifs"] = len(cifs)
    rep["n_npzs"] = len(npzs)
    rep["msa_depth_pdf"] = msa_pdf.exists()
    if len(cifs) != 5:
        rep["issues"].append(f"expected 5 CIFs, found {len(cifs)}")
        rep["ok"] = False
    if len(npzs) != 5:
        rep["issues"].append(f"expected 5 npzs, found {len(npzs)}")
        rep["ok"] = False
    if not msa_pdf.exists():
        rep["issues"].append("missing msa_depth.pdf")

    if not cifs:
        return rep

    # Parse all CIFs (catches malformed output)
    parsed = {}
    for cif in cifs:
        try:
            parsed[cif.name] = chain_summary(cif)
        except Exception as e:
            rep["issues"].append(f"parse failed for {cif.name}: {e}")
            rep["ok"] = False
    if not parsed:
        return rep

    # Sequence integrity vs input — check first model only (all 5 should match)
    out0 = parsed[cifs[0].name]
    pred_protein_seqs = {cid: info["seq"] for cid, info in out0.items() if info["kind"] == "protein"}
    in_protein_seq_list = sorted(in_proteins.values())
    pred_protein_seq_list = sorted(pred_protein_seqs.values())
    if in_protein_seq_list != pred_protein_seq_list:
        rep["issues"].append(
            f"protein sequence mismatch: input has chain lens "
            f"{[len(s) for s in in_protein_seq_list]} vs predicted "
            f"{[len(s) for s in pred_protein_seq_list]}"
        )
        rep["ok"] = False

    # Ligand presence — count non-protein/non-NA chains
    pred_lig_chains = [info for info in out0.values() if info["kind"] == "ligand"]
    rep["n_pred_ligand_chains"] = len(pred_lig_chains)
    if len(pred_lig_chains) != len(in_ligands):
        rep["issues"].append(
            f"ligand count mismatch: input {len(in_ligands)} vs predicted {len(pred_lig_chains)}"
        )
        rep["ok"] = False

    # Per-model atom counts must match (no truncated models)
    atom_counts = {cif.name: sum(c["n_atoms"] for c in info.values()) for cif, info in zip(cifs, parsed.values())}
    if len(set(atom_counts.values())) != 1:
        rep["issues"].append(f"atom counts differ across models: {atom_counts}")
        rep["ok"] = False
    rep["atoms_per_model"] = next(iter(atom_counts.values()))

    # CA-CA bond length sanity (per-chain, model 0)
    bond_issues = []
    for cid, info in out0.items():
        if info["kind"] != "protein":
            continue
        cas = ca_coords(cifs[0])[cid] if cid in ca_coords(cifs[0]) else None
        if cas is None or len(cas) < 2:
            continue
        bls = ca_bond_lengths(cas)
        bad = ((bls < 3.0) | (bls > 4.5)).sum()
        if bad:
            bond_issues.append(f"chain {cid}: {bad}/{len(bls)} CA-CA outside 3.0-4.5Å")
    if bond_issues:
        rep["issues"].append("; ".join(bond_issues))
        # don't flag as not-ok — chains can be discontinuous

    # pLDDT
    plddt = plddt_stats(cifs[0])
    if plddt:
        rep["plddt_mean_m0"] = round(plddt["mean"], 1)
        rep["plddt_min_m0"]  = round(plddt["min"], 1)
        rep["plddt_frac_ge_70_m0"] = round(plddt["frac_ge_70"], 3)
        rep["plddt_frac_ge_90_m0"] = round(plddt["frac_ge_90"], 3)
    else:
        rep["plddt_mean_m0"] = None

    # Inter-model CA RMSD — average pairwise across the 5 models, all protein chains concatenated
    rmsds = []
    coords_per_model = []
    for cif in cifs:
        cm = ca_coords(cif)
        # concatenate all protein chain CA's in a stable chain-name order
        coords_per_model.append(np.concatenate([cm[c] for c in sorted(cm)], axis=0))
    if len({len(c) for c in coords_per_model}) == 1 and len(coords_per_model[0]) > 0:
        for i in range(len(cifs)):
            for j in range(i + 1, len(cifs)):
                rmsds.append(kabsch_rmsd(coords_per_model[i], coords_per_model[j]))
    rep["inter_model_rmsd_mean"] = round(mean(rmsds), 2) if rmsds else None
    rep["inter_model_rmsd_max"]  = round(max(rmsds), 2) if rmsds else None

    # Score consistency
    aggs, ptms, iptms, clashes = [], [], [], []
    for npz in npzs:
        d = np.load(npz)
        aggs.append(float(d["aggregate_score"][0]))
        ptms.append(float(d["ptm"][0]))
        iptms.append(float(d["iptm"][0]))
        clashes.append(bool(d["has_inter_chain_clashes"][0]))
    rep["agg_mean"] = round(mean(aggs), 3)
    rep["agg_std"]  = round(stdev(aggs), 3) if len(aggs) > 1 else 0.0
    rep["ptm_mean"] = round(mean(ptms), 3)
    rep["iptm_mean"] = round(mean(iptms), 3)
    rep["any_clash"] = any(clashes)
    if any(clashes):
        rep["issues"].append(f"inter-chain clash in models {[i for i,c in enumerate(clashes) if c]}")
    return rep


# ---------- output ---------- #

def fmt_md(reports: list[dict]) -> str:
    cols = [
        ("target",       "Target"),
        ("ok",           "OK"),
        ("n_cifs",       "CIFs"),
        ("atoms_per_model", "Atoms"),
        ("n_pred_ligand_chains", "Lig"),
        ("plddt_mean_m0", "pLDDT μ"),
        ("plddt_frac_ge_70_m0", "≥70"),
        ("plddt_frac_ge_90_m0", "≥90"),
        ("inter_model_rmsd_mean", "RMSD μ"),
        ("inter_model_rmsd_max",  "RMSD max"),
        ("agg_mean",     "agg μ"),
        ("agg_std",      "agg σ"),
    ]
    lines = ["| " + " | ".join(c[1] for c in cols) + " | Issues |",
             "|" + "|".join(["---"] * (len(cols) + 1)) + "|"]
    for r in reports:
        cells = []
        for k, _ in cols:
            v = r.get(k)
            if isinstance(v, bool):
                cells.append("✓" if v else "✗")
            elif v is None:
                cells.append("—")
            else:
                cells.append(str(v))
        issues = "; ".join(r.get("issues", [])) or "—"
        lines.append("| " + " | ".join(cells) + f" | {issues} |")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("fasta_dir", type=Path, nargs="?", default=Path("fastas"))
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    target_dirs = sorted(p for p in args.results_dir.iterdir()
                         if p.is_dir() and not p.name.startswith("_"))
    reports = [check_target(d, args.fasta_dir) for d in target_dirs]
    reports.sort(key=lambda r: -r.get("agg_mean", 0))

    if args.json:
        print(json.dumps(reports, indent=2))
    else:
        print(fmt_md(reports))
        bad = [r for r in reports if not r["ok"]]
        if bad:
            print(f"\n⚠️  {len(bad)} target(s) flagged: {[r['target'] for r in bad]}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
