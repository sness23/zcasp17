#!/usr/bin/env python3
"""Compare chai-lab predictions to CASP crystal references.

For each target, picks the best-aggregate-score model and computes:
  - Protein CA RMSD after Kabsch alignment (per matched chain + overall)
  - Ligand centroid distance from the nearest crystal ligand (after the
    protein alignment is applied to the predicted structure)
  - Ligand "footprint overlap": number of crystal ligand heavy atoms with
    any predicted-ligand heavy atom within 4 Å

Reference structures come from
  ../targets_ligand/Targets_ligand/<TARGET>_lig.pdb
which CASP releases as the ground-truth protein+ligand for each target.

Limitations:
  - Chain matching is by longest exact substring of the protein sequence
  - Ligand atom-by-atom RMSD is not computed (would need symmetry-aware
    substructure matching); centroid + footprint give a coarse signal

Usage:
    python3 compare_to_pdb.py results_A40
"""
import argparse
import sys
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent
REF_DIR = ROOT / "targets_ligand" / "Targets_ligand"


# ---------- structure loading ---------- #

def load_protein_ca_by_chain(path: Path) -> dict[str, tuple[str, np.ndarray]]:
    """{chain_id: (sequence, Nx3 CA coords)} for protein chains only."""
    s = gemmi.read_structure(str(path))
    out = {}
    for chain in s[0]:
        seq, cas = [], []
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if not (info and info.is_amino_acid() and info.one_letter_code.isalpha()):
                continue
            ca = next((a for a in res if a.name == "CA"), None)
            if ca is None:
                continue
            seq.append(info.one_letter_code.upper())
            cas.append([ca.pos.x, ca.pos.y, ca.pos.z])
        if cas:
            out[chain.name] = ("".join(seq), np.asarray(cas))
    return out


def load_ligand_atoms(path: Path) -> list[tuple[str, np.ndarray]]:
    """[(residue_name, Nx3 heavy atom coords)] for each non-protein, non-NA residue."""
    s = gemmi.read_structure(str(path))
    out = []
    for chain in s[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            if res.name in ("HOH", "WAT", "DOD"):
                continue
            heavy = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"])
            if len(heavy) > 0:
                out.append((res.name, heavy))
    return out


# ---------- alignment ---------- #

def kabsch_transform(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (R, t) so that (P @ R.T + t) best fits Q."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Qc - Pc @ R.T
    return R, t


def rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))


def match_chains(pred: dict, ref: dict) -> list[tuple[str, str, int]]:
    """Pair predicted chain → reference chain by longest exact protein-seq prefix.

    Returns list of (pred_chain, ref_chain, n_matched_residues).
    Greedy: pairs largest seq match first; each ref chain used at most once.
    """
    pairs = []
    used_ref = set()
    pred_items = sorted(pred.items(), key=lambda kv: -len(kv[1][0]))
    for pcid, (pseq, _pcoords) in pred_items:
        best_ref, best_n = None, 0
        for rcid, (rseq, _rcoords) in ref.items():
            if rcid in used_ref:
                continue
            # Find longest common prefix length
            n = 0
            for a, b in zip(pseq, rseq):
                if a != b:
                    break
                n += 1
            if n > best_n:
                best_n, best_ref = n, rcid
        if best_ref and best_n >= 20:
            pairs.append((pcid, best_ref, best_n))
            used_ref.add(best_ref)
    return pairs


def best_model_path(target_dir: Path) -> Path | None:
    """Return path to the highest-aggregate-score CIF."""
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
    return target_dir / f"pred.model_idx_{best_idx}.cif"


# ---------- per-target comparison ---------- #

def compare(target: str, pred_dir: Path) -> dict:
    rep = {"target": target, "issues": []}
    ref_path = REF_DIR / f"{target}_lig.pdb"
    if not ref_path.exists():
        rep["issues"].append(f"no reference at {ref_path}")
        return rep
    pred_path = best_model_path(pred_dir)
    if pred_path is None:
        rep["issues"].append("no scored CIFs")
        return rep
    rep["best_model"] = pred_path.name

    pred_chains = load_protein_ca_by_chain(pred_path)
    ref_chains  = load_protein_ca_by_chain(ref_path)
    pairs = match_chains(pred_chains, ref_chains)
    rep["matched_chains"] = len(pairs)

    # Build aligned coordinate stacks across all matched chains
    P_stack, Q_stack = [], []
    per_chain_rmsd = []
    for pcid, rcid, n in pairs:
        pseq, pcoords = pred_chains[pcid]
        rseq, rcoords = ref_chains[rcid]
        P, Q = pcoords[:n], rcoords[:n]
        # per-chain pre-aligned RMSD using a chain-local Kabsch
        R, t = kabsch_transform(P, Q)
        per_chain_rmsd.append((pcid, rcid, n, rmsd(P @ R.T + t, Q)))
        P_stack.append(P)
        Q_stack.append(Q)
    if not P_stack:
        rep["issues"].append("no chain match")
        return rep

    P_all = np.concatenate(P_stack, axis=0)
    Q_all = np.concatenate(Q_stack, axis=0)
    R, t = kabsch_transform(P_all, Q_all)
    rep["protein_ca_rmsd_overall"] = round(rmsd(P_all @ R.T + t, Q_all), 2)
    rep["protein_ca_rmsd_per_chain"] = [
        f"{pcid}↔{rcid}({n}): {r:.2f}Å" for pcid, rcid, n, r in per_chain_rmsd
    ]
    rep["matched_residues_total"] = sum(n for _, _, n, _ in per_chain_rmsd)

    # Apply the global protein alignment to the predicted ligand atoms,
    # then compare to crystal ligand atoms.
    pred_ligs = load_ligand_atoms(pred_path)
    ref_ligs  = load_ligand_atoms(ref_path)
    rep["n_pred_ligands"] = len(pred_ligs)
    rep["n_ref_ligands"]  = len(ref_ligs)

    if pred_ligs and ref_ligs:
        # Transform every predicted-ligand coordinate into the reference frame
        ref_centroids = np.array([atoms.mean(axis=0) for _, atoms in ref_ligs])
        deltas = []
        overlaps = []
        for pname, patoms in pred_ligs:
            patoms_aligned = patoms @ R.T + t
            pcentroid = patoms_aligned.mean(axis=0)
            # nearest reference ligand by centroid
            dists = np.linalg.norm(ref_centroids - pcentroid, axis=1)
            nearest = int(np.argmin(dists))
            deltas.append(float(dists[nearest]))
            # footprint overlap: # of crystal heavy atoms within 4Å of any predicted heavy atom
            ref_atoms = ref_ligs[nearest][1]
            d_matrix = np.linalg.norm(
                ref_atoms[:, None, :] - patoms_aligned[None, :, :], axis=2
            )
            overlaps.append(int((d_matrix.min(axis=1) < 4.0).sum()))
        rep["lig_centroid_delta_min"] = round(min(deltas), 2)
        rep["lig_centroid_delta_max"] = round(max(deltas), 2)
        rep["lig_centroid_delta_mean"] = round(sum(deltas) / len(deltas), 2)
        rep["lig_footprint_overlap_total"] = sum(overlaps)
        rep["lig_footprint_overlap_per_lig"] = overlaps
    return rep


def fmt_md(reports: list[dict]) -> str:
    cols = [
        ("target", "Target"),
        ("best_model", "Best CIF"),
        ("matched_chains", "Chains"),
        ("matched_residues_total", "Matched res"),
        ("protein_ca_rmsd_overall", "CA RMSD (Å)"),
        ("n_pred_ligands", "Pred lig"),
        ("n_ref_ligands", "Ref lig"),
        ("lig_centroid_delta_mean", "Lig Δ μ (Å)"),
        ("lig_centroid_delta_min", "Lig Δ min"),
        ("lig_footprint_overlap_total", "Footprint ≤4Å"),
    ]
    lines = ["| " + " | ".join(c[1] for c in cols) + " | Issues |",
             "|" + "|".join(["---"] * (len(cols) + 1)) + "|"]
    for r in reports:
        cells = []
        for k, _ in cols:
            v = r.get(k)
            cells.append("—" if v is None else str(v))
        issues = "; ".join(r.get("issues", [])) or "—"
        lines.append("| " + " | ".join(cells) + f" | {issues} |")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="also print per-chain RMSD breakdowns")
    args = ap.parse_args()

    reports = []
    for sub in sorted(args.results_dir.iterdir()):
        if sub.is_dir() and not sub.name.startswith("_"):
            reports.append(compare(sub.name, sub))
    reports.sort(key=lambda r: r.get("protein_ca_rmsd_overall", 999))
    print(fmt_md(reports))
    if args.verbose:
        print()
        for r in reports:
            print(f"## {r['target']}")
            for line in r.get("protein_ca_rmsd_per_chain", []):
                print(f"  - {line}")
            if r.get("lig_footprint_overlap_per_lig"):
                print(f"  - per-ligand footprint overlap: {r['lig_footprint_overlap_per_lig']}")
            print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
