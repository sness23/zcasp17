# 05 — Reproduction

End-to-end guide to reproduce the CASP16 Protenix-v1.0.0 benchmark
from scratch. Assumes Linux, conda, RunPod or Vast.ai account.

This doc parallels [`../../casp15_ligands_protenix/docs/05_reproduction.md`](../../casp15_ligands_protenix/docs/05_reproduction.md)
with CASP16-specific differences flagged.

## 0. Prerequisites

Source data already unpacked under `~/data/vaults/casp/`:

- `casp16/pharma_ligands/exper_struct/<SET>/<SET>_prepared/<TARGET>/protein_aligned.pdb`
- `casp16/pharma_ligands/exper_struct/<SET>/<SET>_prepared/<TARGET>/ligand_*.pdb`
- `casp16_ligands/<SET>/<TARGET>.tsv` — ligand SMILES for each target
- `casp16_ligands/refs/<T>/` — prepared refs (if absent, rerun
  `../casp16_ligands/prep_references.py`)

## 1. Set up local environment

Same as CASP15 (see
[../../casp15_ligands_protenix/docs/05_reproduction.md §1](../../casp15_ligands_protenix/docs/05_reproduction.md#1-set-up-the-local-environment)).
Shared between both pipelines.

```bash
conda create -n protenix python=3.11 -y
conda activate protenix
pip install protenix
conda install -n protenix -c bioconda hmmer -y
sudo apt-get install -y kalign

conda create -n ost -c conda-forge openstructure -y

export PROTENIX_ROOT_DIR=/workspace/protenix_data
mkdir -p /workspace/protenix_data/{checkpoint,search_database}
```

## 2. Build JSONs

Key difference from CASP15: CASP16 uses SMILES not CCD codes, and
inputs are split (protein PDB + TSV) rather than combined.

```bash
cd ~/data/vaults/casp/casp16_ligands_protenix
python3 build_protenix_json.py --all
# → jsons/<TARGET>.json   (233 files for full set, 44 for the PoC)
```

The builder now (**post-homodimer-fix**) emits per-chain
`proteinChain` entries. Verify on a multi-chain target:

```bash
python3 -c "
import json
d = json.load(open('jsons/L4001.json'))[0]
for s in d['sequences']:
    k = list(s)[0]; v = s[k]
    if 'sequence' in v:
        print(f'{k}: count={v[\"count\"]} len={len(v[\"sequence\"])}')
    else:
        print(f'{k}: {v[\"ligand\"][:50]}')
"
# should show 2 proteinChain entries (chain A 306 + chain B 304), not
# a single 610-aa chimera
```

## 3. Prep locally

```bash
./prep_local.sh manifest_L1L2L4.txt
# 44 targets × ~5-15 min/target of MSA queue (cache hits for Chymase /
# Mpro are more likely since within-set receptors overlap)
```

**Expected wall time:** ~3-5 hr for all 44 (MSA queue-dominated; actual
CPU work is ~10 min total). Run overnight.

MSA server queue behavior is unpredictable; if some targets sit in
`PENDING` for >30 min, see
[prep_phase.md](prep_phase.md) troubleshooting section.

## 4. Stage for GPU pod

```bash
SRC=/home/sness/data/vaults/casp/casp16_ligands_protenix
DST=/workspace/casp16_ligands_protenix
mkdir -p "$DST"
rsync -aL --exclude '_logs' --exclude 'jsons_concat_bug_backup' \
    "$SRC/" "$DST/"
find "$DST/jsons_prepped" -maxdepth 1 -name '*.json' -exec \
    sed -i "s|$SRC|$DST|g" {} \;
```

## 5. Rent GPU and push

```bash
# RunPod L40S or equivalent 48GB
rsync -avzP -e "ssh -p <port>" /workspace/casp16_ligands_protenix/ \
    root@<ip>:/workspace/casp16_ligands_protenix/
```

Or on Vast.ai with an RTX 3090 / 4090:

```bash
scp -P <port> -r /workspace/casp16_ligands_protenix \
    root@<host>:/workspace/
```

## 6. Run on pod

```bash
ssh root@<host> -p <port>
cd /workspace/casp16_ligands_protenix
pip install protenix
apt-get update && apt-get install -y kalign hmmer rsync
export PROTENIX_ROOT_DIR=/workspace/protenix_data
mkdir -p $PROTENIX_ROOT_DIR/{checkpoint,search_database}

./run_batch.sh manifest_L1L2L4.txt
```

**Expected on L40S:** ~2 hours for all 44 targets × 25 samples =
1100 structures. ~$1.80 at $0.86/hr.

Expected on RTX 4090 (24 GB): works for sub-1024-token targets, may
OOM on a few of the larger L4 targets with the post-fix
2-chain representation. Monitor the pod's memory via `nvidia-smi`.

## 7. Pull results

```bash
rsync -avzP -e "ssh -p <port>" \
    root@<host>:/workspace/casp16_ligands_protenix/results/ \
    /workspace/casp16_ligands_protenix/results/
```

## 8. Score

```bash
cd /workspace/casp16_ligands_protenix
python3 score_lddt_pli.py results out/lddt_pli --best-only
column -ts, out/lddt_pli/summary.csv | less -S
```

Expect 3 "scored 0 pair(s)" entries (L1006, L1007, L1008) from OST
graph-isomorphism failures — these are scorer-side, not method-side.
Same targets fail under chai baseline.

## 9. Head-to-head vs chai

Chai's summary CSV is at
`../casp16_ligands/out/lddt_pli_L40S/summary.csv`. Quick comparison:

```python
import csv
paths = {
    "chai":     "../casp16_ligands/out/lddt_pli_L40S/summary.csv",
    "protenix": "out/lddt_pli/summary.csv",
}
def load(p):
    out = {}
    for r in csv.DictReader(open(p)):
        t = r["target"]
        if t not in out:
            out[t] = r
    return out
c = load(paths["chai"])
p = load(paths["protenix"])
targets = sorted(set(c) & set(p))
for t in targets:
    try:
        print(f"{t:8} chai={float(c[t]['lddt_pli']):.2f} prot={float(p[t]['lddt_pli']):.2f}")
    except:
        print(f"{t:8} (one or both unscored)")
```

## 10. Validation spot tests (optional)

If you want to replicate the homodimer-fix validation:

```bash
# Build just L2001 with the new builder
python3 build_protenix_json.py --set L2000 --target L2001

# One-target manifest
echo "L2001" > manifest_spot_L2001.txt

# Rent cheapest 24 GB GPU (RTX 3090 on Vast is ~$0.18/hr)
vastai create instance <OFFER_ID> --image pytorch/pytorch:2.5.1-cuda12.4-cudnn9-devel --disk 40 --ssh

# Push minimal files only
ssh <host> 'mkdir -p /workspace/casp16_ligands_protenix/jsons'
scp -P <port> run_batch.sh manifest_spot_L2001.txt <host>:/workspace/casp16_ligands_protenix/
scp -P <port> jsons/L2001.json <host>:/workspace/casp16_ligands_protenix/jsons/

# Run (prep + pred, ~10 min)
ssh <host>
cd /workspace/casp16_ligands_protenix
pip install protenix
apt-get install -y kalign hmmer rsync
export PROTENIX_ROOT_DIR=/workspace/protenix_data && mkdir -p $PROTENIX_ROOT_DIR/{checkpoint,search_database}
PROTENIX_ROOT_DIR=/workspace/protenix_data ./run_batch.sh manifest_spot_L2001.txt
```

Expect LDDT_pli ≈ 0.75 (was 0.00 pre-fix), RMSD ≈ 1.9 Å (was 50 Å).

## Troubleshooting

Common issues we hit on this benchmark (additional to CASP15
troubleshooting):

| Symptom | Cause | Fix |
|---|---|---|
| L2001 / all-L4 targets scoring 0.00 LDDT_pli with 40+ Å RMSD | Homodimer bug (fixed in current builder; old JSONs may persist) | Rebuild JSONs from source PDB with current `build_protenix_json.py`; verify per-chain entries in the JSON |
| L1006/1007/1008 "scored 0 ligand pair(s)" | OST RDKit graph-iso fails on novel SMILES | Accept (same fate in chai); or implement ref SDF canonicalization via RDKit roundtrip |
| `MISSING JSON` spam when running batch | `jsons/<T>.json` not present even though `jsons_prepped/<T>.json` is | Update `run_batch.sh` to the current version (skips src requirement if prep exists) |
| Slow PENDING on MSA server | Queue overloaded | Wait; or set `MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com MSA_SERVER=colabfold` to shift to public ColabFold |

See also [../../casp15_ligands_protenix/docs/05_reproduction.md](../../casp15_ligands_protenix/docs/05_reproduction.md)
for the shared troubleshooting section.
