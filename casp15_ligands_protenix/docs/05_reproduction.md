# 05 — Reproduction

Step-by-step guide to reproduce the CASP15 Protenix results from
scratch. Assumes Linux, conda, and RunPod or Vast.ai account for GPU
time.

## 0. Prerequisites

Data already available under `~/data/vaults/casp/`:

- `casp15_ligands/targets_ligand/Targets_ligand/<T>_lig.pdb` — the
  CASP15 bundle. Originally from predictioncenter.org via the spider
  at `../casp15/`.
- `casp15_ligands/refs/<T>/` — prepared reference pdb + sdf files. If
  absent, rerun `../casp15_ligands/prep_references.py --all`.

## 1. Set up the local environment

One-time, for building JSONs and running `protenix prep` locally:

```bash
# Conda env for Protenix + its deps
conda create -n protenix python=3.11 -y
conda activate protenix
pip install protenix
conda install -n protenix -c bioconda hmmer -y
sudo apt-get install -y kalign

# Conda env for OST scoring (separate — OST needs its own env)
conda create -n ost -c conda-forge openstructure -y

# Verify
which protenix hmmsearch kalign
which ost     # actually needs: ~/anaconda3/envs/ost/bin/ost
```

One-time, set the Protenix cache dir so the ~1 GB checkpoint and
300 MB template DB don't land in `$HOME`:

```bash
echo 'export PROTENIX_ROOT_DIR=/workspace/protenix_data' >> ~/.zshrc
mkdir -p /workspace/protenix_data/{checkpoint,search_database}
```

## 2. Build the input JSONs

From `~/data/vaults/casp/casp15_ligands_protenix/`:

```bash
python3 build_protenix_json.py --manifest manifest_all.txt
# → jsons/<TARGET>.json  (12 files)
```

Each JSON has one `proteinChain` entry per unique chain sequence with
`count` for homo-oligomers, and one `ligand` entry per unique CCD code
also with `count`.

Spot-check one target:

```bash
python3 -c "
import json
d = json.load(open('jsons/T1124.json'))[0]
for s in d['sequences']:
    k = list(s)[0]; v = s[k]
    print(k, v.get('count'), len(v.get('sequence', v.get('ligand', ''))))
"
```

## 3. Run `protenix prep` locally

MSA search + template search for each target. CPU-bound; MSA server
queue dominates wall time (2–30 min per uncached sequence):

```bash
./prep_local.sh manifest_all.txt
# → jsons_prepped/<T>.json with MSA/template paths filled in
# → jsons_prepped/<T>/msa/0/{pairing,non_pairing,hmmsearch}.a3m
```

`prep_local.sh` is resumable (skips prepped targets) and auto-calls
`clone_prepped.py` after each successful prep to short-circuit
redundant MSAs for identical sequences across the manifest.

**Expected wall time:** 30–90 min for all 12, depending on server
queue. Check the per-target log in another terminal:

```bash
tail -f jsons_prepped/_logs/T1124.log
```

## 4. Stage for GPU pod

Copy the tree to `/workspace/` so the absolute paths baked into
prepped JSONs resolve identically on the pod:

```bash
SRC=/home/sness/data/vaults/casp/casp15_ligands_protenix
DST=/workspace/casp15_ligands_protenix
mkdir -p "$DST"
rsync -aL --exclude '_logs' "$SRC/" "$DST/"
# Rewrite absolute paths (idempotent)
find "$DST/jsons_prepped" -maxdepth 1 -name '*.json' -exec \
    sed -i "s|$SRC|$DST|g" {} \;
```

## 5. Rent a GPU and push

### Option A: RunPod

```bash
# via web UI, rent an L40S 48 GB, note ssh string
# e.g. ssh root@<ip> -p <port>

# push
rsync -avzP -e "ssh -p <port>" /workspace/casp15_ligands_protenix/ \
    root@<ip>:/workspace/casp15_ligands_protenix/
```

### Option B: Vast.ai

```bash
vastai search offers "gpu_ram >= 24 rentable=true verified=true \
    reliability > 0.97 compute_cap >= 800" -o 'dph+' | head
vastai create instance <OFFER_ID> \
    --image pytorch/pytorch:2.5.1-cuda12.4-cudnn9-devel \
    --disk 40 --ssh
vastai ssh-url <INSTANCE_ID>   # → ssh root@ssh5.vast.ai -p <port>

# push
scp -P <port> -r /workspace/casp15_ligands_protenix \
    root@<host>:/workspace/
```

## 6. Run on the pod

```bash
ssh root@<host> -p <port>
cd /workspace/casp15_ligands_protenix

# one-time per fresh pod
pip install protenix
apt-get update && apt-get install -y kalign hmmer rsync
export PROTENIX_ROOT_DIR=/workspace/protenix_data
mkdir -p $PROTENIX_ROOT_DIR/{checkpoint,search_database}

# first run downloads ~1 GB checkpoint + compiles LayerNorm CUDA kernel
# (cached for the pod's lifetime)
./run_batch.sh manifest_all.txt
```

**Expected wall time on L40S:** ~108 min for all 12 (31 min A40 set +
77 min A100 set). Cost ~$1.54.

Watch progress in another SSH:

```bash
tail -f "$(ls -t /workspace/casp15_ligands_protenix/results/_logs/*.log | head -1)"
```

## 7. Pull results

```bash
rsync -avzP -e "ssh -p <port>" \
    root@<host>:/workspace/casp15_ligands_protenix/results/ \
    /workspace/casp15_ligands_protenix/results/
```

## 8. Score

Back on the local workstation, activate the OST env (or just use the
full binary path):

```bash
cd /workspace/casp15_ligands_protenix
conda activate ost      # or: export OST_BIN=~/anaconda3/envs/ost/bin/ost

python3 score_lddt_pli.py results out/lddt_pli --best-only
column -ts, out/lddt_pli/summary.csv | less -S
```

## 9. Optional: score all samples

If you want the per-sample distribution (for checking seed variance):

```bash
python3 score_lddt_pli.py results out/lddt_pli     # no --best-only
# → 12 × 25 = 300 rows if all samples have matches
```

## 10. Sync back to the canonical repo

If `/workspace/` is a scratch location and you want the results in the
source tree:

```bash
SRC=/workspace/casp15_ligands_protenix
DST=/home/sness/data/vaults/casp/casp15_ligands_protenix
rsync -av "$SRC/results/" "$DST/results/"
rsync -av "$SRC/out/" "$DST/out/"
# Skip rsyncing jsons_prepped — it's already canonical in DST
```

## 11. Tear down

- RunPod: stop or terminate the pod via the web UI.
- Vast.ai: `vastai destroy instance <INSTANCE_ID>`.

**Before tearing down, consider copying** the protenix checkpoint
cache out of the pod if you want to reproduce without a re-download
(in case ByteDance's URL goes 403 in the future):

```bash
# from local, before destroying the pod
rsync -avzP -e "ssh -p <port>" \
    root@<host>:/workspace/protenix_data/checkpoint/ \
    ~/.cache/protenix/checkpoint/
```

## Troubleshooting

Common failures seen during this work:

| Symptom | Cause | Fix |
|---|---|---|
| `HTTPError: 403 Forbidden` during checkpoint download | Tried to use `protenix-v2`, which is withheld | Set `MODEL=protenix_base_default_v1.0.0` or update the default in `run_batch.sh` |
| `ImportError: FastaBatchedDataset` | `fair-esm` and `esm` 3.x both installed in same env | Use a fresh conda env with only `fair-esm 2.0.0` |
| `PENDING: 0%` for >30 min in prep | MSA server overloaded | Wait, or try `MSA_SERVER=colabfold MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com` |
| `scored 0 ligand pair(s)` in scoring | OST graph isomorphism failed | Try without `--substructure-match`, or accept as scorer tooling issue (see discussion §2.4) |
| `MISSING JSON` in run_batch.sh | Source JSON not present *and* not prepped | Confirm jsons/ and/or jsons_prepped/ are populated on the pod |
| `OOM` on big target | L40S borderline on 2048-token targets | Rerun on H100 80 GB or set smaller `--cycle`/`--step` |

See also [prep_phase.md](prep_phase.md) for a detailed walkthrough of
what happens during Step 3, and [templates.md](templates.md) for the
template-search path.
