#!/usr/bin/env bash
# run_batch.sh — fold every target in a manifest with chai-lab.
#
# Usage (on the pod):
#   ./run_batch.sh manifest_A40.txt
#   FASTA_DIR=/workspace/fastas OUT_ROOT=/workspace/results ./run_batch.sh manifest_A100.txt
#
# Behavior:
#   - reads target names (one per line, # comments allowed) from the manifest
#   - looks for ${FASTA_DIR}/<target>.fasta
#   - folds into ${OUT_ROOT}/<target>/
#   - skips targets that already have 5 CIFs in their out dir (resumable)
#   - writes per-target log to ${OUT_ROOT}/_logs/<target>.log
#   - keeps going on individual failures; prints a summary at the end

set -uo pipefail

MANIFEST="${1:-}"
if [ -z "$MANIFEST" ] || [ ! -f "$MANIFEST" ]; then
    echo "usage: $0 <manifest.txt>" >&2
    exit 2
fi

FASTA_DIR="${FASTA_DIR:-/workspace/fastas}"
MSA_DIR="${MSA_DIR:-/workspace/msas}"   # if present, use precomputed MSAs per target
OUT_ROOT="${OUT_ROOT:-/workspace/results}"
LOG_DIR="${OUT_ROOT}/_logs"
mkdir -p "$LOG_DIR"

# Drop comments and blank lines from the manifest
mapfile -t TARGETS < <(grep -vE '^\s*(#|$)' "$MANIFEST")
total=${#TARGETS[@]}

echo "manifest:   $MANIFEST"
echo "fasta dir:  $FASTA_DIR"
echo "out root:   $OUT_ROOT"
echo "targets:    $total"
echo "GPU:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "  (nvidia-smi unavailable)"
echo

batch_start=$(date +%s)
ok=0; fail=0; skip=0
i=0

for name in "${TARGETS[@]}"; do
    i=$((i+1))
    fasta="$FASTA_DIR/$name.fasta"
    out="$OUT_ROOT/$name"
    log="$LOG_DIR/$name.log"

    if [ ! -f "$fasta" ]; then
        echo "[$i/$total] $name: MISSING FASTA ($fasta)"
        fail=$((fail+1))
        continue
    fi

    # Resume: skip if 5 CIFs already exist
    if [ -d "$out" ] && [ "$(ls "$out"/pred.model_idx_*.cif 2>/dev/null | wc -l)" -eq 5 ]; then
        echo "[$i/$total] $name: SKIP (already done)"
        skip=$((skip+1))
        continue
    fi

    # chai-lab requires the output dir to be empty
    rm -rf "$out"

    # Pick MSA strategy: precomputed dir if available, else live MSA server
    target_msa_dir="$MSA_DIR/$name"
    if [ -d "$target_msa_dir" ] && compgen -G "$target_msa_dir/*.aligned.pqt" >/dev/null; then
        msa_args=(--msa-directory "$target_msa_dir")
        msa_mode="precomputed"
    else
        msa_args=(--use-msa-server)
        msa_mode="live"
    fi

    echo "[$i/$total] $name: FOLDING ($msa_mode MSA)..."
    t0=$(date +%s)
    if chai-lab fold "${msa_args[@]}" "$fasta" "$out" >"$log" 2>&1; then
        dt=$(( $(date +%s) - t0 ))
        # quick aggregate score from model 0
        agg=$(python3 -c "
import numpy as np
try:
    d = np.load('$out/scores.model_idx_0.npz')
    print(f'agg={float(d[\"aggregate_score\"][0]):.3f} ptm={float(d[\"ptm\"][0]):.3f} iptm={float(d[\"iptm\"][0]):.3f}')
except Exception as e:
    print(f'(scores read failed: {e})')
" 2>/dev/null)
        echo "[$i/$total] $name: OK (${dt}s) $agg"
        ok=$((ok+1))
    else
        dt=$(( $(date +%s) - t0 ))
        echo "[$i/$total] $name: FAIL (${dt}s) — see $log"
        # Tail the log so the failure cause is visible inline
        tail -n 5 "$log" | sed 's/^/    | /'
        fail=$((fail+1))
    fi
done

batch_dur=$(( $(date +%s) - batch_start ))
echo
echo "=== batch done in ${batch_dur}s: $ok ok, $fail failed, $skip skipped (of $total) ==="
echo "results in $OUT_ROOT"
