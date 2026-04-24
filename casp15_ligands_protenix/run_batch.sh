#!/usr/bin/env bash
# run_batch.sh — Protenix-v2 PoC batch for CASP16 pharma ligands.
#
# Usage (on the pod):
#   ./run_batch.sh manifest_L1L2L4.txt
#   JSON_DIR=/workspace/jsons OUT_ROOT=/workspace/results ./run_batch.sh manifest_L1L2L4.txt
#
# Per target:
#   1. protenix prep   → MSA (paired+unpaired) + template search, rewrites JSON
#   2. protenix pred   → 5 seeds × 5 samples, protenix-v2, --use_template true
#
# Resumable: skips targets that already have 25 CIFs under their output dir.

set -uo pipefail

MANIFEST="${1:-}"
if [ -z "$MANIFEST" ] || [ ! -f "$MANIFEST" ]; then
    echo "usage: $0 <manifest.txt>" >&2
    exit 2
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
JSON_DIR="${JSON_DIR:-$HERE/jsons}"
PREP_DIR="${PREP_DIR:-$HERE/jsons_prepped}"
OUT_ROOT="${OUT_ROOT:-$HERE/results}"
LOG_DIR="${OUT_ROOT}/_logs"
# protenix-v2 weights are not publicly released (ByteDance issue #295).
# Default to v1.0.0 — also has template/RNA-MSA support and beats AF3.
MODEL="${MODEL:-protenix_base_default_v1.0.0}"
SEEDS="${SEEDS:-101,102,103,104,105}"     # v2 paper recommendation
MSA_SERVER="${MSA_SERVER:-protenix}"       # parser schema: 'protenix' or 'colabfold'
# NOTE: --msa_server_mode only switches the parser.  To actually hit a different
# host, also set MMSEQS_SERVICE_HOST_URL (defaults to https://protenix-server.com/api/msa).
# To use api.colabfold.com:
#     MSA_SERVER=colabfold MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com ./run_batch.sh ...

mkdir -p "$PREP_DIR" "$LOG_DIR"

mapfile -t TARGETS < <(grep -vE '^\s*(#|$)' "$MANIFEST")
total=${#TARGETS[@]}
expected_cifs=$(( $(awk -F, '{print NF}' <<<"$SEEDS") * 5 ))

echo "manifest:     $MANIFEST"
echo "json dir:     $JSON_DIR"
echo "prep dir:     $PREP_DIR"
echo "out root:     $OUT_ROOT"
echo "model:        $MODEL"
echo "seeds:        $SEEDS   (expecting $expected_cifs CIFs/target)"
echo "msa server:   $MSA_SERVER"
echo "targets:      $total"
echo "GPU:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "  (nvidia-smi unavailable)"
echo

batch_start=$(date +%s)
ok=0; fail=0; skip=0
i=0

for name in "${TARGETS[@]}"; do
    i=$((i+1))
    src_json="$JSON_DIR/$name.json"
    prep_json="$PREP_DIR/$name.json"
    out="$OUT_ROOT/$name"
    log="$LOG_DIR/$name.log"

    # Source JSON is only required if prep hasn't been done yet.
    # If jsons_prepped/<T>.json exists we can go straight to pred.
    if [ ! -f "$src_json" ] && [ ! -f "$prep_json" ]; then
        echo "[$i/$total] $name: MISSING JSON ($src_json)"
        fail=$((fail+1))
        continue
    fi

    # Resume: skip if we already have the full sample set
    if [ -d "$out" ]; then
        n_cif=$(find "$out" -name "*_sample_*.cif" 2>/dev/null | wc -l)
        if [ "$n_cif" -eq "$expected_cifs" ]; then
            echo "[$i/$total] $name: SKIP (already $n_cif CIFs)"
            skip=$((skip+1))
            continue
        fi
    fi

    tstart=$(date +%s)
    echo "[$i/$total] $name: START"

    # Step 1: prep (MSA + template). Skip if already prepped.
    if [ ! -f "$prep_json" ]; then
        echo "  prep..."
        if ! protenix prep \
                --input "$src_json" \
                --out_dir "$PREP_DIR" \
                --msa_server_mode "$MSA_SERVER" \
                >> "$log" 2>&1 ; then
            echo "  prep FAILED — see $log"
            fail=$((fail+1))
            continue
        fi
        # protenix prep writes the final JSON as <JSON_DIR>/<T>-final-updated.json
        # (and an intermediate <T>-update-msa.json that we delete).
        final_src="$JSON_DIR/${name}-final-updated.json"
        if [ -f "$final_src" ]; then
            mv "$final_src" "$prep_json"
        fi
        rm -f "$JSON_DIR/${name}-update-msa.json"
    fi

    if [ ! -f "$prep_json" ]; then
        echo "  no prepped JSON produced — see $log"
        fail=$((fail+1))
        continue
    fi

    # Step 2: predict.  protenix-v2 + templates + v2-recommended 5 seeds.
    rm -rf "$out"
    mkdir -p "$out"
    echo "  pred..."
    if protenix pred \
            -i "$prep_json" \
            -o "$out" \
            -s "$SEEDS" \
            -n "$MODEL" \
            --use_template true \
            --use_default_params true \
            >> "$log" 2>&1 ; then
        tend=$(date +%s)
        echo "  OK ($((tend-tstart))s)"
        ok=$((ok+1))
    else
        echo "  pred FAILED — see $log"
        fail=$((fail+1))
    fi
done

batch_end=$(date +%s)
echo
echo "done: ok=$ok fail=$fail skip=$skip total=$total wall=$((batch_end-batch_start))s"
