#!/usr/bin/env bash
# prep_local.sh — run `protenix prep` locally so MSA queue wait doesn't
# burn GPU pod time.
#
# Usage:
#   ./prep_local.sh manifest_L1L2L4.txt
#   MSA_SERVER=protenix ./prep_local.sh manifest_L2000.txt
#
# Behavior:
#   - for each target in the manifest, runs `protenix prep` on jsons/<T>.json
#   - writes prepped JSON to jsons_prepped/<T>.json
#   - skips targets that already have a prepped JSON (resumable)
#   - log per target: jsons_prepped/_logs/<T>.log
#
# When this finishes, scp jsons_prepped/ to the GPU pod and run_batch.sh
# will skip prep (since the prepped JSON exists) and go straight to pred.

set -uo pipefail
shopt -s expand_aliases || true

MANIFEST="${1:-}"
if [ -z "$MANIFEST" ] || [ ! -f "$MANIFEST" ]; then
    echo "usage: $0 <manifest.txt>" >&2
    exit 2
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
JSON_DIR="${JSON_DIR:-$HERE/jsons}"
PREP_DIR="${PREP_DIR:-$HERE/jsons_prepped}"
LOG_DIR="${PREP_DIR}/_logs"
MSA_SERVER="${MSA_SERVER:-protenix}"   # parser schema: 'protenix' or 'colabfold'
# To actually hit api.colabfold.com (not just change the parser), also set:
#     export MMSEQS_SERVICE_HOST_URL=https://api.colabfold.com
# default host is https://protenix-server.com/api/msa

mkdir -p "$PREP_DIR" "$LOG_DIR"

if ! command -v protenix >/dev/null 2>&1; then
    cat >&2 <<'EOF'
ERROR: `protenix` CLI not on PATH.

Install once:
    python3 -m venv ~/venvs/protenix && source ~/venvs/protenix/bin/activate
    pip install protenix

External tools required on PATH:
    Linux: apt-get install -y hmmer kalign
    Mac:   brew install hmmer kalign

Set a data dir so the pdb_seqres template DB caches locally:
    export PROTENIX_ROOT_DIR=$HOME/.cache/protenix
    mkdir -p "$PROTENIX_ROOT_DIR"/{checkpoint,search_database}

Then re-run:
    ./prep_local.sh manifest_L1L2L4.txt
EOF
    exit 3
fi

mapfile -t TARGETS < <(grep -vE '^\s*(#|$)' "$MANIFEST")
total=${#TARGETS[@]}

echo "manifest:    $MANIFEST"
echo "json dir:    $JSON_DIR"
echo "prep dir:    $PREP_DIR"
echo "msa server:  $MSA_SERVER"
echo "targets:     $total"
echo

batch_start=$(date +%s)
ok=0; fail=0; skip=0; i=0

for name in "${TARGETS[@]}"; do
    i=$((i+1))
    src_json="$JSON_DIR/$name.json"
    prep_json="$PREP_DIR/$name.json"
    log="$LOG_DIR/$name.log"

    if [ ! -f "$src_json" ]; then
        echo "[$i/$total] $name: MISSING JSON ($src_json)"
        fail=$((fail+1))
        continue
    fi

    if [ -f "$prep_json" ]; then
        echo "[$i/$total] $name: SKIP (already prepped)"
        skip=$((skip+1))
        continue
    fi

    tstart=$(date +%s)
    echo "[$i/$total] $name: START"
    # Stream prep output to both the terminal (indented) and the log.
    # With `set -o pipefail`, the `if protenix ...` check sees prep's exit code.
    if protenix prep \
            --input "$src_json" \
            --out_dir "$PREP_DIR" \
            --msa_server_mode "$MSA_SERVER" 2>&1 \
        | tr '\r' '\n' | tee -a "$log" | sed -u 's/^/    /' ; then
        # protenix prep writes:
        #   $JSON_DIR/<T>-final-updated.json   (the final prepped JSON)
        #   $JSON_DIR/<T>-update-msa.json      (intermediate — delete)
        #   $PREP_DIR/<T>/msa/0/{pairing,non_pairing,hmmsearch}.a3m  (data)
        final_src="$JSON_DIR/${name}-final-updated.json"
        if [ -f "$final_src" ]; then
            mv "$final_src" "$prep_json"
        fi
        rm -f "$JSON_DIR/${name}-update-msa.json"
        if [ -f "$prep_json" ]; then
            tend=$(date +%s)
            echo "  OK ($((tend-tstart))s)"
            ok=$((ok+1))
            # auto-propagate to any remaining targets with a matching protein
            # sequence — avoids re-submitting identical MSA searches
            if [ -x "$HERE/clone_prepped.py" ]; then
                n_cloned=$(python3 "$HERE/clone_prepped.py" --manifest "$MANIFEST" 2>&1 \
                    | grep -oP 'cloned=\K[0-9]+' | head -n1)
                if [ -n "$n_cloned" ] && [ "$n_cloned" -gt 0 ]; then
                    echo "  auto-cloned $n_cloned sibling target(s) with same sequence"
                fi
            fi
        else
            echo "  prep returned 0 but no prepped JSON found — see $log"
            fail=$((fail+1))
        fi
    else
        echo "  FAIL — see $log"
        fail=$((fail+1))
    fi
done

batch_end=$(date +%s)
echo
echo "done: ok=$ok fail=$fail skip=$skip total=$total wall=$((batch_end-batch_start))s"
