#!/usr/bin/env bash
# push_to_pod.sh — canonical push flow for a CASP pipeline to a fresh GPU pod.
#
# Pushes:
#   - The project dir (canonical inputs, prepped, results, scripts)
#     with exclusions for contaminated/stale dirs.
#   - /workspace/protenix_data/ if locally cached (saves 3 GB download + LayerNorm
#     JIT compile on the first target).
#
# Excluded from the project push (always):
#   - results_postfix/       — per-incident scratch dir
#   - jsons_concat_bug_backup/
#   - out/                   — local scoring output
#
# Usage:
#   ./push_to_pod.sh <project_name> <ssh_target>
#
# Examples:
#   ./push_to_pod.sh casp16_ligands_protenix "root@103.196.86.67 -p 19757 -i ~/.ssh/id_ed25519"
#   ./push_to_pod.sh casp15_ligands_protenix "root@1.2.3.4 -p 12345"

set -euo pipefail

PROJECT="${1:-}"
SSH_TARGET="${2:-}"

if [ -z "$PROJECT" ] || [ -z "$SSH_TARGET" ]; then
    cat >&2 <<EOF
usage: $0 <project_name> "<ssh_args>"

Examples:
  $0 casp16_ligands_protenix "root@103.196.86.67 -p 19757 -i ~/.ssh/id_ed25519"
  $0 casp15_ligands_protenix "root@ssh5.vast.ai -p 31312"
EOF
    exit 2
fi

LOCAL_PROJECT="/workspace/$PROJECT"
REMOTE_PROJECT="/workspace/$PROJECT"

if [ ! -d "$LOCAL_PROJECT" ]; then
    echo "error: no local project at $LOCAL_PROJECT" >&2
    exit 3
fi

# Parse ssh target — rsync wants -e "ssh ..." with the -p/-i flags,
# and the ssh host as the destination spec. Simplest: treat the input
# as free-form ssh args split into (host, flags).
# We keep it simple: user passes the whole thing, we pass it to rsync's -e.
HOST_SPEC="$(echo "$SSH_TARGET" | awk '{print $1}')"
SSH_FLAGS="$(echo "$SSH_TARGET" | cut -d' ' -f2-)"

echo "=== push $PROJECT → $HOST_SPEC ==="

# 1) Project tree, with exclusions
rsync -avzP -e "ssh $SSH_FLAGS" \
    --exclude='results_postfix/' \
    --exclude='jsons_concat_bug_backup/' \
    --exclude='out/' \
    "$LOCAL_PROJECT/" \
    "$HOST_SPEC:$REMOTE_PROJECT/"

# 2) protenix_data cache, if present locally
if [ -d /workspace/protenix_data ] && [ -n "$(ls -A /workspace/protenix_data 2>/dev/null)" ]; then
    echo
    echo "=== pushing cached /workspace/protenix_data ==="
    rsync -avzP -e "ssh $SSH_FLAGS" \
        /workspace/protenix_data/ \
        "$HOST_SPEC:/workspace/protenix_data/"
else
    echo
    echo "NOTE: no local /workspace/protenix_data cache. Pod will redownload the"
    echo "      Protenix checkpoint + data cache (~3 GB, ~3-5 min) on first run."
    echo "      After the pod completes, pull it back:"
    echo
    echo "    rsync -avzP -e \"ssh $SSH_FLAGS\" \\"
    echo "        $HOST_SPEC:/workspace/protenix_data/ \\"
    echo "        /workspace/protenix_data/"
fi

echo
echo "=== done ==="
echo
echo "On the pod:"
echo "  ssh $SSH_TARGET"
echo "  cd /workspace/$PROJECT"
echo "  pip install protenix"
echo "  apt-get update && apt-get install -y kalign hmmer rsync"
echo "  export PROTENIX_ROOT_DIR=/workspace/protenix_data"
echo "  PROTENIX_ROOT_DIR=/workspace/protenix_data ./run_batch.sh <MANIFEST>.txt"
