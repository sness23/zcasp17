#!/usr/bin/env bash
# Fetch all chai-lab model weights + assets into a local cache directory.
# After this runs once, you can rsync DEST to any RunPod pod's $CHAI_DOWNLOADS_DIR
# and skip the ~10 GB download (and its IncompleteRead lottery) on every new pod.
#
# Usage:
#   ./fetch_chai_downloads.sh                      # default dest: $HOME/data/chai_downloads
#   ./fetch_chai_downloads.sh /path/to/cache_dir
#
# Re-running is safe: curl --continue-at resumes partial downloads.

set -euo pipefail

DEST="${1:-$HOME/data/chai_downloads}"
BASE="https://chaiassets.com/chai1-inference-depencencies"

mkdir -p "$DEST/models_v2" "$DEST/esm"

# path_in_cache  relative_url
FILES=(
  "conformers_v1.apkl                                   conformers_v1.apkl"
  "models_v2/feature_embedding.pt                       models_v2/feature_embedding.pt"
  "models_v2/bond_loss_input_proj.pt                    models_v2/bond_loss_input_proj.pt"
  "models_v2/token_embedder.pt                          models_v2/token_embedder.pt"
  "models_v2/trunk.pt                                   models_v2/trunk.pt"
  "models_v2/diffusion_module.pt                        models_v2/diffusion_module.pt"
  "models_v2/confidence_head.pt                         models_v2/confidence_head.pt"
  "esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt            esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt"
)

echo "cache dir: $DEST"
echo

for entry in "${FILES[@]}"; do
    read -r rel url_tail <<<"$entry"
    out="$DEST/$rel"
    url="$BASE/$url_tail"
    if [ -f "$out" ] && [ -s "$out" ]; then
        size=$(du -h "$out" | cut -f1)
        echo "[skip] $rel  ($size already present)"
        continue
    fi
    echo "[get ] $rel"
    # -L follow redirects, -C - resume, --fail bail on 4xx/5xx,
    # --retry for transient blips, --retry-delay 5
    curl --fail -L --retry 5 --retry-delay 5 --retry-all-errors \
         --continue-at - -o "$out" "$url"
done

echo
echo "done. total on disk:"
du -sh "$DEST"
