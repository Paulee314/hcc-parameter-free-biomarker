#!/usr/bin/env bash
# download_data.sh — fetch all public datasets used in the analysis
# Requires: wget (or curl as fallback)
# Total download: ~500 MB

set -euo pipefail
DATA_DIR="$(cd "$(dirname "$0")/.." && pwd)/data/raw"
mkdir -p "$DATA_DIR"

fetch() {
  local url="$1" dest="$2"
  if [ -f "$dest" ]; then echo "  Already exists, skipping $(basename "$dest")"; return; fi
  if command -v wget &>/dev/null; then wget -q --show-progress -O "$dest" "$url"
  elif command -v curl &>/dev/null; then curl -L -o "$dest" "$url"
  else echo "ERROR: need wget or curl"; exit 1; fi
}

echo "=== 1/3  TCGA-LIHC (UCSC Xena) ==="
fetch "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_Hugo_norm_count.gz" \
      "$DATA_DIR/tcga_RSEM_Hugo_norm_count.gz"

echo "=== 2/3  GSE144269 ==="
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144269&format=file" \
      "$DATA_DIR/GSE144269_RAW.tar"
if [ -f "$DATA_DIR/GSE144269_RAW.tar" ]; then
  echo "  Extracting GSE144269..."
  tar xf "$DATA_DIR/GSE144269_RAW.tar" -C "$DATA_DIR"
fi

echo "=== 3/3  GSE135251 ==="
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135251&format=file&file=GSE135251%5Fcounts%2Ecsv%2Egz" \
      "$DATA_DIR/GSE135251_counts.csv.gz"

echo ""
echo "Done. Files saved to: $DATA_DIR"