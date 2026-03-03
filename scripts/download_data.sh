#!/usr/bin/env bash
# download_data.sh — fetch all public datasets used in this study
# Requires: wget (or curl as fallback)
# Total download: ~800 MB
#
# Datasets are organized into two groups:
#   DISCOVERY (gene selection):    GSE14520, GSE126848
#   VALIDATION (parameter-free):   TCGA-LIHC, GSE144269, GSE135251
#
# The discovery and validation cohorts are completely independent.
# No validation data was used during gene selection.

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

echo "════════════════════════════════════════════════════════════"
echo "  DISCOVERY COHORTS (used for gene selection only)"
echo "════════════════════════════════════════════════════════════"

echo ""
echo "=== 1/5  GSE14520 — HBV-HCC microarray (Affymetrix HG-U133A) ==="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz" \
      "$DATA_DIR/GSE14520_series_matrix.txt.gz"

echo ""
echo "=== 2/5  GSE126848 — MASLD/NASH spectrum (RNA-seq counts) ==="
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126848&format=file&file=GSE126848%5FGene%5Fcounts%5Fraw%2Etxt%2Egz" \
      "$DATA_DIR/GSE126848_Gene_counts_raw.txt.gz"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  VALIDATION COHORTS (parameter-free scoring pipeline)"
echo "════════════════════════════════════════════════════════════"

echo ""
echo "=== 3/5  TCGA-LIHC (UCSC Xena) ==="
fetch "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_Hugo_norm_count.gz" \
      "$DATA_DIR/tcga_RSEM_Hugo_norm_count.gz"

echo ""
echo "=== 4/5  GSE144269 — External HCC validation ==="
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144269&format=file" \
      "$DATA_DIR/GSE144269_RAW.tar"
if [ -f "$DATA_DIR/GSE144269_RAW.tar" ]; then
  echo "  Extracting GSE144269..."
  tar xf "$DATA_DIR/GSE144269_RAW.tar" -C "$DATA_DIR"
fi

echo ""
echo "=== 5/5  GSE135251 — Fibrosis spectrum ==="
fetch "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135251&format=file&file=GSE135251%5Fcounts%2Ecsv%2Egz" \
      "$DATA_DIR/GSE135251_counts.csv.gz"

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  Done. Files saved to: $DATA_DIR"
echo ""
echo "  To reproduce the gene selection:"
echo "    python scripts/derive_signature.py --save-intermediates"
echo ""
echo "  To run the validation pipeline:"
echo "    cd analysis && python signature_reference.py"
echo "════════════════════════════════════════════════════════════"
