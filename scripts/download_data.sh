#!/usr/bin/env bash
# download_data.sh — fetch all public datasets used in the analysis
# Includes both DISCOVERY cohorts (gene selection) and VALIDATION cohorts
# Requires: wget or curl
# Total download: ~800 MB

set -euo pipefail
DATA_DIR="$(cd "$(dirname "$0")/.." && pwd)/data/raw"
mkdir -p "$DATA_DIR"

fetch() {
  local url="$1" dest="$2"
  if [ -f "$dest" ]; then echo "  Already exists, skipping $(basename "$dest")"; return; fi
  echo "  Downloading $(basename "$dest")..."
  if command -v curl &>/dev/null; then
    curl -sS -L --retry 3 -o "$dest" "$url"
  elif command -v wget &>/dev/null; then
    wget -q --show-progress --no-check-certificate -O "$dest" "$url"
  else echo "ERROR: need wget or curl"; exit 1; fi
}

# ============================================================
#  DISCOVERY COHORTS (used for gene selection — Section 2.2)
# ============================================================

echo "=== 1/5  GSE14520 — HBV-HCC discovery cohort (Affymetrix HG-U133A) ==="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz" \
      "$DATA_DIR/GSE14520-GPL3921_series_matrix.txt.gz"

echo "=== 2/5  GSE126848 — MASLD/NASH discovery cohort (RNA-seq) ==="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126848/suppl/GSE126848_Counts.csv.gz" \
      "$DATA_DIR/GSE126848_Counts.csv.gz"

# ============================================================
#  VALIDATION COHORTS (parameter-free scoring only)
# ============================================================

echo "=== 3/5  TCGA-LIHC (UCSC Xena) ==="
fetch "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_Hugo_norm_count.gz" \
      "$DATA_DIR/tcga_RSEM_Hugo_norm_count.gz"

echo "=== 4/5  GSE144269 — External HCC validation ==="
# Use FTP for series matrix (avoids HTTPS redirect issues with GEO)
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144269/matrix/GSE144269_series_matrix.txt.gz" \
      "$DATA_DIR/GSE144269_series_matrix.txt.gz"

echo "=== 5/5  GSE135251 — Fibrosis spectrum validation ==="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135251/suppl/GSE135251_counts.csv.gz" \
      "$DATA_DIR/GSE135251_counts.csv.gz"

echo ""
echo "Done. Files saved to: $DATA_DIR"
echo ""
echo "Discovery cohorts:  GSE14520, GSE126848"
echo "Validation cohorts: TCGA-LIHC, GSE144269, GSE135251"
