# Parameter-Free Transcriptomic Biomarker Framework for HCC

A parameter-free scoring framework that exposes preprocessing leakage as the dominant failure mode in transcriptomic biomarker studies for hepatocellular carcinoma (HCC).

## Key Finding

Standard ML pipelines for HCC biomarker validation collapse on external data (AUC 0.31) because transductive leakage — fitting normalization parameters to the combined train+test matrix — inflates apparent accuracy. A parameter-free approach using within-cohort Z-scoring with a fixed 16-gene directional composite achieves AUC 0.976–1.000 on the same held-out data, a **0.67+ AUC gap** attributable entirely to preprocessing-stage leakage.

## The 16-Gene Signature

| Module | Genes | Function |
|--------|-------|----------|
| **UP** (5 genes) | PRC1, RACGAP1, MCM3, DTYMK, CDKN3 | Proliferation / cell cycle |
| **DOWN** (11 genes) | CYP1A2, LCAT, FCN3, MT1F, CXCL14, FCN2, CLEC4M, MT1X, CLEC1B, CRHBP, GDF2 | Hepatocyte identity / differentiation |

The DOWN module drives discrimination at the clinical boundary (HCC vs. F4 cirrhosis). Four of five UP genes actually score *higher* in cirrhosis than HCC, making proliferation-only panels unreliable for the screening population.

## Method

For each dataset independently:
```
Z_gene = (X_gene - μ_normal) / σ_normal
Score  = mean(Z_up) - mean(Z_down)
```

- μ and σ are computed **only from that dataset's control samples**
- No parameters transfer between datasets
- No fitted model, no learned boundaries, no threshold optimization
- Gene list and directions are fixed a priori

## Repository Structure
```
├── analysis/
│   ├── signature_reference.py      # Canonical 16-gene signature (single source of truth)
│   ├── tcga_lihc_scoring.py        # TCGA-LIHC primary cohort scoring
│   ├── cirrhosis_vs_hcc_analysis.py # Seven-test cirrhosis battery
│   ├── phase7_full_correction.py   # Final corrected validation pipeline
│   ├── verification_full.py        # End-to-end result verification
│   ├── verification_fix.py         # Verification fixes
│   └── cfrna_validation_suite.py   # cfRNA exploratory analysis
├── data/
│   ├── metadata/                   # Sample metadata for each cohort
│   └── processed/                  # Analysis outputs and result summaries
├── figures/
│   ├── fig1_fibrosis_gradient.png  # Score across fibrosis-to-cancer continuum (Figure 1)
│   ├── fig2_method_comparison.png  # ML vs parameter-free pipeline (Figure 2)
│   └── fig3_pergene_asymmetry.png  # Per-gene AUC asymmetry (Figure 3)
├── scripts/
│   └── download_data.sh            # Downloads all public datasets (~500 MB)
└── paper/
    └── build_paper_v6.js           # Generates the manuscript .docx (requires Node.js + docx)
```

## Datasets

All data are publicly available:

| Dataset | Role | Samples | Source |
|---------|------|---------|--------|
| TCGA-LIHC | Primary cohort | 320 HCC + 50 normal | [UCSC Xena](https://xenabrowser.net/) |
| GSE144269 | External validation | 70 HCC + 70 matched normal | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269) |
| GSE135251 | Fibrosis spectrum | 216 MASLD/NASH biopsies (F0–F4) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251) |

## Quick Start
```bash
# 1. Download all datasets (~500 MB)
bash scripts/download_data.sh

# 2. Run the analysis (requires Python 3.10+)
cd analysis
python signature_reference.py
```

## Building the Paper
```bash
cd paper
npm install docx
node build_paper_v6.js
# Outputs: HCC_Biomarker_Framework_Paper_v6.docx
```

## Requirements

- **Python 3.10+**: NumPy, SciPy, pandas, scikit-learn
- **Node.js 16+**: `docx` package (paper generation only)

## Key Results

- **HCC vs. F4 cirrhosis**: AUC 0.997, 6.16-unit score gap, no overlapping samples
- **Fibrosis gradient**: Monotonic F0→F4→HCC (Spearman ρ = 0.54, p = 1.1 × 10⁻¹⁷)
- **ML baseline on external data**: AUC 0.31 (StandardScaler + RandomForest, fitted to full matrix)
- **Parameter-free on same external data**: AUC 0.976–1.000
- **Covariate independence**: Age/sex AUC ≈ 0.56; Score AUC = 1.000; no improvement from adding demographics

## License

This project is provided for academic and research purposes.
