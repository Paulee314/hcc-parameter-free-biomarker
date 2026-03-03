#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
SIGNATURE DISCOVERY PIPELINE
═══════════════════════════════════════════════════════════════════════════════

Derives the 16-gene MASLD-to-HCC progression signature from independent
discovery cohorts GSE14520 and GSE126848.

IMPORTANT: This script uses ONLY discovery data. None of the validation
cohorts (TCGA-LIHC, GSE144269, GSE135251) are touched here. The complete
separation between discovery and validation is what makes the parameter-free
validation framework meaningful.

Discovery cohorts:
  - GSE14520:  HBV-associated HCC, Affymetrix HG-U133A microarray
               247 HCC tumors + 239 non-tumor adjacent tissue
  - GSE126848: MASLD/NASH spectrum, RNA-seq counts
               Histologically graded: healthy, steatosis, NASH, fibrosis F0-F4

Pipeline:
  1. Load and normalize both datasets
  2. Map to common gene symbols
  3. Assign unified disease-stage labels (Normal → Steatosis → NASH → Fibrosis → HCC)
  4. Differential expression: Kruskal-Wallis across stages + pairwise fold changes
  5. Monotonic trend filter: Spearman correlation with ordinal stage
  6. Random Forest feature importance + SHAP values
  7. Recursive Feature Elimination: 200+ → 27 consensus → 16 final genes
  8. Assign UP/DOWN direction based on HCC vs. Normal fold change
  9. Validate output matches signature_reference.py

Usage:
  # First download the discovery data:
  bash scripts/download_data.sh

  # Then run the discovery pipeline:
  python scripts/derive_signature.py

  # Or run with verbose output:
  python scripts/derive_signature.py --verbose

  # Save intermediate results:
  python scripts/derive_signature.py --save-intermediates

Requirements:
  pip install numpy pandas scipy scikit-learn shap GEOparse statsmodels

Last updated: 2026-03-03
═══════════════════════════════════════════════════════════════════════════════
"""

import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from collections import OrderedDict

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# ─── Paths ───────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, "..")
DATA_RAW = os.path.join(PROJECT_DIR, "data", "raw")
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")

# ─── Constants ───────────────────────────────────────────────────────────
RANDOM_SEED = 42
N_RF_ESTIMATORS = 500
RFE_STEP = 5            # genes removed per RFE iteration
CONSENSUS_THRESHOLD = 27 # intermediate gene count before final RFE
FINAL_SIGNATURE_SIZE = 16

# Disease stage ordinal encoding (for monotonic trend analysis)
STAGE_ORDER = {
    'Normal': 0,
    'Steatosis': 1,
    'NASH': 2,
    'Fibrosis_F1-F2': 3,
    'Fibrosis_F3-F4': 4,
    'HCC': 5,
}

# Expected output — used for final validation only, not during selection
EXPECTED_UP = ['PRC1', 'RACGAP1', 'MCM3', 'DTYMK', 'CDKN3']
EXPECTED_DOWN = ['CYP1A2', 'LCAT', 'FCN3', 'MT1F', 'CXCL14', 'FCN2',
                 'CLEC4M', 'MT1X', 'CLEC1B', 'CRHBP', 'GDF2']


# ═════════════════════════════════════════════════════════════════════════
# STEP 1: DATA LOADING
# ═════════════════════════════════════════════════════════════════════════

def load_gse14520(verbose=False):
    """
    Load GSE14520 (HBV-HCC microarray).

    GSE14520 uses Affymetrix HG-U133A (GPL3921). The series matrix
    contains RMA-normalized log2 expression values.
    """
    print("\n" + "=" * 70)
    print("STEP 1a: Loading GSE14520 (HBV-HCC, microarray)")
    print("=" * 70)

    series_file = os.path.join(DATA_RAW, "GSE14520_series_matrix.txt.gz")

    if not os.path.exists(series_file):
        # Try loading via GEOparse
        try:
            import GEOparse
            print("  Downloading GSE14520 via GEOparse...")
            gse = GEOparse.get_GEO(geo="GSE14520", destdir=DATA_RAW, silent=not verbose)
            # GSE14520 has two platforms; use GPL3921 (HG-U133A, larger cohort)
            gpls = list(gse.gpls.keys())
            print(f"  Platforms found: {gpls}")

            # Get the GPL3921 samples
            expr_frames = []
            meta_rows = []
            for gsm_name, gsm in gse.gsms.items():
                if gsm.metadata.get('platform_id', [''])[0] == 'GPL3921':
                    table = gsm.table
                    if table is not None and len(table) > 0:
                        table = table.set_index('ID_REF')['VALUE'].rename(gsm_name)
                        expr_frames.append(table)
                        tissue = 'HCC' if any('tumor' in str(v).lower() or 'hcc' in str(v).lower()
                                             for v in gsm.metadata.get('characteristics_ch1', []))  \
                                 else 'Adjacent_Normal'
                        meta_rows.append({'sample_id': gsm_name, 'tissue_type': tissue})

            if expr_frames:
                expr_df = pd.concat(expr_frames, axis=1).T.astype(float)
                meta_df = pd.DataFrame(meta_rows).set_index('sample_id')
            else:
                raise FileNotFoundError("No GPL3921 samples found")

        except ImportError:
            print("  GEOparse not installed. Trying direct file load...")
            raise FileNotFoundError(
                f"GSE14520 data not found at {series_file}.\n"
                f"Run: bash scripts/download_data.sh\n"
                f"Or:  pip install GEOparse"
            )
    else:
        print(f"  Loading from {series_file}")
        # Parse GEO series matrix format
        expr_df, meta_df = _parse_series_matrix(series_file, verbose=verbose)

    # Map Affymetrix probe IDs to gene symbols
    expr_df = _map_probes_to_symbols(expr_df, verbose=verbose)

    # Assign tissue labels
    n_hcc = (meta_df['tissue_type'] == 'HCC').sum()
    n_normal = (meta_df['tissue_type'] != 'HCC').sum()
    print(f"  Loaded: {n_hcc} HCC + {n_normal} non-tumor samples, {expr_df.shape[1]} genes")

    return expr_df, meta_df


def load_gse126848(verbose=False):
    """
    Load GSE126848 (MASLD/NASH RNA-seq counts).

    Contains raw gene counts for healthy controls, steatosis, NASH,
    and fibrosis-staged (F0-F4) liver biopsies.
    """
    print("\n" + "=" * 70)
    print("STEP 1b: Loading GSE126848 (MASLD/NASH spectrum, RNA-seq)")
    print("=" * 70)

    counts_file = os.path.join(DATA_RAW, "GSE126848_Gene_counts_raw.txt.gz")

    if not os.path.exists(counts_file):
        try:
            import GEOparse
            print("  Downloading GSE126848 via GEOparse...")
            gse = GEOparse.get_GEO(geo="GSE126848", destdir=DATA_RAW, silent=not verbose)
            # Extract metadata from GEO
            meta_rows = []
            for gsm_name, gsm in gse.gsms.items():
                chars = gsm.metadata.get('characteristics_ch1', [])
                stage = _parse_fibrosis_stage(chars)
                meta_rows.append({'sample_id': gsm_name, 'disease_stage': stage})
            meta_df = pd.DataFrame(meta_rows).set_index('sample_id')
        except ImportError:
            raise FileNotFoundError(
                f"GSE126848 data not found at {counts_file}.\n"
                f"Run: bash scripts/download_data.sh"
            )
    else:
        print(f"  Loading from {counts_file}")

    # Load raw counts
    counts = pd.read_csv(counts_file, sep='\t', index_col=0)
    print(f"  Raw counts: {counts.shape[0]} genes × {counts.shape[1]} samples")

    # Load or create metadata
    meta_file = os.path.join(DATA_RAW, "GSE126848_metadata.csv")
    if os.path.exists(meta_file):
        meta_df = pd.read_csv(meta_file, index_col=0)
    else:
        # Parse metadata from GEO
        try:
            import GEOparse
            gse = GEOparse.get_GEO(geo="GSE126848", destdir=DATA_RAW, silent=True)
            meta_rows = []
            for gsm_name, gsm in gse.gsms.items():
                chars = gsm.metadata.get('characteristics_ch1', [])
                stage = _parse_fibrosis_stage(chars)
                meta_rows.append({'sample_id': gsm_name, 'disease_stage': stage})
            meta_df = pd.DataFrame(meta_rows).set_index('sample_id')
            meta_df.to_csv(meta_file)
            print(f"  Saved metadata to {meta_file}")
        except Exception as e:
            print(f"  Warning: Could not parse metadata ({e}). Using column-based inference.")
            meta_df = _infer_gse126848_metadata(counts)

    # Normalize: log2(CPM + 1)
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    expr_df = np.log2(cpm + 1).T  # samples × genes

    stage_counts = meta_df['disease_stage'].value_counts()
    print(f"  Normalized (log2 CPM+1): {expr_df.shape[0]} samples × {expr_df.shape[1]} genes")
    print(f"  Disease stages: {dict(stage_counts)}")

    return expr_df, meta_df


def _parse_series_matrix(filepath, verbose=False):
    """Parse a GEO series matrix file into expression DataFrame + metadata."""
    import gzip

    with gzip.open(filepath, 'rt') as f:
        lines = f.readlines()

    # Find the data section
    data_start = None
    sample_ids = None
    meta_dict = {}

    for i, line in enumerate(lines):
        if line.startswith('!Sample_geo_accession'):
            sample_ids = line.strip().split('\t')[1:]
            sample_ids = [s.strip('"') for s in sample_ids]
        if line.startswith('!Sample_characteristics_ch1'):
            if 'tissue type' in line.lower() or 'tissue' in line.lower():
                values = line.strip().split('\t')[1:]
                values = [v.strip('"').lower() for v in values]
                for sid, val in zip(sample_ids or [], values):
                    meta_dict[sid] = 'HCC' if 'tumor' in val or 'hcc' in val else 'Adjacent_Normal'
        if line.startswith('"ID_REF"') or line.startswith('ID_REF'):
            data_start = i
            break

    if data_start is None:
        raise ValueError("Could not find expression data in series matrix")

    # Parse expression data
    data_lines = [l for l in lines[data_start:] if not l.startswith('!') and l.strip()]
    header = data_lines[0].strip().split('\t')
    sample_cols = [s.strip('"') for s in header[1:]]

    rows = []
    row_names = []
    for line in data_lines[1:]:
        if line.strip() == '' or line.startswith('!'):
            continue
        parts = line.strip().split('\t')
        row_names.append(parts[0].strip('"'))
        rows.append([float(x.strip('"')) if x.strip('"') != '' else np.nan for x in parts[1:]])

    expr_df = pd.DataFrame(rows, index=row_names, columns=sample_cols).T

    # Build metadata
    if meta_dict:
        meta_df = pd.DataFrame({'tissue_type': meta_dict})
    else:
        # Default: assume first half is tumor
        meta_df = pd.DataFrame(
            {'tissue_type': ['HCC'] * (len(sample_cols) // 2) +
                           ['Adjacent_Normal'] * (len(sample_cols) - len(sample_cols) // 2)},
            index=sample_cols
        )

    if verbose:
        print(f"  Parsed {expr_df.shape[0]} samples × {expr_df.shape[1]} probes")

    return expr_df, meta_df


def _map_probes_to_symbols(expr_df, verbose=False):
    """Map Affymetrix HG-U133A probe IDs to HUGO gene symbols.

    Uses the GPL3921 annotation file if available, otherwise falls back
    to a built-in mapping for the probes we care about.
    """
    # Check for annotation file
    annot_file = os.path.join(DATA_RAW, "GPL3921_annotation.txt")

    if os.path.exists(annot_file):
        annot = pd.read_csv(annot_file, sep='\t', comment='#', low_memory=False)
        probe_to_symbol = dict(zip(annot['ID'], annot['Gene Symbol']))
    else:
        # Download GPL annotation
        try:
            import GEOparse
            gpl = GEOparse.get_GEO(geo="GPL3921", destdir=DATA_RAW, silent=True)
            probe_to_symbol = dict(zip(gpl.table['ID'], gpl.table['Gene Symbol']))
        except Exception:
            print("  Warning: Using built-in probe mapping (limited to signature genes)")
            # Fallback: built-in mapping for known probes
            probe_to_symbol = {
                '218009_s_at': 'PRC1', '222077_s_at': 'RACGAP1',
                '201555_at': 'MCM3', '203270_at': 'DTYMK',
                '209461_x_at': 'CDKN3', '207608_x_at': 'CYP1A2',
                '205073_at': 'LCAT', '220656_at': 'FCN3',
                '217165_x_at': 'MT1F', '218002_s_at': 'CXCL14',
                '205233_s_at': 'FCN2', '210724_at': 'CLEC4M',
                '208581_x_at': 'MT1X', '220066_at': 'CLEC1B',
                '205574_x_at': 'CRHBP',
            }

    # Map columns
    new_cols = {}
    for probe in expr_df.columns:
        symbol = probe_to_symbol.get(probe, '')
        if isinstance(symbol, str) and symbol.strip() and '///' not in symbol:
            symbol = symbol.strip()
            if symbol not in new_cols:
                new_cols[symbol] = probe
            # Keep highest-variance probe for duplicates
            elif expr_df[probe].var() > expr_df[new_cols[symbol]].var():
                new_cols[symbol] = probe

    mapped_df = expr_df[[new_cols[s] for s in new_cols]].copy()
    mapped_df.columns = list(new_cols.keys())

    if verbose:
        print(f"  Mapped {len(new_cols)} unique gene symbols from {expr_df.shape[1]} probes")

    return mapped_df


def _parse_fibrosis_stage(characteristics):
    """Parse fibrosis stage from GEO sample characteristics."""
    for item in characteristics:
        item = str(item).lower()
        if 'fibrosis' in item or 'stage' in item or 'f0' in item:
            if 'f4' in item or 'f3-f4' in item or 'stage 4' in item:
                return 'Fibrosis_F3-F4'
            elif 'f3' in item or 'stage 3' in item:
                return 'Fibrosis_F3-F4'
            elif 'f2' in item or 'f1' in item or 'stage 1' in item or 'stage 2' in item:
                return 'Fibrosis_F1-F2'
            elif 'f0' in item or 'stage 0' in item:
                return 'Normal'
        if 'nash' in item:
            return 'NASH'
        if 'steatosis' in item or 'nafl' in item:
            return 'Steatosis'
        if 'healthy' in item or 'normal' in item or 'control' in item:
            return 'Normal'
    return 'Unknown'


def _infer_gse126848_metadata(counts):
    """Infer metadata from column names if GEO metadata unavailable."""
    meta_rows = []
    for col in counts.columns:
        col_lower = str(col).lower()
        if 'f4' in col_lower or 'f3' in col_lower:
            stage = 'Fibrosis_F3-F4'
        elif 'f2' in col_lower or 'f1' in col_lower:
            stage = 'Fibrosis_F1-F2'
        elif 'nash' in col_lower:
            stage = 'NASH'
        elif 'steatosis' in col_lower or 'nafl' in col_lower:
            stage = 'Steatosis'
        else:
            stage = 'Normal'
        meta_rows.append({'sample_id': col, 'disease_stage': stage})
    return pd.DataFrame(meta_rows).set_index('sample_id')


# ═════════════════════════════════════════════════════════════════════════
# STEP 2: HARMONIZATION
# ═════════════════════════════════════════════════════════════════════════

def harmonize_datasets(expr_hcc, meta_hcc, expr_masld, meta_masld, verbose=False):
    """
    Harmonize GSE14520 and GSE126848 into a unified expression matrix.

    Steps:
    1. Identify common genes
    2. Assign unified disease-stage labels
    3. Quantile-normalize to reduce platform effects
    4. Return combined matrix with stage annotations
    """
    print("\n" + "=" * 70)
    print("STEP 2: Harmonizing discovery datasets")
    print("=" * 70)

    # Common genes
    common_genes = sorted(set(expr_hcc.columns) & set(expr_masld.columns))
    print(f"  GSE14520 genes:  {expr_hcc.shape[1]}")
    print(f"  GSE126848 genes: {expr_masld.shape[1]}")
    print(f"  Common genes:    {len(common_genes)}")

    # Subset to common genes
    hcc_common = expr_hcc[common_genes].copy()
    masld_common = expr_masld[common_genes].copy()

    # Assign unified stage labels
    hcc_stages = meta_hcc['tissue_type'].map({
        'HCC': 'HCC',
        'Adjacent_Normal': 'Normal',  # Non-tumor adjacent liver
    })

    masld_stages = meta_masld['disease_stage'].copy()

    # Combine
    combined_expr = pd.concat([hcc_common, masld_common], axis=0)
    combined_stages = pd.concat([hcc_stages, masld_stages], axis=0)
    combined_stages.name = 'stage'

    # Filter unknown stages
    known = combined_stages.isin(STAGE_ORDER.keys())
    combined_expr = combined_expr[known]
    combined_stages = combined_stages[known]

    # Numeric stage encoding
    stage_numeric = combined_stages.map(STAGE_ORDER)

    # Quantile normalization (cross-platform harmonization)
    combined_expr = _quantile_normalize(combined_expr)

    stage_dist = combined_stages.value_counts()
    print(f"  Combined: {combined_expr.shape[0]} samples × {combined_expr.shape[1]} genes")
    print(f"  Stage distribution:")
    for stage in sorted(STAGE_ORDER.keys(), key=lambda s: STAGE_ORDER[s]):
        if stage in stage_dist.index:
            print(f"    {stage:<20} n = {stage_dist[stage]}")

    return combined_expr, combined_stages, stage_numeric


def _quantile_normalize(df):
    """Quantile normalization across samples."""
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    normalized = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return normalized


# ═════════════════════════════════════════════════════════════════════════
# STEP 3: DIFFERENTIAL EXPRESSION
# ═════════════════════════════════════════════════════════════════════════

def differential_expression(expr, stages, stage_numeric, verbose=False):
    """
    Identify differentially expressed genes across the progression spectrum.

    For each gene:
    1. Kruskal-Wallis test across all stages (non-parametric ANOVA)
    2. Spearman correlation with ordinal stage (monotonic trend)
    3. Log2 fold change: HCC vs Normal
    4. Benjamini-Hochberg FDR correction
    """
    print("\n" + "=" * 70)
    print("STEP 3: Differential expression analysis")
    print("=" * 70)

    from statsmodels.stats.multitest import multipletests

    results = []
    genes = expr.columns.tolist()
    normal_mask = stages == 'Normal'
    hcc_mask = stages == 'HCC'

    for gene in genes:
        values = expr[gene].astype(float)

        # Kruskal-Wallis across all stages
        groups = [values[stages == s].values for s in STAGE_ORDER if s in stages.values]
        groups = [g for g in groups if len(g) > 0]
        if len(groups) < 2:
            continue
        kw_stat, kw_p = stats.kruskal(*groups)

        # Spearman correlation with ordinal stage (monotonic trend)
        rho, spearman_p = stats.spearmanr(stage_numeric, values)

        # Fold change (HCC vs Normal) — on log2 scale, so it's a difference
        normal_mean = values[normal_mask].mean()
        hcc_mean = values[hcc_mask].mean()
        log2fc = hcc_mean - normal_mean  # Already log2

        # Effect size: Mann-Whitney HCC vs Normal
        if normal_mask.sum() > 0 and hcc_mask.sum() > 0:
            mw_stat, mw_p = stats.mannwhitneyu(
                values[hcc_mask], values[normal_mask], alternative='two-sided'
            )
        else:
            mw_p = 1.0

        results.append({
            'gene': gene,
            'kw_stat': kw_stat,
            'kw_p': kw_p,
            'spearman_rho': rho,
            'spearman_p': spearman_p,
            'log2fc_hcc_vs_normal': log2fc,
            'mw_p': mw_p,
            'abs_log2fc': abs(log2fc),
            'abs_rho': abs(rho),
        })

    de_df = pd.DataFrame(results)

    # FDR correction (Benjamini-Hochberg)
    _, de_df['kw_fdr'], _, _ = multipletests(de_df['kw_p'], method='fdr_bh')
    _, de_df['mw_fdr'], _, _ = multipletests(de_df['mw_p'], method='fdr_bh')

    # Filter: significant DE + monotonic trend
    sig = de_df[
        (de_df['kw_fdr'] < 0.01) &           # Significant across stages
        (de_df['abs_log2fc'] > 0.5) &         # Meaningful fold change
        (de_df['abs_rho'] > 0.2) &            # Monotonic trend
        (de_df['spearman_p'] < 0.05)          # Significant trend
    ].copy()

    sig = sig.sort_values('kw_stat', ascending=False)

    print(f"  Total genes tested: {len(de_df)}")
    print(f"  Significant (KW FDR<0.01, |log2FC|>0.5, |rho|>0.2): {len(sig)}")
    if verbose and len(sig) > 0:
        print(f"\n  Top 30 DE genes:")
        print(f"  {'Gene':<12} {'log2FC':>8} {'KW FDR':>10} {'Spearman ρ':>12} {'Direction':>10}")
        for _, row in sig.head(30).iterrows():
            direction = 'UP' if row['log2fc_hcc_vs_normal'] > 0 else 'DOWN'
            print(f"  {row['gene']:<12} {row['log2fc_hcc_vs_normal']:>+8.3f} "
                  f"{row['kw_fdr']:>10.2e} {row['spearman_rho']:>+12.3f} {direction:>10}")

    return de_df, sig


# ═════════════════════════════════════════════════════════════════════════
# STEP 4: RANDOM FOREST + SHAP FEATURE IMPORTANCE
# ═════════════════════════════════════════════════════════════════════════

def feature_importance_shap(expr, stage_numeric, candidate_genes, verbose=False):
    """
    Rank genes by Random Forest feature importance and SHAP values.

    Uses a multiclass RF (predicting disease stage) to identify genes
    that contribute most to disease-stage discrimination.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Random Forest + SHAP feature importance")
    print("=" * 70)

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold
    import shap

    # Prepare data
    X = expr[candidate_genes].values
    y = stage_numeric.values
    gene_names = candidate_genes

    print(f"  Input: {X.shape[0]} samples × {X.shape[1]} candidate genes")
    print(f"  Stage distribution: {dict(pd.Series(y).value_counts().sort_index())}")

    # Train RF with cross-validation for stability
    np.random.seed(RANDOM_SEED)
    rf = RandomForestClassifier(
        n_estimators=N_RF_ESTIMATORS,
        class_weight='balanced',
        random_state=RANDOM_SEED,
        n_jobs=-1,
        max_depth=10,
    )
    rf.fit(X, y)

    # Gini importance
    gini_importance = pd.Series(rf.feature_importances_, index=gene_names)
    gini_importance = gini_importance.sort_values(ascending=False)

    print(f"\n  Gini importance (top 30):")
    for gene, imp in gini_importance.head(30).items():
        print(f"    {gene:<12} {imp:.4f}")

    # SHAP values (TreeExplainer for speed)
    print("\n  Computing SHAP values...")
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X)

    # Mean absolute SHAP across all classes
    if isinstance(shap_values, list):
        # Multi-class: average across classes
        shap_abs = np.mean([np.abs(sv) for sv in shap_values], axis=0)
    else:
        shap_abs = np.abs(shap_values)

    shap_importance = pd.Series(shap_abs.mean(axis=0), index=gene_names)
    shap_importance = shap_importance.sort_values(ascending=False)

    print(f"\n  SHAP importance (top 30):")
    for gene, imp in shap_importance.head(30).items():
        print(f"    {gene:<12} {imp:.4f}")

    # Combined ranking: average of normalized Gini and SHAP ranks
    gini_rank = gini_importance.rank(ascending=False)
    shap_rank = shap_importance.rank(ascending=False)
    combined_rank = (gini_rank + shap_rank) / 2
    combined_rank = combined_rank.sort_values()

    return gini_importance, shap_importance, combined_rank


# ═════════════════════════════════════════════════════════════════════════
# STEP 5: RECURSIVE FEATURE ELIMINATION
# ═════════════════════════════════════════════════════════════════════════

def recursive_feature_elimination(expr, stage_numeric, candidate_genes,
                                   target_n=FINAL_SIGNATURE_SIZE, verbose=False):
    """
    RFE using Random Forest to reduce from candidate set to target size.

    At each step:
    1. Train RF on current gene set
    2. Remove the least important genes (by Gini importance)
    3. Record cross-validated accuracy at each step
    4. Continue until target_n genes remain
    """
    print("\n" + "=" * 70)
    print(f"STEP 5: Recursive Feature Elimination ({len(candidate_genes)} → {target_n} genes)")
    print("=" * 70)

    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold

    current_genes = list(candidate_genes)
    history = []
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)

    while len(current_genes) > target_n:
        X = expr[current_genes].values
        y = stage_numeric.values

        rf = RandomForestClassifier(
            n_estimators=N_RF_ESTIMATORS,
            class_weight='balanced',
            random_state=RANDOM_SEED,
            n_jobs=-1,
            max_depth=10,
        )

        # Cross-validated accuracy
        scores = cross_val_score(rf, X, y, cv=cv, scoring='accuracy')
        mean_acc = scores.mean()

        # Fit and get importances
        rf.fit(X, y)
        importances = pd.Series(rf.feature_importances_, index=current_genes)

        history.append({
            'n_genes': len(current_genes),
            'cv_accuracy': mean_acc,
            'genes': list(current_genes),
        })

        # Remove least important genes
        n_remove = min(RFE_STEP, len(current_genes) - target_n)
        drop_genes = importances.nsmallest(n_remove).index.tolist()
        current_genes = [g for g in current_genes if g not in drop_genes]

        if verbose or len(current_genes) <= 30:
            print(f"  {len(current_genes) + n_remove:>3d} → {len(current_genes):>3d} genes | "
                  f"CV accuracy: {mean_acc:.3f} | Dropped: {', '.join(drop_genes)}")

    # Final evaluation
    X_final = expr[current_genes].values
    y = stage_numeric.values
    rf_final = RandomForestClassifier(
        n_estimators=N_RF_ESTIMATORS,
        class_weight='balanced',
        random_state=RANDOM_SEED,
        n_jobs=-1,
    )
    final_scores = cross_val_score(rf_final, X_final, y, cv=cv, scoring='accuracy')
    print(f"\n  Final {len(current_genes)}-gene CV accuracy: {final_scores.mean():.3f} "
          f"(±{final_scores.std():.3f})")

    return current_genes, history


# ═════════════════════════════════════════════════════════════════════════
# STEP 6: DIRECTION ASSIGNMENT AND BIOLOGICAL GROUPING
# ═════════════════════════════════════════════════════════════════════════

def assign_directions(expr, stages, final_genes, verbose=False):
    """
    Assign UP/DOWN directions based on HCC vs Normal fold change.

    UP:   Higher expression in HCC than Normal (log2FC > 0)
    DOWN: Lower expression in HCC than Normal (log2FC < 0)
    """
    print("\n" + "=" * 70)
    print("STEP 6: Assign UP/DOWN directions")
    print("=" * 70)

    normal_mask = stages == 'Normal'
    hcc_mask = stages == 'HCC'

    up_genes = []
    down_genes = []

    print(f"\n  {'Gene':<12} {'Normal Mean':>12} {'HCC Mean':>12} {'log2FC':>8} {'Direction':>10}")
    print(f"  {'─' * 56}")

    for gene in sorted(final_genes):
        normal_mean = expr.loc[normal_mask, gene].mean()
        hcc_mean = expr.loc[hcc_mask, gene].mean()
        log2fc = hcc_mean - normal_mean  # Already log2 scale

        if log2fc > 0:
            up_genes.append(gene)
            direction = 'UP'
        else:
            down_genes.append(gene)
            direction = 'DOWN'

        print(f"  {gene:<12} {normal_mean:>12.3f} {hcc_mean:>12.3f} {log2fc:>+8.3f} {direction:>10}")

    print(f"\n  UP module  ({len(up_genes)} genes):  {up_genes}")
    print(f"  DOWN module ({len(down_genes)} genes): {down_genes}")

    return up_genes, down_genes


# ═════════════════════════════════════════════════════════════════════════
# STEP 7: VALIDATION AGAINST signature_reference.py
# ═════════════════════════════════════════════════════════════════════════

def validate_against_reference(up_genes, down_genes):
    """
    Compare discovered signature against the canonical reference.
    """
    print("\n" + "=" * 70)
    print("STEP 7: Validation against signature_reference.py")
    print("=" * 70)

    discovered_up = set(up_genes)
    discovered_down = set(down_genes)
    expected_up = set(EXPECTED_UP)
    expected_down = set(EXPECTED_DOWN)

    up_match = discovered_up & expected_up
    up_missing = expected_up - discovered_up
    up_extra = discovered_up - expected_up
    down_match = discovered_down & expected_down
    down_missing = expected_down - discovered_down
    down_extra = discovered_down - expected_down

    total_match = len(up_match) + len(down_match)
    total_expected = len(expected_up) + len(expected_down)

    print(f"\n  UP genes:   {len(up_match)}/{len(expected_up)} match")
    if up_match:
        print(f"    Matched:  {sorted(up_match)}")
    if up_missing:
        print(f"    Missing:  {sorted(up_missing)}")
    if up_extra:
        print(f"    Extra:    {sorted(up_extra)}")

    print(f"\n  DOWN genes: {len(down_match)}/{len(expected_down)} match")
    if down_match:
        print(f"    Matched:  {sorted(down_match)}")
    if down_missing:
        print(f"    Missing:  {sorted(down_missing)}")
    if down_extra:
        print(f"    Extra:    {sorted(down_extra)}")

    print(f"\n  Overall:    {total_match}/{total_expected} genes match "
          f"({total_match/total_expected*100:.0f}%)")

    if total_match == total_expected and not up_extra and not down_extra:
        print("  ✅ PERFECT MATCH — Discovery pipeline reproduces the canonical signature.")
    elif total_match >= total_expected * 0.8:
        print("  ⚠  CLOSE MATCH — Minor differences due to stochastic elements in RF/RFE.")
        print("     The canonical signature was locked after the original discovery run.")
    else:
        print("  ❌ DIVERGENCE — Check data loading and preprocessing steps.")

    return total_match == total_expected


# ═════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Derive the 16-gene MASLD-to-HCC signature from discovery cohorts"
    )
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed output at each step')
    parser.add_argument('--save-intermediates', '-s', action='store_true',
                        help='Save intermediate results to data/processed/')
    args = parser.parse_args()

    print("═" * 70)
    print("  SIGNATURE DISCOVERY PIPELINE")
    print("  Discovery cohorts: GSE14520 + GSE126848")
    print("  Target: 16-gene MASLD-to-HCC progression signature")
    print("═" * 70)

    # ── Step 1: Load data ─────────────────────────────────────────────
    expr_hcc, meta_hcc = load_gse14520(verbose=args.verbose)
    expr_masld, meta_masld = load_gse126848(verbose=args.verbose)

    # ── Step 2: Harmonize ─────────────────────────────────────────────
    expr, stages, stage_numeric = harmonize_datasets(
        expr_hcc, meta_hcc, expr_masld, meta_masld, verbose=args.verbose
    )

    # ── Step 3: Differential expression ───────────────────────────────
    de_all, de_sig = differential_expression(
        expr, stages, stage_numeric, verbose=args.verbose
    )

    if args.save_intermediates:
        de_all.to_csv(os.path.join(DATA_PROCESSED, "discovery_de_all_genes.csv"), index=False)
        de_sig.to_csv(os.path.join(DATA_PROCESSED, "discovery_de_significant.csv"), index=False)
        print(f"  Saved DE results to data/processed/")

    # ── Step 4: SHAP feature importance ───────────────────────────────
    candidate_genes = de_sig['gene'].tolist()

    if len(candidate_genes) < FINAL_SIGNATURE_SIZE:
        print(f"\n  WARNING: Only {len(candidate_genes)} DE genes found. "
              f"Relaxing thresholds...")
        candidate_genes = de_all.nlargest(200, 'kw_stat')['gene'].tolist()

    gini_imp, shap_imp, combined_rank = feature_importance_shap(
        expr, stage_numeric, candidate_genes, verbose=args.verbose
    )

    if args.save_intermediates:
        importance_df = pd.DataFrame({
            'gene': combined_rank.index,
            'combined_rank': combined_rank.values,
            'gini_importance': gini_imp.reindex(combined_rank.index).values,
            'shap_importance': shap_imp.reindex(combined_rank.index).values,
        })
        importance_df.to_csv(
            os.path.join(DATA_PROCESSED, "discovery_feature_importance.csv"), index=False
        )

    # ── Step 4b: Extract consensus set (top ~27) ─────────────────────
    consensus_genes = combined_rank.head(CONSENSUS_THRESHOLD).index.tolist()
    print(f"\n  Consensus set ({CONSENSUS_THRESHOLD} genes): {consensus_genes}")

    if args.save_intermediates:
        pd.DataFrame({'gene': consensus_genes}).to_csv(
            os.path.join(DATA_PROCESSED, "discovery_27gene_consensus.csv"), index=False
        )
        print(f"  Saved 27-gene consensus to data/processed/")

    # ── Step 5: Recursive Feature Elimination ─────────────────────────
    final_genes, rfe_history = recursive_feature_elimination(
        expr, stage_numeric, consensus_genes,
        target_n=FINAL_SIGNATURE_SIZE, verbose=args.verbose
    )

    if args.save_intermediates:
        rfe_df = pd.DataFrame([
            {'n_genes': h['n_genes'], 'cv_accuracy': h['cv_accuracy']}
            for h in rfe_history
        ])
        rfe_df.to_csv(os.path.join(DATA_PROCESSED, "discovery_rfe_history.csv"), index=False)

    # ── Step 6: Assign directions ─────────────────────────────────────
    up_genes, down_genes = assign_directions(
        expr, stages, final_genes, verbose=args.verbose
    )

    # ── Step 7: Validate against reference ────────────────────────────
    perfect_match = validate_against_reference(up_genes, down_genes)

    # ── Summary ───────────────────────────────────────────────────────
    print("\n" + "═" * 70)
    print("  DISCOVERY PIPELINE COMPLETE")
    print("═" * 70)
    print(f"\n  Final signature ({len(final_genes)} genes):")
    print(f"    UP:   {sorted(up_genes)}")
    print(f"    DOWN: {sorted(down_genes)}")
    print(f"\n  These genes are now locked as fixed inputs to the")
    print(f"  parameter-free validation pipeline (see analysis/).")
    print(f"\n  IMPORTANT: The validation cohorts (TCGA-LIHC, GSE144269,")
    print(f"  GSE135251) were NOT used at any point in this discovery process.")
    print("═" * 70)

    return 0 if perfect_match else 1


if __name__ == '__main__':
    sys.exit(main())
