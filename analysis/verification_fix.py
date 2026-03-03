#!/usr/bin/env python3
"""
Fix checks 6 and 7 from the verification pipeline.
Check 6: Single-cell uses 'site' column (Normal vs Tumor), not 'stage'
Check 7: Actually compute mouse directional concordance
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.metrics import roc_auc_score

DATA = "/sessions/beautiful-nifty-allen/mnt/NIH GEO"

UP = ['PRC1', 'RACGAP1', 'MCM3', 'DTYMK', 'CDKN3']
DOWN = ['CYP1A2', 'LCAT', 'FCN3', 'MT1F', 'CXCL14', 'FCN2', 'CLEC4M',
        'MT1X', 'CLEC1B', 'CRHBP', 'GDF2']

print("=" * 80)
print("  CHECK 6 (FIXED): SINGLE-CELL PAIRED TEST")
print("  Using 'site' column: Normal vs Tumor")
print("=" * 80)

sc_counts = pd.read_csv(f'{DATA}/GSE149614_signature_counts.csv', index_col=0)
sc_meta = pd.read_csv(f'{DATA}/GSE149614_metadata.csv', index_col=0)
sc_meta.columns = sc_meta.columns.str.strip().str.strip('\r')
for col in sc_meta.columns:
    if sc_meta[col].dtype == 'object':
        sc_meta[col] = sc_meta[col].str.strip().str.strip('\r')

sc_up = [g for g in UP if g in sc_counts.index]
sc_dn = [g for g in DOWN if g in sc_counts.index]
print(f"Genes in scRNA: up={sc_up}, down={sc_dn}")
print(f"Site values: {sc_meta['site'].value_counts().to_dict()}\n")

for celltype in ['Hepatocyte', 'Endothelial']:
    print(f"\n--- {celltype} ---")
    ct_mask = sc_meta['celltype'] == celltype
    ct_meta = sc_meta[ct_mask]
    ct_ids = [c for c in ct_meta.index if c in sc_counts.columns]
    ct_counts = sc_counts[ct_ids]
    ct_meta_f = ct_meta.loc[ct_ids]

    print(f"  Sites: {ct_meta_f['site'].value_counts().to_dict()}")

    # Patient-level pseudo-bulk
    paired_norm = []
    paired_tum = []
    pair_pats = []

    for pat in sorted(ct_meta_f['patient'].unique()):
        pat_meta = ct_meta_f[ct_meta_f['patient'] == pat]

        # Normal (adjacent) pseudo-bulk
        norm_ids = pat_meta[pat_meta['site'] == 'Normal'].index.tolist()
        tum_ids = pat_meta[pat_meta['site'] == 'Tumor'].index.tolist()

        if len(norm_ids) < 5 or len(tum_ids) < 5:
            continue

        pb_norm = ct_counts[norm_ids].mean(axis=1)
        pb_tum = ct_counts[tum_ids].mean(axis=1)

        score_norm = pb_norm[sc_up].mean() - pb_norm[sc_dn].mean()
        score_tum = pb_tum[sc_up].mean() - pb_tum[sc_dn].mean()

        paired_norm.append(score_norm)
        paired_tum.append(score_tum)
        pair_pats.append(pat)

        print(f"  {pat}: Normal={score_norm:.4f} (n={len(norm_ids)}), "
              f"Tumor={score_tum:.4f} (n={len(tum_ids)}), "
              f"Δ={score_tum-score_norm:+.4f} {'✓' if score_tum > score_norm else '✗'}")

    if len(paired_norm) >= 3:
        stat, pval = wilcoxon(paired_tum, paired_norm)
        n_correct = sum(t > n for t, n in zip(paired_tum, paired_norm))
        print(f"\n  Paired Wilcoxon: n={len(paired_norm)} pairs")
        print(f"  Statistic: {stat:.1f}")
        print(f"  p-value: {pval:.6f}")
        print(f"  Direction: {n_correct}/{len(paired_norm)} patients tumor > normal")
        if pval < 0.05 and n_correct == len(paired_norm):
            print(f"  ✓ PASS: p<0.05 and all patients correct direction")
        elif pval < 0.05:
            print(f"  ~ PARTIAL: p<0.05 but {len(paired_norm)-n_correct} patients wrong")
        else:
            print(f"  ✗ FAIL: p≥0.05")
    else:
        print(f"  Only {len(paired_norm)} paired patients — insufficient")

# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("  CHECK 7 (FIXED): MOUSE DIRECTIONAL CONCORDANCE")
print("=" * 80)

ortho = pd.read_csv(f'{DATA}/mouse_human_orthologs.csv')
ens2sym = pd.read_csv(f'{DATA}/mouse_ensembl_to_symbol.csv')
ortho_map = dict(zip(ortho['human_symbol'], ortho['mouse_symbol']))
mapped = {g: ortho_map[g] for g in UP + DOWN if g in ortho_map}
print(f"Mapped: {len(mapped)}/16 genes\n")

# ── GSE227620 ──
print("--- GSE227620 (DEN/ALIOS NASH-HCC model) ---")
mouse_counts = pd.read_csv(f'{DATA}/GSE227620_counts.csv', index_col=0)
mouse_meta = pd.read_csv(f'{DATA}/GSE227620_metadata.csv', index_col=0)

# Map Ensembl to symbol
ens_map = dict(zip(ens2sym.iloc[:, 0], ens2sym.iloc[:, 1]))
mouse_counts.index = mouse_counts.index.map(lambda x: ens_map.get(x, x))
mouse_counts = mouse_counts[~mouse_counts.index.str.startswith('__')]

# Identify normal liver vs tumor
# From metadata: source_name_ch1 has tissue type, characteristics_ch1 has tissue source
# characteristics_ch1_3 has cell type including "Whole Tissue", "Large Primary Tumour", etc.
# characteristics_ch1_4 has treatment
# We want "Whole Tissue" samples with "Normal" treatment vs tumor samples

cell_type_col = 'characteristics_ch1_3'
treatment_col = 'characteristics_ch1_4'
tissue_col = 'characteristics_ch1'

print(f"Cell types: {mouse_meta[cell_type_col].unique()}")
print(f"Treatments: {mouse_meta[treatment_col].unique()}")
print(f"Tissues: {mouse_meta[tissue_col].unique()}")

# Normal: Whole Tissue + Normal treatment + Liver tissue
normal_mask = (
    (mouse_meta[cell_type_col] == 'cell type: Whole Tissue') &
    (mouse_meta[treatment_col] == 'treatment: Normal') &
    (mouse_meta[tissue_col] == 'source tissue: Liver')
)

# Tumor: Large or Small Primary Tumour
tumor_mask = mouse_meta[cell_type_col].isin([
    'cell type: Large Primary Tumour',
    'cell type: Small Primary Tumour'
])

normal_ids = mouse_meta[normal_mask].index
tumor_ids = mouse_meta[tumor_mask].index
print(f"\nNormal liver: {len(normal_ids)} samples")
print(f"Tumor: {len(tumor_ids)} samples")

if len(normal_ids) > 0 and len(tumor_ids) > 0:
    n_concordant = 0
    n_tested = 0
    print(f"\nGene-level direction concordance:")
    for human_gene, mouse_gene in mapped.items():
        if mouse_gene not in mouse_counts.index:
            print(f"  {human_gene}→{mouse_gene}: NOT FOUND in counts")
            continue

        norm_vals = mouse_counts.loc[mouse_gene, normal_ids].values.astype(float)
        tum_vals = mouse_counts.loc[mouse_gene, tumor_ids].values.astype(float)

        # Log2(count+1) for raw counts
        norm_log = np.log2(norm_vals + 1)
        tum_log = np.log2(tum_vals + 1)

        fc = tum_log.mean() - norm_log.mean()
        expected = 'up' if human_gene in UP else 'down'
        actual = 'up' if fc > 0 else 'down'
        concordant = expected == actual

        _, p = mannwhitneyu(tum_log, norm_log, alternative='two-sided')

        mark = "✓" if concordant else "✗"
        if concordant:
            n_concordant += 1
        n_tested += 1
        print(f"  {human_gene:10s}→{mouse_gene:10s}: expect={expected:4s} actual={actual:4s} "
              f"Δ={fc:+.3f} p={p:.2e} {mark}")

    print(f"\nConcordance: {n_concordant}/{n_tested} genes match expected direction")
    if n_concordant >= 7:
        print(f"✓ PASS: ≥7 genes concordant")
    else:
        print(f"✗ FAIL: <7 genes concordant")

# ── GSE67679 ──
print("\n--- GSE67679 (diet-induced NASH, microarray) ---")
mouse_expr2 = pd.read_csv(f'{DATA}/GSE67679_expression.csv', index_col=0)
mouse_meta2 = pd.read_csv(f'{DATA}/GSE67679_metadata.csv', index_col=0)
probe_map = pd.read_csv(f'{DATA}/GPL6887_probe_to_gene.csv')

probe_dict = dict(zip(probe_map.iloc[:, 0].astype(str), probe_map.iloc[:, 1]))
mouse_expr2.index = mouse_expr2.index.astype(str).map(lambda x: probe_dict.get(x, x))

# Average duplicate genes
mouse_expr2 = mouse_expr2.groupby(mouse_expr2.index).mean()

# Identify groups from metadata
print(f"Meta columns: {list(mouse_meta2.columns)}")
for col in mouse_meta2.columns:
    vals = mouse_meta2[col].unique()
    if len(vals) < 20:
        print(f"  {col}: {vals}")

# Look for diet/tissue columns
# Typical: CD (control diet) liver vs WD (western diet) liver or tumor
diet_col = None
tissue_col2 = None
for col in mouse_meta2.columns:
    vals_str = ' '.join(str(v) for v in mouse_meta2[col].unique())
    if 'CD' in vals_str or 'diet' in vals_str.lower():
        diet_col = col
    if 'liver' in vals_str.lower() or 'tumor' in vals_str.lower():
        tissue_col2 = col

if diet_col:
    print(f"\nUsing diet column: {diet_col}: {mouse_meta2[diet_col].unique()}")
if tissue_col2:
    print(f"Using tissue column: {tissue_col2}: {mouse_meta2[tissue_col2].unique()}")

# Check which mapped genes are in this dataset
found_67679 = {h: m for h, m in mapped.items() if m in mouse_expr2.index}
print(f"\nGenes found in GSE67679: {len(found_67679)}/{len(mapped)}")
for h, m in found_67679.items():
    print(f"  {h} → {m}: ✓")
