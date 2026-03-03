#!/usr/bin/env python3
"""
TCGA-LIHC: Method B (leak-free) scoring with 16-gene MASLD-HCC signature.
Stage-stratified analysis (Early vs Late) and overall AUC.
"""

import pandas as pd
import numpy as np
import gzip
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WORK = "/sessions/beautiful-nifty-allen"
OUT = f"{WORK}/mnt/NIH GEO"

# ── 16-gene signature ──
UP_GENES = ['PRC1', 'RACGAP1', 'MCM3', 'DTYMK', 'CDKN3']
DOWN_GENES = ['CYP1A2', 'LCAT', 'FCN3', 'MT1F', 'CXCL14', 'FCN2',
              'CLEC4M', 'MT1X', 'CLEC1B', 'CRHBP', 'GDF2']

# ── Load data ──
print("Loading TCGA-LIHC data...")
with gzip.open(f"{WORK}/TCGA-LIHC.htseq_fpkm.tsv.gz", 'rt') as f:
    expr = pd.read_csv(f, sep='\t', index_col=0)

with gzip.open(f"{WORK}/TCGA-LIHC.GDC_phenotype.tsv.gz", 'rt') as f:
    pheno = pd.read_csv(f, sep='\t', index_col=0)

print(f"Expression: {expr.shape[0]} samples × {expr.shape[1]} genes")
print(f"Phenotype: {pheno.shape}")

# Align
common = expr.index.intersection(pheno.index)
expr = expr.loc[common]
pheno = pheno.loc[common]
print(f"Aligned: {len(common)} samples")

# ── Check gene availability ──
found_up = [g for g in UP_GENES if g in expr.columns]
found_down = [g for g in DOWN_GENES if g in expr.columns]
missing_up = [g for g in UP_GENES if g not in expr.columns]
missing_down = [g for g in DOWN_GENES if g not in expr.columns]
print(f"\nUP genes found: {len(found_up)}/{len(UP_GENES)} {found_up}")
print(f"DOWN genes found: {len(found_down)}/{len(DOWN_GENES)} {found_down}")
if missing_up: print(f"  Missing UP: {missing_up}")
if missing_down: print(f"  Missing DOWN: {missing_down}")

# ── Method B: Z-score using ONLY normal samples ──
normal_mask = pheno['sample_type'] == 'Normal'
tumor_mask = pheno['sample_type'] == 'Tumor'
normal_idx = pheno[normal_mask].index
tumor_idx = pheno[tumor_mask].index
print(f"\nNormal: {len(normal_idx)}, Tumor: {len(tumor_idx)}")

# Z-score each gene using normal mean/sd
up_z = pd.DataFrame(index=common, columns=found_up, dtype=float)
down_z = pd.DataFrame(index=common, columns=found_down, dtype=float)

for gene in found_up:
    ref = expr.loc[normal_idx, gene]
    mu, sd = ref.mean(), ref.std()
    up_z[gene] = (expr[gene] - mu) / sd if sd > 0 else 0

for gene in found_down:
    ref = expr.loc[normal_idx, gene]
    mu, sd = ref.mean(), ref.std()
    down_z[gene] = (expr[gene] - mu) / sd if sd > 0 else 0

# Composite score = mean(UP z-scores) - mean(DOWN z-scores)
up_mean = up_z.mean(axis=1)
down_mean = down_z.mean(axis=1)
scores = up_mean - down_mean

# ── Overall AUC (explicit Tumor vs Normal only) ──
eval_idx = normal_idx.tolist() + tumor_idx.tolist()
labels = pd.Series(0, index=eval_idx)
labels.loc[tumor_idx] = 1
eval_scores = scores.loc[eval_idx]
auc_overall = roc_auc_score(labels, eval_scores)
print(f"\n{'='*60}")
print(f"Overall AUC (Tumor vs Normal): {auc_overall:.3f}")
print(f"{'='*60}")

# ── Stage-stratified analysis ──
print("\n── Stage-Stratified Analysis ──")
stage_map = pheno['stage'].copy()

# Define early (I, II) vs late (III, IV)
early_stages = ['Stage I', 'Stage II']
late_stages = ['Stage III', 'Stage IV']

early_mask = stage_map.isin(early_stages) & tumor_mask
late_mask = stage_map.isin(late_stages) & tumor_mask

early_idx = pheno[early_mask].index
late_idx = pheno[late_mask].index
print(f"Early-stage (I+II): {len(early_idx)}")
print(f"Late-stage (III+IV): {len(late_idx)}")

# AUC: Normal vs Early-stage
early_labels = pd.Series(0, index=normal_idx.tolist() + early_idx.tolist())
early_labels.loc[early_idx] = 1
early_scores = scores.loc[early_labels.index]
auc_early = roc_auc_score(early_labels, early_scores)

# AUC: Normal vs Late-stage
late_labels = pd.Series(0, index=normal_idx.tolist() + late_idx.tolist())
late_labels.loc[late_idx] = 1
late_scores = scores.loc[late_labels.index]
auc_late = roc_auc_score(late_labels, late_scores)

# Individual stages
results = {}
for stage_name in ['Stage I', 'Stage II', 'Stage III', 'Stage IV']:
    s_mask = (stage_map == stage_name) & tumor_mask
    s_idx = pheno[s_mask].index
    if len(s_idx) >= 5:
        s_labels = pd.Series(0, index=normal_idx.tolist() + s_idx.tolist())
        s_labels.loc[s_idx] = 1
        s_scores = scores.loc[s_labels.index]
        s_auc = roc_auc_score(s_labels, s_scores)
        results[stage_name] = {'n': len(s_idx), 'auc': s_auc}
        print(f"  {stage_name}: n={len(s_idx)}, AUC={s_auc:.3f}")

print(f"\n  Early (I+II): n={len(early_idx)}, AUC = {auc_early:.3f}")
print(f"  Late (III+IV): n={len(late_idx)}, AUC = {auc_late:.3f}")

# ── Bootstrap CIs ──
def bootstrap_auc(labels, scores, n_boot=2000):
    aucs = []
    idx = np.arange(len(labels))
    for _ in range(n_boot):
        boot = np.random.choice(idx, size=len(idx), replace=True)
        l, s = labels.iloc[boot], scores.iloc[boot]
        if len(set(l)) < 2: continue
        aucs.append(roc_auc_score(l, s))
    return np.percentile(aucs, [2.5, 97.5])

np.random.seed(42)
ci_overall = bootstrap_auc(labels, eval_scores)
ci_early = bootstrap_auc(early_labels, early_scores)
ci_late = bootstrap_auc(late_labels, late_scores)

print(f"\n── 95% Bootstrap CIs ──")
print(f"  Overall: {auc_overall:.3f} [{ci_overall[0]:.3f}-{ci_overall[1]:.3f}]")
print(f"  Early:   {auc_early:.3f} [{ci_early[0]:.3f}-{ci_early[1]:.3f}]")
print(f"  Late:    {auc_late:.3f} [{ci_late[0]:.3f}-{ci_late[1]:.3f}]")

# ── Score distributions by stage ──
print("\n── Score Distributions ──")
print(f"  Normal:     mean={scores[normal_idx].mean():.3f}, sd={scores[normal_idx].std():.3f}")
print(f"  Early HCC:  mean={scores[early_idx].mean():.3f}, sd={scores[early_idx].std():.3f}")
print(f"  Late HCC:   mean={scores[late_idx].mean():.3f}, sd={scores[late_idx].std():.3f}")
print(f"  All Tumor:  mean={scores[tumor_idx].mean():.3f}, sd={scores[tumor_idx].std():.3f}")

# ── Mann-Whitney U tests ──
u1, p1 = stats.mannwhitneyu(scores[tumor_idx], scores[normal_idx], alternative='greater')
u2, p2 = stats.mannwhitneyu(scores[early_idx], scores[normal_idx], alternative='greater')
u3, p3 = stats.mannwhitneyu(scores[late_idx], scores[normal_idx], alternative='greater')
print(f"\n── Mann-Whitney U ──")
print(f"  All Tumor vs Normal:   p = {p1:.2e}")
print(f"  Early HCC vs Normal:   p = {p2:.2e}")
print(f"  Late HCC vs Normal:    p = {p3:.2e}")

# ── Grade analysis (bonus) ──
print("\n── Tumor Grade Distribution ──")
grade_col = 'tumor_grade'
if grade_col in pheno.columns:
    for grade in sorted(pheno.loc[tumor_mask, grade_col].dropna().unique()):
        g_idx = pheno[(pheno[grade_col] == grade) & tumor_mask].index
        if len(g_idx) >= 5:
            g_labels = pd.Series(0, index=normal_idx.tolist() + g_idx.tolist())
            g_labels.loc[g_idx] = 1
            g_scores = scores.loc[g_labels.index]
            g_auc = roc_auc_score(g_labels, g_scores)
            print(f"  {grade}: n={len(g_idx)}, AUC={g_auc:.3f}, mean_score={scores[g_idx].mean():.3f}")

# ── Demographics audit ──
print("\n── Demographics Audit ──")
if 'gender' in pheno.columns:
    for grp, mask in [('Tumor', tumor_mask), ('Normal', normal_mask)]:
        g = pheno.loc[mask, 'gender'].value_counts()
        print(f"  {grp}: {dict(g)}")

if 'age_at_diagnosis' in pheno.columns:
    for grp, mask in [('Tumor', tumor_mask), ('Normal', normal_mask)]:
        ages = pheno.loc[mask, 'age_at_diagnosis'].dropna()
        if len(ages) > 0:
            print(f"  {grp} age: mean={ages.mean():.1f}, sd={ages.std():.1f}, n={len(ages)}")

# ════════════════════════════════════════════════
# FIGURES
# ════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: ROC curves by stage
ax = axes[0]
fpr0, tpr0, _ = roc_curve(labels, eval_scores)
ax.plot(fpr0, tpr0, 'k-', lw=2, label=f'All Tumor (n={len(tumor_idx)}) AUC={auc_overall:.3f}')

fpr_e, tpr_e, _ = roc_curve(early_labels, early_scores)
ax.plot(fpr_e, tpr_e, 'b-', lw=2, label=f'Early I+II (n={len(early_idx)}) AUC={auc_early:.3f}')

fpr_l, tpr_l, _ = roc_curve(late_labels, late_scores)
ax.plot(fpr_l, tpr_l, 'r-', lw=2, label=f'Late III+IV (n={len(late_idx)}) AUC={auc_late:.3f}')

ax.plot([0,1],[0,1],'--', color='gray', alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('TCGA-LIHC: ROC by Stage\n(Method B, Normal as Reference)', fontsize=13)
ax.legend(loc='lower right', fontsize=10)

# Panel B: Individual stage ROCs
ax = axes[1]
colors = {'Stage I': '#2196F3', 'Stage II': '#4CAF50', 'Stage III': '#FF9800', 'Stage IV': '#F44336'}
for stage_name, info in sorted(results.items()):
    s_mask = (stage_map == stage_name) & tumor_mask
    s_idx = pheno[s_mask].index
    s_labels = pd.Series(0, index=normal_idx.tolist() + s_idx.tolist())
    s_labels.loc[s_idx] = 1
    s_scores = scores.loc[s_labels.index]
    fpr_s, tpr_s, _ = roc_curve(s_labels, s_scores)
    ax.plot(fpr_s, tpr_s, color=colors[stage_name], lw=2,
            label=f'{stage_name} (n={info["n"]}) AUC={info["auc"]:.3f}')

ax.plot([0,1],[0,1],'--', color='gray', alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('TCGA-LIHC: ROC by Individual Stage', fontsize=13)
ax.legend(loc='lower right', fontsize=10)

# Panel C: Score distributions
ax = axes[2]
parts = ax.violinplot(
    [scores[normal_idx].values, scores[early_idx].values, scores[late_idx].values],
    positions=[0, 1, 2], showmeans=True, showmedians=True
)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(['#4CAF50', '#2196F3', '#F44336'][i])
    pc.set_alpha(0.7)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels([f'Normal\n(n={len(normal_idx)})',
                     f'Early I+II\n(n={len(early_idx)})',
                     f'Late III+IV\n(n={len(late_idx)})'])
ax.set_ylabel('Composite Score (Method B)', fontsize=12)
ax.set_title('Score Distributions by Stage Group', fontsize=13)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUT}/fig_tcga_lihc_stage_stratified.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved figure: fig_tcga_lihc_stage_stratified.png")

# ── Save scores ──
score_df = pd.DataFrame({
    'Sample': common,
    'SampleType': pheno.loc[common, 'sample_type'],
    'Stage': pheno.loc[common, 'stage'],
    'Gender': pheno.loc[common, 'gender'] if 'gender' in pheno.columns else 'NA',
    'Score': scores.loc[common].values
})
score_df.to_csv(f"{WORK}/TCGA-LIHC_method_b_scores.csv", index=False)
print(f"Saved scores: TCGA-LIHC_method_b_scores.csv")

print("\n" + "="*60)
print("TCGA-LIHC ANALYSIS COMPLETE")
print("="*60)
