#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
16-GENE MASLD-TO-HCC PROGRESSION SIGNATURE — CANONICAL REFERENCE
═══════════════════════════════════════════════════════════════════════════════

This file is the SINGLE SOURCE OF TRUTH for the 16-gene signature.
All analysis scripts should import from here. Never hardcode gene lists elsewhere.

Validated across:
  - TCGA-LIHC (STAR FPKM):  AUC = 0.9974, 16/16 concordant, 50 paired samples
  - TCGA-LIHC (Xena HTSeq): AUC = 1.0000, 16/16 concordant, 50 paired samples
  - GSE144269 (RNA-seq):     AUC = 0.9745, 16/16 concordant, 70 paired samples
  - GSE94660  (RNA-seq):     AUC = 1.0000, 16/16 concordant, 21 paired samples
  - GSE14520  (Affy U133A):  14/15 concordant (CLEC1B below microarray detection)

Ensembl IDs verified against GSE144269 ENSG|Symbol pipe-delimited index.
GRCh38 coordinates verified from UCSC/Ensembl and BigWig analysis.

CAUTION: A prior audit script used WRONG Ensembl IDs for 7/9 genes.
         Always use the IDs in THIS file.

Last updated: 2026-03-02
═══════════════════════════════════════════════════════════════════════════════
"""

from collections import OrderedDict

# ═══════════════════════════════════════════════════════════════════════════
# GENE DEFINITIONS — the canonical 16-gene signature
# ═══════════════════════════════════════════════════════════════════════════

SIGNATURE_GENES = OrderedDict({
    # ── UP-REGULATED IN HCC (5 genes) ──────────────────────────────────
    'PRC1': {
        'direction': 'UP',
        'ensembl_id': 'ENSG00000198901',
        'grch38_coords': 'chr15:90981043-91026813',
        'affy_hgu133a_probe': '218009_s_at',
        'tier': 'B',
        'description': 'Protein Regulator of Cytokinesis 1; mitotic spindle midzone protein',
    },
    'RACGAP1': {
        'direction': 'UP',
        'ensembl_id': 'ENSG00000161800',
        'grch38_coords': 'chr12:12524082-12561455',
        'affy_hgu133a_probe': '222077_s_at',
        'tier': 'B',
        'description': 'Rac GTPase Activating Protein 1; cytokinesis regulator',
    },
    'MCM3': {
        'direction': 'UP',
        'ensembl_id': 'ENSG00000112118',
        'grch38_coords': 'chr6:52128554-52149820',
        'affy_hgu133a_probe': '201555_at',
        'tier': 'A',
        'description': 'Minichromosome Maintenance Complex Component 3; DNA replication licensing',
    },
    'DTYMK': {
        'direction': 'UP',
        'ensembl_id': 'ENSG00000168393',
        'grch38_coords': 'chr2:241487843-241499762',
        'affy_hgu133a_probe': '203270_at',
        'tier': 'B',
        'description': 'Deoxythymidylate Kinase; thymidine nucleotide biosynthesis',
    },
    'CDKN3': {
        'direction': 'UP',
        'ensembl_id': 'ENSG00000100526',
        'grch38_coords': 'chr14:54397767-54421082',
        'affy_hgu133a_probe': '209461_x_at',
        'tier': 'A',
        'description': 'Cyclin Dependent Kinase Inhibitor 3; CDK2 dephosphorylation',
    },

    # ── DOWN-REGULATED IN HCC (11 genes) ──────────────────────────────
    'CYP1A2': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000140505',
        'grch38_coords': 'chr15:74748846-74756607',
        'affy_hgu133a_probe': '207608_x_at',
        'tier': 'A+',
        'description': 'Cytochrome P450 1A2; hepatocyte-specific drug metabolism enzyme',
    },
    'LCAT': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000213398',
        'grch38_coords': 'chr16:67905816-67910062',
        'affy_hgu133a_probe': '205073_at',
        'tier': 'A',
        'description': 'Lecithin-Cholesterol Acyltransferase; HDL metabolism',
    },
    'FCN3': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000142748',
        'grch38_coords': 'chr1:27393510-27400458',
        'affy_hgu133a_probe': '220656_at',
        'tier': 'A+',
        'description': 'Ficolin 3; complement lectin pathway, liver-secreted',
    },
    'MT1F': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000198417',
        'grch38_coords': 'chr16:56677847-56679048',
        'affy_hgu133a_probe': '217165_x_at',
        'tier': 'A',
        'description': 'Metallothionein 1F; heavy metal binding, oxidative stress response',
    },
    'CXCL14': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000145824',
        'grch38_coords': 'chr5:135892710-135901548',
        'affy_hgu133a_probe': '218002_s_at',
        'tier': 'A',
        'description': 'C-X-C Motif Chemokine Ligand 14; constitutively expressed chemokine',
    },
    'FCN2': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000160339',
        'grch38_coords': 'chr9:134906949-134915987',
        'affy_hgu133a_probe': '205233_s_at',
        'tier': 'A+',
        'description': 'Ficolin 2; complement lectin pathway, liver-secreted opsonin',
    },
    'CLEC4M': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000104938',
        'grch38_coords': 'chr19:7741931-7750461',
        'affy_hgu133a_probe': '210724_at',
        'tier': 'A+',
        'description': 'C-Type Lectin Domain Family 4 Member M; liver sinusoidal endothelial cell marker',
    },
    'MT1X': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000187193',
        'grch38_coords': 'chr16:56648804-56650100',
        'affy_hgu133a_probe': '208581_x_at',
        'tier': 'A',
        'description': 'Metallothionein 1X; heavy metal detoxification',
    },
    'CLEC1B': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000165682',
        'grch38_coords': 'chr12:10025537-10034697',
        'affy_hgu133a_probe': '220066_at',
        'tier': 'A+',
        'description': 'C-Type Lectin Domain Family 1 Member B; platelet activation receptor',
    },
    'CRHBP': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000145708',
        'grch38_coords': 'chr5:76267138-76280210',
        'affy_hgu133a_probe': '205574_x_at',
        'tier': 'A+',
        'description': 'Corticotropin Releasing Hormone Binding Protein; stress axis modulator',
    },
    'GDF2': {
        'direction': 'DOWN',
        'ensembl_id': 'ENSG00000263761',
        'grch38_coords': 'chr10:46998736-47006675',
        'affy_hgu133a_probe': None,  # Not on HGU133A
        'tier': 'A+',
        'description': 'Growth Differentiation Factor 2 (BMP9); liver sinusoidal homeostasis',
    },
})


# ═══════════════════════════════════════════════════════════════════════════
# CONVENIENCE ACCESSORS
# ═══════════════════════════════════════════════════════════════════════════

ALL_GENES    = list(SIGNATURE_GENES.keys())
UP_GENES     = [g for g, d in SIGNATURE_GENES.items() if d['direction'] == 'UP']
DOWN_GENES   = [g for g, d in SIGNATURE_GENES.items() if d['direction'] == 'DOWN']

ENSEMBL_TO_SYMBOL = {d['ensembl_id']: g for g, d in SIGNATURE_GENES.items()}
SYMBOL_TO_ENSEMBL = {g: d['ensembl_id'] for g, d in SIGNATURE_GENES.items()}

AFFY_TO_SYMBOL = {d['affy_hgu133a_probe']: g for g, d in SIGNATURE_GENES.items()
                  if d['affy_hgu133a_probe'] is not None}

TIER_AP = [g for g, d in SIGNATURE_GENES.items() if d['tier'] == 'A+']
TIER_A  = [g for g, d in SIGNATURE_GENES.items() if d['tier'] == 'A']
TIER_B  = [g for g, d in SIGNATURE_GENES.items() if d['tier'] == 'B']


# ═══════════════════════════════════════════════════════════════════════════
# KNOWN-WRONG ENSEMBL IDs (from the flawed audit script)
# Documented here so we never repeat the mistake
# ═══════════════════════════════════════════════════════════════════════════

WRONG_ENSEMBL_IDS = {
    # These were used in the prior audit script.  7/9 are wrong.
    'ENSG00000117724': {'claimed': 'PRC1',   'actual': 'CENPF',  'correct_ensg': 'ENSG00000198901'},
    'ENSG00000112118': {'claimed': 'MCM3',   'actual': 'MCM3',   'correct_ensg': 'ENSG00000112118'},  # ✅ correct
    'ENSG00000103126': {'claimed': 'DTYMK',  'actual': 'AXIN1',  'correct_ensg': 'ENSG00000168393'},
    'ENSG00000111206': {'claimed': 'CDKN3',  'actual': 'FOXM1',  'correct_ensg': 'ENSG00000100526'},
    'ENSG00000160200': {'claimed': 'LCAT',   'actual': 'CBS',    'correct_ensg': 'ENSG00000213398'},
    'ENSG00000125144': {'claimed': 'MT1F',   'actual': 'MT1G',   'correct_ensg': 'ENSG00000198417'},
    'ENSG00000145824': {'claimed': 'CXCL14', 'actual': 'CXCL14', 'correct_ensg': 'ENSG00000145824'},  # ✅ correct
    'ENSG00000125148': {'claimed': 'MT1X',   'actual': 'MT2A',   'correct_ensg': 'ENSG00000187193'},
    'ENSG00000132693': {'claimed': 'CRHBP',  'actual': 'CRP',    'correct_ensg': 'ENSG00000145708'},
}


# ═══════════════════════════════════════════════════════════════════════════
# METHOD B: AUDITABLE Z-SCORE SCORING
# ═══════════════════════════════════════════════════════════════════════════

import numpy as np
import pandas as pd


def method_b_score(expr_df, normal_idx, up_genes=None, down_genes=None,
                   verbose=False, audit_file=None):
    """
    Compute Method B composite score (leak-free z-score normalization).

    Method:
      1. For each gene, compute mean and SD from NORMAL samples only
      2. Z-score ALL samples using those normal-derived parameters
      3. Composite = mean(UP z-scores) - mean(DOWN z-scores)

    Parameters
    ----------
    expr_df : pd.DataFrame
        Samples × Genes expression matrix (already log2-transformed)
    normal_idx : list
        Sample IDs of normal/control samples (used for z-score reference)
    up_genes : list, optional
        Genes expected UP in tumor. Defaults to signature UP_GENES.
    down_genes : list, optional
        Genes expected DOWN in tumor. Defaults to signature DOWN_GENES.
    verbose : bool
        If True, print step-by-step computation details to stdout
    audit_file : str, optional
        If provided, write full audit trail to this file path

    Returns
    -------
    pd.Series
        Composite score for each sample (higher = more tumor-like)
    """
    if up_genes is None:
        up_genes = UP_GENES
    if down_genes is None:
        down_genes = DOWN_GENES

    audit_lines = []

    def log(msg):
        if verbose:
            print(msg)
        audit_lines.append(msg)

    log("=" * 90)
    log("METHOD B SCORING — COMPUTATION AUDIT TRAIL")
    log("=" * 90)
    log(f"Total samples: {len(expr_df)}")
    log(f"Normal samples used for z-score reference: {len(normal_idx)}")
    log(f"UP genes requested: {up_genes}")
    log(f"DOWN genes requested: {down_genes}")

    available_up = [g for g in up_genes if g in expr_df.columns]
    available_down = [g for g in down_genes if g in expr_df.columns]
    missing = [g for g in up_genes + down_genes if g not in expr_df.columns]

    log(f"\nAvailable UP genes:   {available_up}")
    log(f"Available DOWN genes: {available_down}")
    if missing:
        log(f"⚠ MISSING genes:     {missing}")

    # ── Step 1: Compute normal-derived reference statistics ──────────
    log(f"\n{'─' * 90}")
    log("STEP 1: Normal-cohort reference statistics (mean, SD)")
    log(f"{'─' * 90}")
    log(f"{'Gene':<12} {'Direction':<6} {'Normal Mean':>12} {'Normal SD':>12} {'Normal N':>10} {'SD OK?':>8}")
    log(f"{'─' * 62}")

    gene_params = {}
    for gene in available_up + available_down:
        ref_vals = expr_df.loc[normal_idx, gene].astype(float)
        mu = ref_vals.mean()
        sd = ref_vals.std()
        direction = 'UP' if gene in available_up else 'DOWN'
        sd_ok = '✅' if sd > 0.1 else '⚠ LOW'

        gene_params[gene] = {'mu': mu, 'sd': sd, 'direction': direction}
        log(f"  {gene:<12} {direction:<6} {mu:>12.6f} {sd:>12.6f} {len(ref_vals):>10d} {sd_ok:>8}")

    # ── Step 2: Z-score all samples ─────────────────────────────────
    log(f"\n{'─' * 90}")
    log("STEP 2: Z-score each sample using normal reference")
    log(f"  Formula: z = (sample_value - normal_mean) / normal_sd")
    log(f"{'─' * 90}")

    up_z = pd.DataFrame(index=expr_df.index)
    down_z = pd.DataFrame(index=expr_df.index)

    for gene in available_up:
        p = gene_params[gene]
        if p['sd'] > 0:
            up_z[gene] = (expr_df[gene].astype(float) - p['mu']) / p['sd']
        else:
            up_z[gene] = 0.0
            log(f"  ⚠ {gene}: SD=0, z-scores set to 0")

    for gene in available_down:
        p = gene_params[gene]
        if p['sd'] > 0:
            down_z[gene] = (expr_df[gene].astype(float) - p['mu']) / p['sd']
        else:
            down_z[gene] = 0.0
            log(f"  ⚠ {gene}: SD=0, z-scores set to 0")

    # Show a few sample z-scores for audit
    show_samples = list(normal_idx[:2]) + [s for s in expr_df.index if s not in normal_idx][:2]
    if show_samples:
        log(f"\n  Sample z-scores (first 2 normal + first 2 tumor):")
        header = f"  {'Sample':<25}"
        for gene in available_up + available_down:
            header += f" {gene:>10}"
        log(header)

        for sample in show_samples:
            row = f"  {str(sample)[:25]:<25}"
            for gene in available_up:
                if gene in up_z.columns:
                    row += f" {up_z.loc[sample, gene]:>+10.4f}"
            for gene in available_down:
                if gene in down_z.columns:
                    row += f" {down_z.loc[sample, gene]:>+10.4f}"
            log(row)

    # ── Step 3: Composite score ─────────────────────────────────────
    log(f"\n{'─' * 90}")
    log("STEP 3: Composite score = mean(UP z-scores) - mean(DOWN z-scores)")
    log(f"{'─' * 90}")

    up_mean = up_z.mean(axis=1) if up_z.shape[1] > 0 else pd.Series(0, index=expr_df.index)
    down_mean = down_z.mean(axis=1) if down_z.shape[1] > 0 else pd.Series(0, index=expr_df.index)
    composite = up_mean - down_mean

    if show_samples:
        log(f"\n  {'Sample':<25} {'UP_mean':>10} {'DOWN_mean':>10} {'Composite':>10}")
        for sample in show_samples:
            log(f"  {str(sample)[:25]:<25} {up_mean.loc[sample]:>+10.4f} {down_mean.loc[sample]:>+10.4f} {composite.loc[sample]:>+10.4f}")

    log(f"\n  Score distribution:")
    log(f"    Normal samples: mean={composite.loc[normal_idx].mean():+.4f}, "
        f"SD={composite.loc[normal_idx].std():.4f}, "
        f"range=[{composite.loc[normal_idx].min():+.4f}, {composite.loc[normal_idx].max():+.4f}]")

    non_normal = [s for s in expr_df.index if s not in normal_idx]
    if non_normal:
        log(f"    Tumor samples:  mean={composite.loc[non_normal].mean():+.4f}, "
            f"SD={composite.loc[non_normal].std():.4f}, "
            f"range=[{composite.loc[non_normal].min():+.4f}, {composite.loc[non_normal].max():+.4f}]")

    log("=" * 90)

    # Write audit file if requested
    if audit_file:
        with open(audit_file, 'w') as f:
            f.write('\n'.join(audit_lines))
        log(f"\nAudit trail written to: {audit_file}")

    return composite


# ═══════════════════════════════════════════════════════════════════════════
# SELF-TEST: Verify this module loads correctly
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("16-Gene MASLD-HCC Signature Reference")
    print("=" * 60)
    print(f"\nTotal genes: {len(ALL_GENES)}")
    print(f"  UP in HCC:   {len(UP_GENES)} → {UP_GENES}")
    print(f"  DOWN in HCC: {len(DOWN_GENES)} → {DOWN_GENES}")

    print(f"\nTier A+ ({len(TIER_AP)} genes): {TIER_AP}")
    print(f"Tier A  ({len(TIER_A)} genes):  {TIER_A}")
    print(f"Tier B  ({len(TIER_B)} genes):  {TIER_B}")

    print(f"\nEnsembl ID mapping:")
    for gene, info in SIGNATURE_GENES.items():
        print(f"  {gene:<10} → {info['ensembl_id']:<20} {info['direction']:<5} {info['grch38_coords']}")

    print(f"\nAffy HGU133A probe mapping:")
    for probe, symbol in AFFY_TO_SYMBOL.items():
        print(f"  {probe:<16} → {symbol}")

    # Integrity checks
    assert len(ALL_GENES) == 16, f"Expected 16 genes, got {len(ALL_GENES)}"
    assert len(UP_GENES) == 5, f"Expected 5 UP genes, got {len(UP_GENES)}"
    assert len(DOWN_GENES) == 11, f"Expected 11 DOWN genes, got {len(DOWN_GENES)}"
    assert len(set(ENSEMBL_TO_SYMBOL.values())) == 16, "Duplicate Ensembl IDs!"
    assert len(set(SYMBOL_TO_ENSEMBL.values())) == 16, "Duplicate symbols!"
    print("\n✅ All integrity checks passed.")
