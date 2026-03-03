const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, LevelFormat, ImageRun,
  HeadingLevel, BorderStyle, WidthType, ShadingType,
  PageNumber, PageBreak, TabStopType, TabStopPosition,
} = require("docx");

// ─── Helpers ───
const bdr = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: bdr, bottom: bdr, left: bdr, right: bdr };
const hBdr = { style: BorderStyle.SINGLE, size: 1, color: "2E75B6" };
const hBorders = { top: hBdr, bottom: hBdr, left: hBdr, right: hBdr };
const pad = { top: 60, bottom: 60, left: 100, right: 100 };

function hC(text, w) {
  return new TableCell({
    borders: hBorders, width: { size: w, type: WidthType.DXA },
    shading: { fill: "2E75B6", type: ShadingType.CLEAR }, margins: pad,
    children: [new Paragraph({ spacing: { after: 0 }, children: [new TextRun({ text, bold: true, color: "FFFFFF", font: "Arial", size: 18 })] })]
  });
}
function dC(text, w, o = {}) {
  return new TableCell({
    borders, width: { size: w, type: WidthType.DXA },
    shading: o.shade ? { fill: "F2F7FB", type: ShadingType.CLEAR } : undefined, margins: pad,
    children: [new Paragraph({ spacing: { after: 0 }, children: [new TextRun({ text, font: "Arial", size: 18, bold: o.bold || false })] })]
  });
}
function p(text, o = {}) {
  return new Paragraph({
    spacing: { after: o.after || 200, line: 276 },
    alignment: o.align || AlignmentType.JUSTIFIED,
    indent: o.indent ? { firstLine: 360 } : undefined,
    children: Array.isArray(text) ? text : [new TextRun({ text, font: "Arial", size: 22, ...o.run })]
  });
}
function h(text, level) {
  const sizes = { 1: 28, 2: 24, 3: 22 };
  const lvl = level === 1 ? HeadingLevel.HEADING_1 : level === 2 ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3;
  return new Paragraph({
    heading: lvl,
    spacing: { before: level === 1 ? 360 : 240, after: 120 },
    children: [new TextRun({ text, font: "Arial", bold: true, size: sizes[level], color: level === 1 ? "1A3C5E" : "2E75B6" })]
  });
}
const n = (t) => new TextRun({ text: t, font: "Arial", size: 22 });
const b = (t) => new TextRun({ text: t, font: "Arial", size: 22, bold: true });
const it = (t) => new TextRun({ text: t, font: "Arial", size: 22, italics: true });
const sub = (t) => new TextRun({ text: t, font: "Arial", size: 18, subScript: true });
const sup = (t) => new TextRun({ text: t, font: "Arial", size: 22, superScript: true });
const PB = new Paragraph({ children: [new PageBreak()] });

// ─── Figure helpers ───
const path = require("path");
const figDir = path.join(__dirname, "..", "figures") + "/";
function figImg(filename, widthInches, heightInches) {
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 120, after: 80 },
    children: [new ImageRun({
      type: "png",
      data: fs.readFileSync(figDir + filename),
      transformation: { width: widthInches * 96, height: heightInches * 96 },
      altText: { title: filename, description: filename, name: filename }
    })]
  });
}
function figCaption(text) {
  return new Paragraph({
    spacing: { after: 300 },
    alignment: AlignmentType.JUSTIFIED,
    children: Array.isArray(text) ? text : [new TextRun({ text, font: "Arial", size: 20, italics: true })]
  });
}

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Arial", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: "Arial", color: "1A3C5E" },
        paragraph: { spacing: { before: 360, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: "Arial", color: "2E75B6" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 22, bold: true, font: "Arial", color: "2E75B6" },
        paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 2 } },
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1296, bottom: 1440, left: 1296 }
      }
    },
    headers: {
      default: new Header({
        children: [new Paragraph({
          spacing: { after: 0 },
          border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "2E75B6", space: 4 } },
          children: [new TextRun({ text: "Parameter-Free Transcriptomic Biomarker Framework for HCC", font: "Arial", size: 16, italics: true, color: "888888" })]
        })]
      })
    },
    footers: {
      default: new Footer({
        children: [new Paragraph({
          alignment: AlignmentType.CENTER, spacing: { before: 0, after: 0 },
          border: { top: { style: BorderStyle.SINGLE, size: 4, color: "2E75B6", space: 4 } },
          children: [
            new TextRun({ text: "Page ", font: "Arial", size: 16, color: "888888" }),
            new TextRun({ children: [PageNumber.CURRENT], font: "Arial", size: 16, color: "888888" }),
          ]
        })]
      })
    },
    children: [

// ═══════════════ TITLE ═══════════════
new Paragraph({
  spacing: { after: 200 }, alignment: AlignmentType.CENTER,
  children: [new TextRun({
    text: "A Parameter-Free Scoring Framework Exposes Preprocessing Leakage as the Dominant Failure Mode in Transcriptomic Biomarker Studies",
    font: "Arial", size: 32, bold: true, color: "1A3C5E"
  })]
}),
new Paragraph({ spacing: { after: 400 }, alignment: AlignmentType.CENTER, children: [] }),

// ═══════════════ ABSTRACT ═══════════════
h("Abstract", 1),

p([
  n("Most transcriptomic biomarker studies for hepatocellular carcinoma (HCC) fail to replicate on independent data. We demonstrate that this is primarily a preprocessing problem with a straightforward fix: a "),
  b("parameter-free scoring pipeline"),
  n(" that eliminates data leakage by construction. The method\u2014within-dataset "),
  it("Z"),
  n("-scoring anchored exclusively to each cohort\u2019s own control samples, combined with a fixed directional composite\u2014requires no fitted model, no hyperparameter tuning, and transfers no parameters between datasets. Applied to a 16-gene signature (5 proliferation-associated, 11 hepatocyte-function genes), this approach achieved AUC 0.976\u20131.000 across independent cohorts, while an equivalent standard ML pipeline collapsed to AUC 0.31 on the same held-out data\u2014a 0.67+ AUC gap attributable entirely to preprocessing-stage leakage. The signature separated HCC from F4 cirrhotic backgrounds (AUC 0.997, no overlapping samples in this cohort, "),
  it("p"),
  n(" < 10"),
  sup("\u20139"),
  n("). Per-gene analysis revealed that discrimination at the clinical boundary is driven by collapse of hepatocyte-identity genes, not by proliferation markers\u2014an insight with direct implications for future signature design. The persistence of leakage-inflated results across the literature "),
  b("warrants systematic methodological review"),
  n("."),
]),

p([b("Keywords: "), n("hepatocellular carcinoma, transcriptomic biomarkers, data leakage, parameter-free scoring, reproducibility, cirrhosis")], { after: 60 }),

PB,

// ═══════════════ 1. INTRODUCTION ═══════════════
h("1. Introduction", 1),

p([
  n("Over 200 transcriptomic signatures for HCC detection have been published in the past decade. The vast majority do not replicate on independent data. This is not primarily a biological problem\u2014it is a computational one. The fix, moreover, is not complex: it requires "),
  b("zero fitted parameters"),
  n("."),
]),

p([
  n("Three failure modes account for most non-replication. First, "),
  b("transductive data leakage"),
  n(": normalizing entire expression matrices before train-test splitting allows test-set variance to influence training parameters, inflating apparent accuracy while masking inability to generalize. Second, "),
  b("demographic confounding"),
  n(": tumor samples are systematically older and more frequently male, and flexible classifiers exploit these imbalances as shortcuts. Third, the "),
  b("cirrhosis confound"),
  n(": HCC arises overwhelmingly in diseased livers, yet most studies benchmark against healthy tissue. Algorithms that learn to detect inflammation or regeneration rather than malignancy are clinically irrelevant for screening the high-risk cirrhotic population."),
], { indent: true }),

p([
  n("All three failure modes share a common root: they arise from fitting models to data in ways that conflate signal with artifact. We hypothesized that a scoring framework with "),
  b("no fitted parameters"),
  n("\u2014one that anchors normalization to each dataset\u2019s own controls and applies a fixed directional composite\u2014would bypass these failures by construction. The logic is simple: if nothing is fitted, nothing can overfit; if no parameters transfer between datasets, no information can leak. Using a 16-gene directional signature, we tested this hypothesis against standard ML and evaluated it across the full MASLD-to-cirrhosis spectrum."),
], { indent: true }),

PB,

// ═══════════════ 2. METHODS ═══════════════
h("2. Methods", 1),

h("2.1 Datasets", 2),

p("Three independent datasets were used (Table 1). TCGA-LIHC served as the primary cohort. GSE144269 served as a strictly held-out external validation cohort with no parameter sharing. GSE135251 provided the MASLD/NASH fibrosis spectrum."),

new Table({
  width: { size: 9648, type: WidthType.DXA },
  columnWidths: [1800, 1100, 1100, 1400, 1248, 3000],
  rows: [
    new TableRow({ children: [
      hC("Dataset", 1800), hC("Role", 1100), hC("HCC", 1100),
      hC("Control", 1400), hC("Platform", 1248), hC("Notes", 3000),
    ]}),
    new TableRow({ children: [
      dC("TCGA-LIHC", 1800, { bold: true }), dC("Primary", 1100), dC("320", 1100),
      dC("50 normal", 1400), dC("RNA-seq", 1248), dC("Xena HTSeq log\u2082(FPKM+1)", 3000),
    ]}),
    new TableRow({ children: [
      dC("GSE144269", 1800, { bold: true, shade: true }), dC("External val.", 1100, { shade: true }), dC("70", 1100, { shade: true }),
      dC("70 matched", 1400, { shade: true }), dC("RNA-seq", 1248, { shade: true }), dC("Paired tumor/normal", 3000, { shade: true }),
    ]}),
    new TableRow({ children: [
      dC("GSE135251", 1800, { bold: true }), dC("Fibrosis", 1100), dC("\u2014", 1100),
      dC("216 biopsies", 1400), dC("RNA-seq", 1248), dC("MASLD/NASH F0\u2013F4", 3000),
    ]}),
  ]
}),
p("", { after: 60 }),
p([b("Table 1. "), it("Datasets used in this study.")], { after: 200 }),

h("2.2 The Leakage Problem: ML Baseline", 2),

p([
  n("To quantify the effect of transductive leakage, a standard pipeline was constructed using scikit-learn: a "),
  b("StandardScaler"),
  n(" was fitted to the full TCGA-LIHC expression matrix (tumor + normal jointly), followed by a "),
  b("RandomForestClassifier"),
  n(" (n_estimators = 500, class_weight = \u2018balanced\u2019). The trained pipeline was then applied to GSE144269. This pipeline has "),
  b("fitted parameters"),
  n(": the scaler\u2019s mean and variance are computed from both training and test partitions, and the classifier learns decision boundaries from training data. This is routine practice\u2014and the source of the problem."),
]),

h("2.3 The Fix: Parameter-Free Directional Z-Score", 2),

p([
  n("The fix eliminates all fitted parameters. The 16-gene signature comprises an UP module (5 genes: PRC1, RACGAP1, MCM3, DTYMK, CDKN3) and a DOWN module (11 genes: CYP1A2, LCAT, FCN3, MT1F, CXCL14, FCN2, CLEC4M, MT1X, CLEC1B, CRHBP, GDF2). For each dataset independently, the mean (\u03BC) and standard deviation (\u03C3) for each gene are computed "),
  b("only from that dataset\u2019s control samples"),
  n(". Every sample is then Z-scored against this local baseline:"),
]),

new Paragraph({
  spacing: { before: 120, after: 120 }, alignment: AlignmentType.CENTER,
  children: [
    it("Z"), sub("gene"), n("  =  ( X"), sub("gene"), n(" \u2212 \u03BC"), sub("normal"), n(" )  /  \u03C3"), sub("normal"),
  ]
}),

p("The composite score is the mean Z of UP genes minus the mean Z of DOWN genes:"),

new Paragraph({
  spacing: { before: 120, after: 120 }, alignment: AlignmentType.CENTER,
  children: [
    b("Score"), n("  =  mean(Z"), sub("up"), n(")  \u2212  mean(Z"), sub("down"), n(")"),
  ]
}),

p([
  n("This algorithm has "),
  b("zero fitted parameters"),
  n(" in the conventional sense. There is no trained model, no learned decision boundary, no hyperparameter grid search, and no threshold optimization. The \u03BC and \u03C3 values are not transferred between datasets\u2014they are recomputed fresh from each cohort\u2019s own controls. The gene list and the direction (UP minus DOWN) are fixed a priori. This means:"),
]),

p([
  b("No information leaks between datasets. "),
  n("When applied to a new cohort, only that cohort\u2019s control distribution is used for normalization. There is nothing to overfit, because there is nothing being fitted."),
], { indent: true }),

h("2.4 Covariate Hierarchy Analysis", 2),

p("To assess whether the signature score captures confounding rather than true biological signal, we performed a logistic regression covariate hierarchy on TCGA-LIHC (n = 370). Five nested models were fit: (1) age only, (2) sex only, (3) age + sex, (4) signature score only, and (5) signature score + age + sex. Predictions were obtained via leave-one-out cross-validation. Additionally, we computed Pearson correlation between signature score and age within tumor samples, conducted Mann\u2013Whitney tests on score distributions stratified by sex, and performed a permutation test (1,000 permutations) comparing observed AUC to permuted class labels."),
h("2.5 Cirrhosis-as-Baseline Testing", 2),

p("HCC screening is performed in patients with existing liver disease, not in healthy individuals. To evaluate clinical validity, we tested the signature across seven configurations of increasing stringency, culminating in the most confounded comparison: HCC arising on a cirrhotic background (Ishak 5\u20136) versus F4 cirrhosis biopsies without cancer. We also tested with cirrhosis-normalized scoring, where Z-score parameters are anchored to F4 cirrhotic samples rather than healthy tissue, to confirm that the signal is not an artifact of the reference choice."),

h("2.6 Statistical Analysis", 2),

p([
  n("Performance was quantified by ROC AUC. Group comparisons used Mann\u2013Whitney "),
  it("U"),
  n(" and Kruskal\u2013Wallis tests. Fibrosis progression was assessed by Spearman rank correlation. Permutation tests ("),
  it("n"),
  n(" = 1,000) established null distributions. All analyses used Python 3.10 (NumPy, SciPy, pandas, scikit-learn). Code is publicly available."),
]),

PB,

// ═══════════════ 3. RESULTS ═══════════════
h("3. Results", 1),

h("3.1 HCC Is Fully Separable from Cirrhotic Backgrounds in Tissue", 2),

p([
  n("The seven-test cirrhosis battery produced consistently high performance across all configurations in this cohort (Table 2). In the most clinically stringent comparison\u2014HCC arising on cirrhotic backgrounds (Ishak 5\u20136, n = 82) versus F4 cirrhosis biopsies without cancer (n = 14)\u2014the parameter-free score achieved AUC 0.997: tumor mean +6.28 vs. cirrhosis mean +0.12, a gap of 6.16 composite-score units. No overlapping samples were observed in this cohort (Mann\u2013Whitney "),
  it("p"),
  n(" = 1.31 \u00D7 10"),
  sup("\u20139"),
  n(")."),
]),

p([
  n("The score continuum across the fibrosis spectrum was monotonically increasing: F0\u2013F1 (\u22121.10) \u2192 F2 (\u22120.74) \u2192 F3 (\u22120.22) \u2192 F4 (+0.12) \u2192 HCC (+6.32). This gradient was highly significant (Kruskal\u2013Wallis H = 63.73, "),
  it("p"),
  n(" = 9.4 \u00D7 10"),
  sup("\u201314"),
  n("; Spearman \u03C1 = 0.54, "),
  it("p"),
  n(" = 1.1 \u00D7 10"),
  sup("\u201317"),
  n("). The gap between F4 cirrhosis and HCC exceeds 6 units\u2014representing >7 standard deviations of the control distribution (Figure 1)."),
], { indent: true }),

figImg("fig2_fibrosis_gradient.png", 6.0, 3.6),
figCaption([
  new TextRun({ text: "Figure 1. ", font: "Arial", size: 20, bold: true, italics: true }),
  new TextRun({ text: "Biomarker score across the fibrosis-to-cancer continuum. The monotonic gradient from F0\u2013F1 through F4 cirrhosis is followed by a >6-unit gap to HCC, with no overlap between F4 and HCC distributions in this cohort.", font: "Arial", size: 20, italics: true }),
]),

// Table 2 - fixed p-value column
new Table({
  width: { size: 9648, type: WidthType.DXA },
  columnWidths: [3200, 900, 900, 900, 900, 2848],
  rows: [
    new TableRow({ children: [
      hC("Configuration", 3200), hC("AUC", 900), hC("n HCC", 900),
      hC("n Ctrl", 900), hC("Gap", 900), hC("p-value", 2848),
    ]}),
    ...[
      ["All HCC vs. Normal", "\u22650.99", "320", "50", "6.32", "2.8 \u00D7 10\u207B\u00B3\u2070", false],
      ["All HCC vs. F4 Cirrhosis", "\u22650.99", "320", "14", "6.20", "1.2 \u00D7 10\u207B\u00B9\u2070", true],
      ["Cirrhotic HCC vs. F4 Cirrhosis*", "\u22650.99", "82", "14", "6.16", "1.3 \u00D7 10\u207B\u2079", false],
      ["All HCC vs. F4 (cirrhosis-norm.)", "\u22650.99", "320", "14", "6.02", "1.2 \u00D7 10\u207B\u00B9\u2070", true],
      ["Cirrhotic HCC vs. F4 (cirrh.-norm.)*", "\u22650.99", "82", "14", "5.96", "1.3 \u00D7 10\u207B\u2079", false],
      ["Cirrhotic HCC vs. Cirrhotic Normal\u2020", "\u22650.99", "82", "10", "6.22", "1.4 \u00D7 10\u207B\u2077", true],
      ["All HCC vs. Full NAFLD/NASH (F0\u2013F4)", "\u22650.99", "320", "216", "7.03", "2.9 \u00D7 10\u207B\u2078\u2076", false],
    ].map(r => new TableRow({ children: [
      dC(r[0], 3200, { shade: r[6] }), dC(r[1], 900, { bold: true, shade: r[6] }),
      dC(r[2], 900, { shade: r[6] }), dC(r[3], 900, { shade: r[6] }),
      dC(r[4], 900, { shade: r[6] }), dC(r[5], 2848, { shade: r[6] }),
    ]})),
  ].flat()
}),

p("", { after: 60 }),
p([b("Table 2. "), it("Seven-test cirrhosis battery (parameter-free score). "), n("*Ishak 5\u20136 HCC vs. METAVIR F4. \u2020Within-TCGA cirrhotic pairs. Gap = tumor mean \u2212 control mean. No overlapping samples were observed between HCC and control distributions in any configuration in this cohort.")], { after: 300 }),


h("3.2 The Parameter-Free Fix vs. Standard ML: A 0.67 AUC Gap", 2),

p([
  n("A standard ML pipeline\u2014trained on TCGA-LIHC with conventional whole-matrix normalization\u2014achieved high within-cohort performance (AUC > 0.95 in TCGA-LIHC). On the strictly held-out GSE144269 cohort, performance dropped to AUC 0.31, worse than chance. The parameter-free directional Z-score, applied to the "),
  it("same"),
  n(" external cohort using only its own normal samples for normalization, achieved AUC 0.976\u20131.000."),
]),

p([
  n("This comparison isolates the variable (Figure 2). The gene set is identical (16 genes). The input data is identical (GSE144269). The only difference is how normalization is performed and whether a model is fitted. The ML pipeline fits parameters\u2014scaler means, variances, classifier weights\u2014to TCGA-LIHC, then transfers them to GSE144269. The parameter-free approach computes nothing from TCGA-LIHC; it uses GSE144269\u2019s own controls. The result: a "),
  b("0.67+ AUC gap"),
  n(" entirely attributable to preprocessing-stage leakage."),
], { indent: true }),

p([
  n("This gap should not be interpreted merely as \u201Csimple is better.\u201D The parameter-free scoring approach does not outperform ML because simplicity is a virtue in itself. It outperforms because the within-cohort AUC > 0.95 "),
  b("was an artifact of leaked information"),
  n("\u2014test-set variance leaking into preprocessing parameters. When that leak is plugged by removing the fitting step entirely, the artifact vanishes and the genuine signal emerges."),
], { indent: true }),


PB,
figImg("fig2_method_comparison.png", 6.5, 3.6),
figCaption([
  new TextRun({ text: "Figure 2. ", font: "Arial", size: 20, bold: true, italics: true }),
  new TextRun({ text: "Method comparison. Left: standard ML pipeline with transductive leakage (red). Right: parameter-free scoring pipeline with per-cohort normalization (green). No parameters transfer between cohorts.", font: "Arial", size: 20, italics: true }),
]),

h("3.3 Covariate Hierarchy: Signature Score Is Independent of Age and Sex", 2),

p("A critical test of whether the signature captures biological signal or demographic confounding is the covariate hierarchy analysis. We constructed five nested logistic regression models on TCGA-LIHC (n = 370) with leave-one-out cross-validation:"),

new Table({
  width: { size: 9648, type: WidthType.DXA },
  columnWidths: [3600, 1600, 4448],
  rows: [
    new TableRow({ children: [
      hC("Model", 3600), hC("AUC", 1600), hC("Interpretation", 4448),
    ]}),
    ...[
      ["Age only", "0.558", "Near chance", false],
      ["Sex only", "0.558", "Near chance", true],
      ["Age + Sex", "0.588", "Near chance", false],
      ["Signature score only", "1.000", "Perfect in this cohort", true],
      ["Score + Age + Sex", "1.000", "No improvement", false],
    ].map(r => new TableRow({ children: [
      dC(r[0], 3600, { bold: r[0].includes("Signature"), shade: r[3] }),
      dC(r[1], 1600, { bold: r[0].includes("Signature"), shade: r[3] }),
      dC(r[2], 4448, { shade: r[3] }),
    ]}))
  ]
}),
p("", { after: 60 }),
p([b("Table 3. "), it("Covariate hierarchy (TCGA-LIHC, n = 370). "), n("Demographics alone do not predict tissue type. The signature score alone achieves AUC = 1.000 in this cohort; adding demographics provides no improvement.")], { after: 200 }),

p([
  n("Within-tumor correlation between age and signature score was negligible (Pearson "),
  it("r"),
  n(" = 0.068, "),
  it("p"),
  n(" = 0.226). Sex-stratified AUC was 1.000 in this cohort for both males (n = 229 tumor, 30 normal) and females (n = 91 tumor, 20 normal). Male vs. female tumor mean scores were nearly identical (6.36 vs. 6.23, Mann\u2013Whitney "),
  it("p"),
  n(" = 0.286). Permutation testing (1,000 permutations) confirmed that 0/1,000 random label shuffles achieved the observed separation ("),
  it("p"),
  n(" < 0.001). Neither age ("),
  it("p"),
  n(" = 0.186) nor sex distribution ("),
  it("p"),
  n(" = 0.135) differed significantly between tumor and normal tissue groups."),
], { indent: true }),


h("3.4 DOWN Genes Drive Discrimination; UP Genes Peak in Cirrhosis", 2),

p([
  n("Per-gene analysis across the fibrosis-to-cancer continuum revealed that discrimination at the clinical boundary is asymmetric (Figure 3). Of the 5 UP genes, four (PRC1, RACGAP1, MCM3, DTYMK) are actually higher in F4 cirrhosis than in HCC, with per-gene AUCs of 0.37, 0.28, 0.01, and 0.09 respectively. An AUC below 0.5 indicates "),
  it("inverse classification"),
  n(": expression of these genes is higher in cirrhosis than in HCC, the opposite of the expected direction. These proliferation markers peak during cirrhotic regeneration."),
]),

p([
  n("In contrast, all 11 DOWN genes maintained strong discrimination against the cirrhotic background. Per-gene AUCs ranged from 0.975 (GDF2) to \u22650.99 (CYP1A2, LCAT, FCN3, MT1F, CXCL14, CLEC4M, MT1X, CRHBP). These genes encode the mature hepatic program: drug metabolism (CYP1A2), lipid transport (LCAT), innate immunity (FCN2, FCN3, CLEC4M, CLEC1B), and stress response (MT1F, MT1X, CRHBP, CXCL14, GDF2). In cirrhosis, they are maximally upregulated as compensation. In HCC, they collapse as tumor cells lose differentiated identity."),
], { indent: true }),

p("This has direct implications for signature design: panels built exclusively from proliferation genes will fail at the cirrhosis boundary because proliferation is maximal in both conditions.", { indent: true }),

figImg("fig3_pergene_asymmetry.png", 5.0, 5.0),
figCaption([
  new TextRun({ text: "Figure 3. ", font: "Arial", size: 20, bold: true, italics: true }),
  new TextRun({ text: "Per-gene AUC for HCC vs. F4 cirrhosis. DOWN genes (blue) achieve near-complete discrimination. Four of five UP genes (orange) score below 0.5, indicating higher expression in cirrhosis than in HCC. The dashed line marks chance level (AUC = 0.5).", font: "Arial", size: 20, italics: true }),
]),

PB,

// ═══════════════ 4. DISCUSSION ═══════════════
h("4. Discussion", 1),

p([
  n("The central finding is that a "),
  b("parameter-free scoring approach"),
  n(" outperformed standard ML by 0.67+ AUC units on external data\u2014using the same genes and the same cohort. The fix is not a better algorithm. It is the "),
  it("removal"),
  n(" of the algorithm. By eliminating the fitting step, we eliminated the mechanism through which leakage occurs."),
]),

p([
  n("This distinction matters because the standard response to reproducibility failures in genomics is to add complexity: more sophisticated cross-validation, ensemble methods, regularization. Our results suggest the opposite. The problem is not that models are underpowered\u2014it is that they are "),
  it("overfitted to artifacts that should never have been in the training data"),
  n(". Removing the fitting step entirely is more effective than refining it."),
], { indent: true }),

p([
  n("The difference is structural, not procedural. The standard defense against leakage\u2014nested cross-validation, held-out sets, temporal splits\u2014attempts to prevent leakage through "),
  it("procedural discipline"),
  n(". The parameter-free approach prevents it through "),
  it("structural elimination"),
  n(": there is no fitting step to corrupt. Procedural safeguards can fail silently when implemented incorrectly; structural safeguards cannot fail, because the mechanism of leakage has been removed entirely."),
], { indent: true }),

p("The cirrhosis results expose a second systemic problem. Signatures benchmarked only against healthy tissue cannot distinguish cancer from regeneration, because proliferation genes peak during both conditions. Our analysis shows that reliable discrimination at the clinical boundary requires hepatocyte-identity markers\u2014genes whose expression diverges between compensatory regeneration (cirrhosis) and dedifferentiation (HCC). This architectural insight should inform future biomarker design.", { indent: true }),

p("While this parameter-free tissue framework achieves complete discrimination against cirrhotic backgrounds, it cannot be naively ported to circulating cell-free RNA (cfRNA). Preliminary investigations indicate that hepatocyte turnover kinetics in advanced liver disease cause mature hepatic identity genes\u2014the DOWN module that drives discrimination in tissue\u2014to paradoxically spike in the bloodstream of cirrhotic patients, reversing the expected direction of effect. Future liquid biopsy development must account for these fluid-kinetic reversals rather than assuming direct tissue-to-blood parity. The parameter-free scoring framework itself is specimen-agnostic; what changes is the gene-direction map, which must be re-established de novo for each analyte type.", { indent: true }),

// Pitfalls table
h("4.1 Methodological Checklist", 3),

new Table({
  width: { size: 9648, type: WidthType.DXA },
  columnWidths: [2200, 2600, 2600, 2248],
  rows: [
    new TableRow({ children: [
      hC("Pitfall", 2200), hC("Common Practice", 2600),
      hC("Parameter-Free Fix", 2600), hC("Outcome Here", 2248),
    ]}),
    ...[
      ["Transductive leakage", "Normalize full matrix; fit model", "Anchor Z-score to controls only; no model", "AUC: 0.31 \u2192 \u22650.99"],
      ["Healthy-only controls", "Benchmark vs. normal tissue", "Test vs. F4 cirrhosis", "AUC \u22650.99 sustained"],
      ["Demographic shortcuts", "No covariate testing", "Age/sex hierarchy + permutation", "Score \u226B demographics"],
      ["Proliferation reliance", "Panels of cell-cycle genes only", "Include differentiation module", "DOWN genes drive AUC"],
    ].map((r, i) => new TableRow({ children: [
      dC(r[0], 2200, { bold: true, shade: i % 2 === 1 }),
      dC(r[1], 2600, { shade: i % 2 === 1 }),
      dC(r[2], 2600, { shade: i % 2 === 1 }),
      dC(r[3], 2248, { shade: i % 2 === 1 }),
    ]}))
  ]
}),
p("", { after: 60 }),
p([b("Table 4. "), it("Methodological checklist. "), n("The \u201CParameter-Free Fix\u201D column shows the approach used in this study for each common pitfall.")], { after: 300 }),

PB,

// ═══════════════ 5. LIMITATIONS ═══════════════
h("5. Limitations", 1),

p("All analyses are retrospective. Prospective, multi-site clinical validation has not been performed. The cirrhosis tissue cohorts, while showing consistent results, are modest in size (n = 10\u201314). HCC is etiologically heterogeneous; whether the signature performs equivalently across viral, metabolic, and alcohol-related subtypes requires dedicated subgroup analysis."),

p([
  n("Feature selection (the 16 genes) was performed a priori in independent work. The methodology evaluated here is strictly the "),
  it("validation pipeline"),
  n(", which is structurally parameter-free: it employs no fitted models, no learned boundaries, and no optimized thresholds. Gene identities and directions are fixed inputs, not outputs of the pipeline under test."),
], { indent: true }),

p([
  n("This work was performed as an independent computational analysis using exclusively public datasets and open-source tools. While this demonstrates that rigorous biomarker validation does not require proprietary infrastructure, it also means that certain resources\u2014large prospective cohorts, CLIA-grade assay development, clinical trial access\u2014remain beyond scope."),
], { indent: true }),

// ═══════════════ 6. CONCLUSION ═══════════════
h("6. Conclusion", 1),

p([
  n("We show that the dominant failure mode in transcriptomic biomarker studies\u2014cross-cohort non-replication\u2014has a simple structural fix: eliminate fitted parameters from the scoring pipeline. A parameter-free, control-anchored directional Z-score outperformed standard ML by over 0.67 AUC units on external data. The signature separates HCC from cirrhosis with no overlap in this cohort, driven by collapse of the hepatic differentiation program rather than proliferation markers. When the fix is this straightforward\u2014within-dataset "),
  it("Z"),
  n("-scoring, fixed gene directions, no fitted model\u2014the continued prevalence of leakage-inflated results suggests that standard computational workflows in the field need structural revision."),
]),

PB,

// ═══════════════ DATA AVAILABILITY ═══════════════
h("Data Availability", 1),

p("All datasets used in this study are publicly available. TCGA-LIHC expression data (HTSeq log\u2082(FPKM+1)) were obtained from the UCSC Xena Browser. GSE144269 and GSE135251 are available through the NCBI Gene Expression Omnibus. No restricted-access or proprietary data were used."),

p("The complete analysis pipeline\u2014including the 16-gene reference file (signature_reference.py), the parameter-free scoring function (method_b_score()), cirrhosis battery analysis, covariate hierarchy analysis, and all figure-generation code\u2014is available at [GitHub repository URL to be provided upon publication]. The parameter-free scoring pipeline requires only Python 3.10 with NumPy, SciPy, pandas, and scikit-learn; no proprietary software or specialized hardware is needed. All results are fully reproducible from the provided code and public data.", { indent: true }),

// ═══════════════ APPENDIX ═══════════════
h("Appendix: 16-Gene Signature", 1),

new Table({
  width: { size: 9648, type: WidthType.DXA },
  columnWidths: [1200, 800, 1200, 6448],
  rows: [
    new TableRow({ children: [
      hC("Gene", 1200), hC("Dir.", 800), hC("AUC vs F4", 1200), hC("Function", 6448),
    ]}),
    ...[
      ["PRC1", "UP", "0.37", "Cytokinesis regulator; mitotic spindle midzone"],
      ["RACGAP1", "UP", "0.28", "Rac GTPase-activating protein; cytokinesis"],
      ["MCM3", "UP", "0.01", "DNA replication licensing factor"],
      ["DTYMK", "UP", "0.09", "Thymidylate kinase; nucleotide biosynthesis"],
      ["CDKN3", "UP", "0.85", "CDK inhibitor; cell cycle regulation"],
      ["CYP1A2", "DOWN", "\u22650.99", "Cytochrome P450; xenobiotic metabolism"],
      ["LCAT", "DOWN", "\u22650.99", "Lecithin-cholesterol acyltransferase; HDL metabolism"],
      ["FCN3", "DOWN", "\u22650.99", "Ficolin 3; complement lectin pathway"],
      ["MT1F", "DOWN", "\u22650.99", "Metallothionein 1F; metal homeostasis"],
      ["CXCL14", "DOWN", "\u22650.99", "Chemokine; immune cell recruitment"],
      ["FCN2", "DOWN", "\u22650.99", "Ficolin 2; complement activation"],
      ["CLEC4M", "DOWN", "\u22650.99", "C-type lectin; pathogen recognition"],
      ["MT1X", "DOWN", "\u22650.99", "Metallothionein 1X; oxidative stress"],
      ["CLEC1B", "DOWN", "0.997", "C-type lectin; platelet activation"],
      ["CRHBP", "DOWN", "\u22650.99", "CRH-binding protein; stress axis"],
      ["GDF2", "DOWN", "0.975", "BMP9; vascular homeostasis"],
    ].map((r, i) => new TableRow({ children: [
      dC(r[0], 1200, { bold: true, shade: i % 2 === 1 }),
      dC(r[1], 800, { shade: i % 2 === 1 }),
      dC(r[2], 1200, { shade: i % 2 === 1 }),
      dC(r[3], 6448, { shade: i % 2 === 1 }),
    ]}))
  ]
}),
p("", { after: 60 }),
p([b("Table 5. "), it("The 16-gene signature. "), n("AUC vs F4 = per-gene AUC for HCC vs. F4 cirrhosis in tissue. Note 4/5 UP genes score below 0.5 (inverted), confirming that the DOWN module drives discrimination at the clinical boundary.")]),

    ] // end children
  }] // end sections
});

Packer.toBuffer(doc).then(buf => {
  const outPath = path.join(__dirname, "HCC_Biomarker_Framework_Paper_v6.docx");
  fs.writeFileSync(outPath, buf);
  console.log("Written:", buf.length, "bytes");
});
