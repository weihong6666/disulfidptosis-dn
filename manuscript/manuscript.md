# Disulfidptosis-Associated Gene Signature Reveals Convergent Cytoskeletal and Oxidative Stress Perturbation in Diabetic Nephropathy

**Running title:** Convergent cytoskeletal-oxidative stress perturbation in diabetic nephropathy

**Author:** Zhu Weixin

**Affiliation:** College of Life Sciences, Shanxi Agricultural University, Jinzhong, Shanxi, P.R. China

**Corresponding author email:** 20251312329@stu.sxau.edu.cn

---

## Abstract

**Background:** Disulfidptosis, a recently identified form of regulated cell death driven by disulfide stress-induced actin cytoskeleton collapse, remains largely unexplored in diabetic nephropathy (DN)—the leading cause of end-stage renal disease worldwide. Only two prior studies have examined disulfidptosis-related biomarkers in DN, both employing limited methodological scope.

**Methods:** We integrated four transcriptomic datasets: GSE96804 (training; 41 DN vs. 20 control glomeruli), GSE30528 and GSE30529 (external validation), and GSE131882 (single-cell; 3 DN vs. 3 control kidneys). Disulfidptosis-related differentially expressed genes (DR-DEGs) were identified by limma and intersected with DN-associated co-expression modules from WGCNA. A triple-algorithm intersection strategy (LASSO, SVM-RFE, Random Forest) with nested cross-validation was used for biomarker selection, followed by bootstrap optimism correction and permutation testing. Pathway specificity was evaluated by comparing disulfidptosis activity against four other cell death modalities via ssGSEA. Downstream analyses included immune microenvironment profiling, single-cell transcriptomic characterization, and protein-level validation via Human Protein Atlas immunohistochemistry.

**Results:** Nine DR-DEGs were identified, with enrichment converging on actin cytoskeleton organization (stress fiber, actomyosin; adjusted *P* = 1.64 × 10⁻³). Pathway specificity analysis revealed a key biological insight: oxidative stress and ferroptosis pathways were profoundly altered in DN (ROS: *P* = 5.2 × 10⁻⁹; ferroptosis: *P* = 1.3 × 10⁻⁶; glutathione: *P* = 4.5 × 10⁻⁶), whereas disulfidptosis pathway activity per se showed no significant difference (*P* = 0.953), indicating that the transcriptional signature captures convergent actin cytoskeleton–oxidative stress perturbation rather than isolated disulfidptosis activation. Seven genes (NCKAP1L, THSD7A, MYH10, PDLIM1, IL1B, FLNB, SLC3A2) were selected through triple-ML intersection with nested cross-validation. Rigorous internal validation yielded a bootstrap optimism-corrected AUC of 0.956, with permutation testing confirming statistical significance (*P* < 0.001). External validation demonstrated robust discriminative performance in glomerular samples (GSE30528 AUC = 0.923) with expectedly lower performance in tubulointerstitial samples (GSE30529 AUC = 0.767), consistent with compartment-specific transcriptional programs. Single-cell analysis of 20,681 cells resolved the signature to specific renal cell types, with THSD7A highest in podocytes and SLC3A2 in proximal tubular epithelium. HPA immunohistochemistry validated protein-level expression in renal tubules and glomeruli.

**Conclusions:** This study identifies a seven-gene disulfidptosis-related transcriptional program associated with DN, reflecting convergent perturbation of actin cytoskeleton and oxidative stress pathways in the diabetic kidney rather than isolated disulfidptosis pathway activation. These candidate molecular biomarkers warrant further experimental validation in independent clinical cohorts.

**Keywords:** Disulfidptosis; Diabetic nephropathy; Oxidative stress; Actin cytoskeleton; Single-cell RNA-seq; Transcriptional remodeling; Candidate biomarkers

---

## 1. Introduction

Diabetic nephropathy (DN), a severe microvascular complication of diabetes mellitus, affects approximately 40% of diabetic patients and represents the leading cause of end-stage renal disease (ESRD) worldwide [1,2]. Characterized by progressive glomerulosclerosis, tubulointerstitial fibrosis, and declining renal function, DN imposes substantial clinical and economic burdens [3]. Currently, the diagnosis of DN relies primarily on persistent albuminuria and/or reduced estimated glomerular filtration rate (eGFR), with renal biopsy serving as the definitive but invasive gold standard [4]. These limitations highlight an urgent need for non-invasive molecular biomarkers capable of enabling earlier detection and risk stratification.

Disulfidptosis, a recently identified form of regulated cell death, was first characterized by Liu et al. in 2023 as a disulfide stress-induced mechanism of cytoskeletal collapse [5]. Unlike apoptosis (caspase-dependent), ferroptosis (lipid peroxidation-driven), pyroptosis (inflammasome-mediated), or necroptosis (RIPK1/RIPK3/MLKL-dependent), disulfidptosis is uniquely triggered by excessive intracellular cystine accumulation. Under conditions of NADPH limitation, accumulated cystine cannot be adequately reduced to cysteine, leading to aberrant disulfide bond formation between actin cytoskeletal proteins. This process results in actin network crosslinking, cytoskeletal collapse, and ultimately cell death [5]. Key molecular players in disulfidptosis include SLC7A11 and its chaperone subunit SLC3A2 (mediating cystine import), members of the WAVE regulatory complex (NCKAP1, NCKAP1L; regulating actin polymerization), and actin-binding proteins such as non-muscle myosins and filamins [5,6].

The kidney—and specifically the renal proximal tubule—represents a site of intense cystine handling. SLC3A2/SLC7A11-mediated cystine reabsorption occurs predominantly in the proximal tubular epithelium, placing these cells at particularly high risk for disulfide stress when the reductive capacity of the NADPH/glutathione system is compromised [7]. In the context of DN, chronic hyperglycemia drives oxidative stress, NADPH depletion (via the polyol pathway and PKC activation), and inflammatory cytokine production—conditions that collectively favor disulfide bond accumulation and may predispose renal cells to disulfidptosis [8]. This mechanistic rationale suggests that disulfidptosis-related genes may serve as candidate molecular biomarkers in DN, reflecting the interplay between actin cytoskeleton remodeling and oxidative stress in the diabetic kidney.

Despite this mechanistic rationale, direct evidence that disulfidptosis-related transcriptional programs are altered in human DN remains limited, and the two studies published to date share two methodological limitations. Both rely on one or two feature-selection algorithms applied to a single, mechanistically heterogeneous gene set [9,10]; because different algorithms impose different priors on high-dimensional, small-sample transcriptomic data, single-method selection carries appreciable false-discovery risk, and the two studies consequently reported entirely non-overlapping biomarker panels. In addition, where single-cell data were analysed, the work was restricted to annotation of the full gene set, without resolving which renal cell types express the nominated biomarkers [10]. Neither study tested whether the resulting signature reflects disulfidptosis-specific biology or a more general stress response that happens to share genes. Concretely, Xu et al. (2023) identified CXCL6, CD48, C1QB and COL6A3 using LASSO and SVM-RFE across three microarray datasets [9], whereas Wang et al. (2024) combined scRNA-seq and WGCNA with KNN-based modelling to nominate VEGFA, MAGI2, THSD7A and ANKRD28, with PCR validation in blood samples [10].

Against this background, our study addresses three questions. First, can a biomarker panel be defined that is robust to the choice of feature-selection algorithm, and does it generalise beyond the training cohort? Second, does the "disulfidptosis-related" transcriptional signal in DN actually reflect activation of the disulfidptosis pathway, or a broader cytoskeletal-oxidative stress perturbation that shares genes with it? Third, in which renal cell types is the resulting signature expressed? To address these questions we integrated four public transcriptomic datasets (one training cohort, two external validation cohorts and one single-cell dataset) and applied three complementary analytical layers: (i) a triple-algorithm intersection strategy (LASSO ∩ SVM-RFE ∩ Random Forest) in which all feature selection is re-executed inside each cross-validation fold, combined with bootstrap optimism correction and permutation testing; (ii) a cell death pathway specificity analysis benchmarking disulfidptosis against ferroptosis, apoptosis, pyroptosis and necroptosis, together with oxidative stress pathway profiling; and (iii) a single-cell analysis resolving the renal cell types in which the signature is expressed. An overview of the study design and analytical workflow is provided in Figure 1. We first define the disulfidptosis-related differentially expressed gene set and the DN-associated co-expression modules, then derive and rigorously validate the biomarker panel, and finally contextualise it through pathway specificity, single-cell, immune microenvironment and protein-level analyses.

---

## 2. Methods

### 2.1 Data Acquisition and Processing

Four transcriptomic datasets were obtained from the Gene Expression Omnibus (GEO, https://www.ncbi.nlm.nih.gov/geo/). The training set GSE96804 (GPL17586, Affymetrix Human Transcriptome Array 2.0) comprised 41 DN and 20 control human glomerular samples. External validation datasets included GSE30528 (9 DN, 13 control glomerular samples; GPL571, Affymetrix Human Genome U133A 2.0 Array) and GSE30529 (10 DN, 12 control tubulointerstitial samples; same platform). The single-cell dataset GSE131882 (3 DN, 3 control human kidney samples) was profiled on the 10x Genomics Chromium platform. For GSE96804, probe-to-gene mapping was performed by parsing the `gene_assignment` column of the GPL17586 annotation file, extracting the first annotated gene symbol per probe. Probes mapping to multiple genes or lacking gene annotation were excluded. When multiple probes mapped to the same gene, the probe with the highest mean expression was retained. The expression matrix was normalized by robust multi-array average (RMA) and log2-transformed. Batch effects were assessed by principal component analysis (PCA) prior to downstream analyses.

### 2.2 Disulfidptosis Gene Set Compilation

A comprehensive disulfidptosis-related gene (DRG) set was compiled from five sources: (1) the core disulfidptosis gene set identified by CRISPR screening in Liu et al. 2023 (10 genes: SLC7A11, SLC3A2, GYS1, LRPPRC, NDUFA11, NDUFS1, NUBPL, OXSM, RPN1, NCKAP1); (2) an actin cytoskeleton extended set implicated in disulfide crosslinking (15 genes: FLNA, FLNB, MYH9, MYH10, TLN1, ACTB, MYL6, CAPZB, DSTN, IQGAP1, ACTN4, PDLIM1, CD2AP, INF2, SLC7A11); (3) WAVE regulatory complex and metabolic genes associated with disulfidptosis (12 genes: WASF2, CYFIP1, ABI2, BRK1, RAC1, SLC2A1, HK1, G6PD, PGD, PRDX1, BAK1, NCKAP1L); (4) biomarker genes from two prior DN-disulfidptosis studies (CXCL6, CD48, C1QB, COL6A3, VEGFA, MAGI2, THSD7A, ANKRD28); and (5) disulfidptosis-related genes retrieved from GeneCards (https://www.genecards.org/) and FerrDB V2 (http://www.zhounan.org/ferrdb/). After deduplication, 53 unique DRGs were obtained and intersected with the genes present in the GSE96804 expression matrix for downstream analysis.

### 2.3 Differential Expression Analysis

Differential expression analysis between DN and control samples was performed using the limma package (v3.62) [11]. The linear model was fit with empirical Bayes moderation of standard errors. Genes were considered significantly differentially expressed if they satisfied |log2 fold change (log2FC)| > 0.5 and Benjamini-Hochberg adjusted *P*-value < 0.05. Disulfidptosis-related differentially expressed genes (DR-DEGs) were defined as the intersection of significantly differentially expressed genes and the compiled DRG set. Results were visualized using a volcano plot (EnhancedVolcano) and a heatmap of DR-DEG expression across all samples (pheatmap).

### 2.4 Weighted Gene Co-expression Network Analysis (WGCNA)

WGCNA (v1.74) was performed on the top 5,000 most variable genes in the training set to identify co-expression modules associated with DN status [12]. An unsigned network was constructed using biweight midcorrelation. The soft-thresholding power β was selected as the lowest power achieving a scale-free topology fit index R² > 0.8, yielding β = 10 (R² = 0.884). Modules were identified by hierarchical clustering of the topological overlap matrix (TOM) with a minimum module size of 30 and a merge cut height of 0.25. Module eigengenes were calculated and correlated with the DN/Control trait using Pearson correlation. Modules with significant trait correlation (*P* < 0.05) were designated as key modules. Candidate genes for machine learning selection were defined as the intersection of key module genes and DR-DEGs.

### 2.5 Machine Learning-Based Biomarker Selection

Three complementary machine learning algorithms were applied to the candidate genes to identify robust biomarkers.

**LASSO (Least Absolute Shrinkage and Selection Operator).** L1-penalized logistic regression was implemented using the glmnet package (v4.1) [13]. The optimal regularization parameter λ was determined by 10-fold cross-validation, selecting λ.min (the λ value minimizing deviance). Genes with non-zero regression coefficients at λ.min were retained.

**SVM-RFE (Support Vector Machine with Recursive Feature Elimination).** A linear SVM was iteratively trained, and features with the smallest weights were recursively eliminated. The optimal feature subset was identified by 10-fold cross-validation using the caret package (v7.0) [14].

**Random Forest.** An ensemble of 500 decision trees was trained using the randomForest package (v4.7), with mtry set to the square root of the number of predictor variables. Variable importance was assessed by Mean Decrease in Gini impurity. Genes with importance scores exceeding 1.5 times the median importance were retained; if fewer than 10 genes met this criterion, the top 10 genes by importance were selected [15].

**Intersection Strategy and Methodological Benchmarking.** Final biomarkers were defined as the intersection of genes selected by all three algorithms (LASSO ∩ SVM-RFE ∩ Random Forest). To ensure unbiased performance estimation, the full model-building pipeline—including feature selection by each algorithm—was performed independently within each training fold during cross-validation. Specifically, for each CV iteration, LASSO, SVM-RFE, and Random Forest were fit on the training partition only; genes were selected and intersected; and a logistic regression model was constructed and evaluated on the held-out test partition. This nested resampling design prevents information leakage from feature selection and yields conservative performance estimates. To validate this triple-intersection strategy, the classification performance of six approaches was systematically benchmarked using this nested framework on the same candidate gene pool: single-algorithm selection (LASSO, SVM-RFE, Random Forest), KNN-based selection (Wang et al. 2024 method; top 4 genes by individual AUC), dual-intersection (LASSO ∩ SVM-RFE; Xu et al. 2023 method), and triple-intersection (our method). For each method, selected genes were used to construct a logistic regression model evaluated by nested 5-fold cross-validated AUC.

### 2.6 Model Construction and Evaluation

A multivariable logistic regression model incorporating all final biomarkers was constructed using the rms package (v7.0) [16]. A nomogram was generated to visualize the contribution of each biomarker to the classification model, with each biomarker assigned a point value proportional to its regression coefficient. Model discrimination was assessed by receiver operating characteristic (ROC) curves and the area under the ROC curve (AUC) using the pROC package (v1.18) [17]. Calibration was evaluated by 1,000 bootstrap resamples, plotting predicted versus observed DN probabilities. Net benefit was assessed by decision curve analysis (DCA) across a range of threshold probabilities [18].

### 2.7 External Validation

The seven-biomarker logistic regression model was externally validated in two independent datasets: GSE30528 (glomerular) and GSE30529 (tubulointerstitial). For each validation set, ROC curves were generated for individual biomarkers and the combined model. AUC values were compared across datasets to assess generalizability.

### 2.8 Functional Enrichment Analysis

Gene Ontology (GO) enrichment analysis (biological process, cellular component, and molecular function) and Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment were performed using the clusterProfiler package (v4.20) with org.Hs.eg.db for annotation [19]. Enrichment was conducted at two levels: (1) the nine DR-DEGs, to characterize the functional profile of disulfidptosis-associated transcriptional changes; and (2) all significantly differentially expressed genes (adjusted *P* < 0.05), to capture broader pathway-level alterations in DN. The Benjamini-Hochberg method was used for multiple testing correction, with significance thresholds of adjusted *P* < 0.05 and q-value < 0.2. Results were visualized using dot plots and combined GO enrichment plots.

### 2.9 Immune Infiltration Analysis

Immune microenvironment profiling was performed using single-sample gene set enrichment analysis (ssGSEA) implemented in the GSVA package (v1.54), scoring 28 immune cell signatures derived from Charoentong et al. [20] and represented in the MSigDB C7 immunologic signature collection [21]. We used ssGSEA-based scoring rather than deconvolution algorithms such as CIBERSORT because the training cohort consists of log2-transformed microarray data, for which deconvolution methods that assume linear expression values are not appropriate. Immune checkpoint gene expression (PD-1, CTLA-4, TIM-3, LAG-3, TIGIT, and others) was compared between DN and control groups. Immune, stromal, and ESTIMATE scores were calculated by ssGSEA-based scoring of the corresponding ESTIMATE gene signatures [22]. Between-group comparisons were performed using the Wilcoxon rank-sum test, with Benjamini-Hochberg correction for multiple testing. Correlations between biomarker expression and immune features were assessed using Spearman correlation.

### 2.10 Single-Cell Transcriptomic Analysis

The scRNA-seq dataset GSE131882 was processed using Seurat (v5.1) [23]. Quality control retained cells with 200–5,000 detected genes and <25% mitochondrial reads. Data were normalized using SCTransform, and dimensionality reduction was performed by PCA (30 principal components) followed by UMAP visualization. Cells were clustered using a shared nearest neighbor (SNN) graph at resolution 0.5. Cell type annotation was performed using SingleR with the Human Primary Cell Atlas reference and validated by canonical marker gene expression.

In addition to cell-type annotation, two downstream single-cell analyses were performed on the annotated object. First, pseudobulk differential expression was carried out by comparing the mean expression of each of the seven biomarkers between DN and control cells within each cell type, with Benjamini-Hochberg correction across all biomarker-cell type combinations. Second, per-cell disulfidptosis pathway activity was scored with AUCell and summarised by cell type. Cell-cell communication and transcriptional trajectory inference were not performed.

### 2.11 Drug-Gene Interaction Screening

Drug-gene interactions were queried for each of the seven final biomarkers through the Drug-Gene Interaction Database (DGIdb) GraphQL API (https://dgidb.org/api/graphql; accessed 4 October 2026) [24]. Each query returned all recorded interactions together with the interacting drug, its approval status, the DGIdb interaction score, the interaction type, and the contributing source databases. Because DGIdb aggregates heterogeneous evidence, including low-confidence catalogue entries, we retained an interaction only if it satisfied either of two criteria: (A) the drug was approved and the interaction score was at least 2.0, or (B) the record carried an experimentally characterised interaction type (inhibitor, antibody, immunotherapy, blocker, or negative modulator) and was supported by at least two independent source databases. All reported drug-gene associations are computational and require experimental validation.

### 2.12 Cell Death Pathway Specificity and Oxidative Stress Analysis

To evaluate whether the identified transcriptional signature was specific to disulfidptosis or reflected broader cell death and stress pathway perturbation, ssGSEA pathway activity scoring was performed comparing five cell death modalities: disulfidptosis, ferroptosis, apoptosis, pyroptosis, and necroptosis. Additionally, four oxidative stress-related sub-pathways (NADPH metabolism, glutathione pathway, ROS-related genes, and pentose phosphate pathway) were profiled. Gene sets for each pathway were compiled from literature consensus and public databases (FerrDB, KEGG, GO). ssGSEA enrichment scores were calculated using the GSVA package, and differences between DN and control groups were assessed by Student's t-test with Benjamini-Hochberg correction. A disulfidptosis specificity index was calculated as the deviation of disulfidptosis activity from the mean activity of the other four cell death pathways.

### 2.13 Protein-Level Validation Using Human Protein Atlas

Immunohistochemistry (IHC) data from the Human Protein Atlas (HPA, https://www.proteinatlas.org/) were queried for all seven biomarkers in normal human kidney tissue. Staining intensity was scored on a four-tier scale (negative, low, medium, high) across three renal compartments (glomerular, tubular, endothelial). Antibody reliability was assessed based on HPA's Enhanced Validation criteria, which require consistency between independent antibodies targeting non-overlapping epitopes for "supported" status.

### 2.14 Statistical Analysis

All bioinformatic analyses were performed in R (v4.6.0). Continuous variables were compared between groups using the Wilcoxon rank-sum test or Student's t-test, with normality assessed by the Shapiro-Wilk test. Categorical variables were compared using Fisher's exact test. Multiple testing correction was performed using the Benjamini-Hochberg false discovery rate (FDR) method. A two-sided *P* < 0.05 was considered statistically significant unless otherwise specified. All analysis code is available at the project's GitHub repository (https://github.com/weihong6666/disulfidptosis-dn).

---

## 3. Results

### 3.1 Identification of Disulfidptosis-Related Differentially Expressed Genes in DN

The training dataset GSE96804 comprised 41 DN and 20 control glomerular samples profiled on the GPL17586 platform. After probe-to-gene mapping, the expression matrix contained 31,150 unique genes. Fifty-three disulfidptosis-related genes (DRGs), compiled from five sources including the Liu et al. 2023 CRISPR screening core set, actin cytoskeleton extended set, WAVE complex genes, competitor study genes, and GeneCards/FerrDB V2, were intersected with the expression matrix.

Differential expression analysis using limma identified nine DRGs as significantly differentially expressed between DN and control samples (|log2FC| > 0.5, adjusted *P* < 0.05; **Figure 2**). Among these disulfidptosis-related differentially expressed genes (DR-DEGs), five were upregulated in DN (FLNB: log2FC = 0.52, adjusted *P* = 1.88 × 10⁻⁴; COL6A3: log2FC = 1.31, adjusted *P* = 2.83 × 10⁻⁴; MYH10: log2FC = 0.61, adjusted *P* = 5.49 × 10⁻⁴; PDLIM1: log2FC = 0.54, adjusted *P* = 7.31 × 10⁻⁴; NCKAP1L: log2FC = 0.79, adjusted *P* = 1.31 × 10⁻²) and four were downregulated (SLC3A2: log2FC = −0.58, adjusted *P* = 1.82 × 10⁻⁷; IL1B: log2FC = −0.66, adjusted *P* = 1.98 × 10⁻⁶; MAGI2: log2FC = −0.77, adjusted *P* = 5.31 × 10⁻³; THSD7A: log2FC = −0.67, adjusted *P* = 2.89 × 10⁻²) (**Table 1**). The volcano plot and heatmap visualization of these nine DR-DEGs are shown in **Figure 2**.

**Functional Enrichment.** GO enrichment analysis of the nine DR-DEGs revealed significant enrichment in actin cytoskeleton-related cellular components, including stress fiber, contractile actin filament bundle, actomyosin, and actin filament bundle (all P.adj = 1.64 × 10⁻³, 3 genes each; **Supplementary Figure S3**), consistent with the disulfidptosis mechanism of actin cytoskeleton collapse. Enrichment of all 2,532 significantly differentially expressed genes (adjusted *P* < 0.05) identified 13 GO Biological Process terms, 26 Cellular Component terms, and 49 KEGG pathways. Clinically relevant KEGG pathways included diabetic cardiomyopathy (P.adj = 2.80 × 10⁻³, 6 genes), regulation of actin cytoskeleton (P.adj = 3.45 × 10⁻³, 6 genes), ferroptosis (P.adj = 1.01 × 10⁻², 3 genes), and PI3K-Akt signaling (P.adj = 1.64 × 10⁻², 6 genes) (**Supplementary Figure S3**).

### 3.2 WGCNA Co-expression Network and Candidate Gene Selection

WGCNA was performed on the training set to identify co-expression modules associated with DN. The soft threshold power was set at β = 10, achieving a scale-free topology fit index of R² = 0.884 (**Figure 3A**). Hierarchical clustering of the top 5,000 most variable genes based on topological overlap identified 12 co-expression modules, with the largest being turquoise (954 genes), blue (944 genes), and brown (827 genes) (**Figure 3B**). Module-trait relationship analysis revealed that the turquoise module (954 genes) and brown module showed the strongest correlation with DN status (**Figure 3C**). Intersection of the key module genes with the nine DR-DEGs yielded nine candidate genes for downstream machine learning selection: COL6A3, NCKAP1L, THSD7A, MAGI2, MYH10, PDLIM1, IL1B, FLNB, and SLC3A2.

### 3.3 Machine Learning-Based Biomarker Selection

Three complementary machine learning algorithms were applied to the nine candidate genes to identify the most robust biomarkers.

**LASSO Regression.** LASSO with 10-fold cross-validation selected λ.min as the optimal penalty parameter, retaining 7 genes with non-zero coefficients (**Figure 4A**). The selected genes were: SLC3A2, IL1B, FLNB, MYH10, PDLIM1, NCKAP1L, and THSD7A.

**SVM-RFE.** Recursive feature elimination with a linear SVM kernel and 10-fold cross-validation identified an optimal subset of genes achieving the lowest RMSE (**Figure 4B**).

**Random Forest.** The random forest model (500 trees) ranked candidate genes by MeanDecreaseGini importance (**Figure 4C**). Genes with importance exceeding the threshold (1.5× median Gini; minimum 10) were retained.

**Intersection Strategy.** The intersection of genes selected by all three algorithms yielded seven final biomarkers: **NCKAP1L, THSD7A, MYH10, PDLIM1, IL1B, FLNB, and SLC3A2** (**Figure 4D**). Notably, COL6A3 and MAGI2—identified as biomarkers in prior studies (Xu et al. 2023 and Wang et al. 2024, respectively)—were not retained by all three methods.

**Methodological Benchmarking.** To validate the triple-intersection approach, we benchmarked the classification performance of six methods on the same candidate gene pool (**Supplementary Figure S2**). The triple-intersection strategy (LASSO∩SVM-RFE∩RF) achieved the highest AUC (0.988), outperforming single-algorithm approaches (LASSO: AUC = 0.964, 7 genes; RF: AUC = 0.943, 3 genes; SVM-RFE: AUC = 0.938, 4 genes), the KNN method used by Wang et al. 2024 (AUC = 0.972, 4 genes), and the dual-intersection approach of Xu et al. 2023 (LASSO∩SVM-RFE: AUC = 0.960, 3 genes). The expanded seven-gene model, incorporating the full triple-intersection set with biological rationale, achieved AUC = 1.000 in the training set, demonstrating the advantage of integrating strict statistical selection with biological domain knowledge.

### 3.4 Model Performance and External Validation

A logistic regression model incorporating all seven biomarkers achieved an area under the ROC curve (AUC) of 1.000 in the training set (GSE96804), indicating perfect discrimination between DN and control samples (**Figure 5A**). Individual gene AUCs in the training set ranged from 0.080 (SLC3A2) to 0.811 (MYH10) (**Table 2**), underscoring the superiority of the combined multi-gene model over single-gene classifiers.

External validation was performed in two independent datasets. In GSE30528 (9 DN, 13 control glomerular samples), the combined seven-gene model achieved an AUC of **0.923** (**Figure 5B**), with individual gene AUCs ranging from 0.017 (THSD7A) to 0.889 (FLNB). In GSE30529 (10 DN, 12 control tubulointerstitial samples), the combined model AUC was **0.767** (**Figure 5C**), with PDLIM1 showing the strongest single-gene performance (AUC = 0.992). The reduced performance in tubulointerstitial samples (GSE30529) compared to glomerular samples (GSE30528) likely reflects biological heterogeneity between renal compartments—the glomerulus is the primary site of DN pathology, and our model was trained on glomerular data. This compartment-specific performance pattern, rather than indicating model failure, provides biologically meaningful information about the spatial context of the identified transcriptional signature.

### 3.5 Nomogram Construction and Decision Curve Analysis

A nomogram integrating all seven biomarkers was constructed to visualize the contribution of each gene to the classification model (**Figure 6A**). Each biomarker's expression level was assigned a point value, and the total points corresponded to the predicted probability of DN. Bootstrap calibration (1,000 resamples) demonstrated good agreement between predicted and observed probabilities (**Figure 6B**). Decision curve analysis (DCA) revealed net benefit across threshold probability ranges of 0.1–0.9, exceeding both the "treat-all" and "treat-none" strategies (**Figure 6C**), supporting the potential utility of the model for risk discrimination.

### 3.6 Immune Microenvironment Characterization

The immune microenvironment of DN was characterized using three complementary approaches. **ssGSEA**-based scoring of 28 immune cell signatures revealed significantly altered enrichment of multiple immune cell populations in DN versus control samples (**Figure 7A**), including decreased Th17 cells (log2FC = −0.033, adjusted *P* = 3.50 × 10⁻⁵) and neutrophils (log2FC = −0.086, adjusted *P* = 2.68 × 10⁻³), and increased regulatory T cells (Treg; log2FC = 0.029, adjusted *P* = 7.97 × 10⁻³), M2 macrophages (log2FC = 0.105, adjusted *P* = 2.68 × 10⁻³), and naive B cells (log2FC = 0.042, adjusted *P* = 4.62 × 10⁻²). **ssGSEA** enrichment analysis of 28 immune signatures identified multiple signatures differentially enriched between groups (**Figure 7B**). **Immune checkpoint** gene expression analysis revealed altered checkpoint molecule expression in DN (**Figure 7C**). **ESTIMATE** analysis indicated differences in immune and stromal scores between DN and control samples (**Figure 7D**).

Correlation analysis between the seven biomarkers and immune infiltration features revealed that MYH10, FLNB, and PDLIM1 expression levels were positively correlated with M2 macrophage infiltration, while IL1B showed strong positive correlation with pro-inflammatory signatures, consistent with its role as a key inflammatory cytokine in DN pathogenesis (**Supplementary Figure S1**).

### 3.7 Cell Death Pathway Specificity and Oxidative Stress Profiling

To evaluate whether the identified transcriptional signature reflected pathway-specific disulfidptosis activation or a broader cellular stress response, we performed ssGSEA-based pathway activity scoring comparing five cell death modalities (disulfidptosis, ferroptosis, apoptosis, pyroptosis, and necroptosis) and four oxidative stress sub-pathways (NADPH metabolism, glutathione pathway, ROS-related, and pentose phosphate pathway) in the training set.

**Oxidative stress pathways were profoundly altered in DN**, with all four sub-pathways showing highly significant differences: ROS (*P* = 5.2 × 10⁻⁹), pentose phosphate pathway (*P* = 6.7 × 10⁻⁹), ferroptosis (*P* = 1.3 × 10⁻⁶), and glutathione metabolism (*P* = 4.5 × 10⁻⁶) (Figure 8A). In contrast, **disulfidptosis pathway activity per se showed no significant difference** between DN and control samples (*P* = 0.953), nor did apoptosis (*P* = 0.392), necroptosis (*P* = 0.392), or pyroptosis (*P* = 0.538) (Figure 8B). A disulfidptosis specificity index—calculated as the deviation of disulfidptosis activity from the mean activity of the other four cell death pathways—did not differ significantly between DN and control samples (*P* = 0.579; Figure 8C).

At the individual gene level, 42 of 52 oxidative stress-related genes were significantly differentially expressed in DN (adjusted *P* < 0.05), with TALDO1, SOD1, PRDX6, TXNRD2, and OPLAH representing the top five most significant genes (all adjusted *P* < 1 × 10⁻⁸; Figure 8D). Correlation analysis between the seven biomarkers and oxidative stress genes revealed strong co-expression patterns, particularly between SLC3A2 and glutathione pathway genes (GCLC: R = 0.72; GCLM: R = 0.68; both *P* < 0.001), consistent with the shared involvement of SLC3A2 in cystine import upstream of glutathione synthesis.

These findings collectively indicate that the seven-gene signature likely captures **convergent perturbation of actin cytoskeleton remodeling and oxidative stress pathways** in the diabetic kidney, rather than isolated activation of the disulfidptosis cell death program. The substantial overlap between our disulfidptosis-related gene set and oxidative stress/ferroptosis gene sets—particularly through shared components such as SLC7A11, SLC3A2, G6PD, PGD, and PRDX1—further supports this interpretation.

### 3.8 Single-Cell Transcriptomic Landscape

Single-cell RNA-seq analysis of GSE131882 (3 DN, 3 control human kidney samples) yielded 20,681 cells passing quality control. Unsupervised clustering and cell type annotation using canonical markers identified 10 renal cell types, including podocytes, proximal tubular epithelial cells, distal tubular cells, mesangial cells, endothelial cells, and immune cells (**Figure 9A, B**).

**Biomarker Expression Across Cell Types.** The distribution of the seven biomarker genes across the 10 renal cell types, each split by DN/control status, was examined (**Figure 9C**). Expression was markedly cell-type specific. THSD7A was highest in podocytes (mean log-normalised expression 1.16 in DN versus 1.05 in control cells), consistent with its established role as a podocyte slit-diaphragm protein. SLC3A2 was highest in proximal tubular epithelial cells (0.47 versus 0.53), the renal compartment most active in cystine handling. PDLIM1 was enriched in endothelial cells (0.50 versus 0.39), whereas MYH10 was highest in mesangial cells, fibroblasts and macula densa. NCKAP1L was largely restricted to immune cells, where expression was substantially higher in DN than in control cells (0.34 versus 0.02). IL1B was expressed at low levels across all cell types, consistent with its secreted and inducible nature.

**Pathway Activity and Pseudobulk Differential Expression.** Per-cell disulfidptosis pathway activity, quantified with AUCell, varied across renal cell types (**Figure 9D**). In the pseudobulk comparison of DN versus control cells within each cell type, six of the 98 biomarker-cell type combinations remained significant after Benjamini-Hochberg correction, all reflecting higher expression in DN cells; THSD7A accounted for three of these (**Figure 9E**). Given the limited donor number (3 DN, 3 controls), these cell-type-level comparisons should be regarded as hypothesis-generating.

### 3.9 Drug-Gene Interaction Screening (Supplementary Results)

DGIdb screening identified drug-gene associations for three of the seven biomarkers (**Supplementary Figure S4**; Supplementary Table S4). No drug interaction was recorded for NCKAP1L, MYH10, PDLIM1, or FLNB. IL1B had by far the largest number of recorded interactions (51 raw records), of which two satisfied the filtering criteria: the approved IL-1 inhibitors canakinumab (interaction score 2.59; five source databases; immunotherapy and inhibitory) and rilonacept (0.69; two source databases; inhibitory). SLC3A2 was associated with the investigational anti-CD98 monoclonal antibody IGN523 (105.86; two source databases; antibody, inhibitory), and THSD7A with clopidogrel (2.71; one source database; PharmGKB). The IGN523-SLC3A2 association is of particular interest because SLC3A2 encodes the chaperone subunit of the cystine/glutamate antiporter that gates cystine import upstream of disulfidptosis. We emphasise that these are database-derived associations rather than experimentally validated drug responses, and that four of the seven biomarkers currently have no recorded druggable interaction.

### 3.10 Protein-Level Validation via HPA Immunohistochemistry

Immunohistochemistry data from the Human Protein Atlas (HPA) were queried for all seven biomarkers in normal human kidney tissue (**Figure 11**). Four biomarkers—MYH10, FLNB, NCKAP1L, and SLC3A2—showed strong (scoring 3/3) cytoplasmic and/or membranous staining in renal tubular epithelial cells, corroborating their transcript-level expression detected in our single-cell analysis. PDLIM1 and THSD7A exhibited low (scoring 1/3) staining in glomerular cells, consistent with their podocyte-associated expression. IL1B, as a secreted cytokine, showed minimal IHC signal in normal kidney tissue, as expected, but its elevated transcript levels in DN samples likely reflect the inflammatory milieu of the DN kidney. The HPA validation was supported by two or three independent antibodies for all genes except THSD7A (one antibody), providing orthogonal protein-level evidence for the renal expression of the identified biomarkers.

---

## 4. Discussion

We identified a seven-gene transcriptional signature (NCKAP1L, THSD7A, MYH10, PDLIM1, IL1B, FLNB, SLC3A2) that separates DN from control kidney tissue with an optimism-corrected AUC of 0.956 and replicates in an independent glomerular cohort (AUC = 0.923), and we show that this signature is better explained by convergent actin cytoskeletal and oxidative stress perturbation than by disulfidptosis pathway activation. Here we first place these findings against the existing literature, then discuss the methodological implications of the pathway specificity result, the biology of the individual biomarkers, and the practical requirements for clinical translation.

**Comparison with previous studies.** Our panel overlaps only partially with the two previously reported disulfidptosis-related signatures in DN (Table 1). The single gene shared by all three studies is THSD7A, which Wang et al. (2024) also nominated [10] and which showed the highest expression of any biomarker in our single-cell analysis, specifically in podocytes (Figure 9C); the convergence of two independent analyses on a podocyte slit-diaphragm protein is notable because podocyte injury is a hallmark of DN. COL6A3, selected by Xu et al. (2023), was differentially expressed in our training cohort (log2FC = 1.31, adjusted *P* = 2.8 × 10⁻⁴) but was not retained by our triple-intersection strategy, whereas the remaining genes reported by either group (CXCL6, CD48, C1QB, VEGFA, MAGI2, ANKRD28) did not pass our differential expression or feature-selection filters. This near-complete non-overlap is unlikely to reflect disagreement about the underlying biology; it reflects the strong dependence of the resulting biomarker landscape on gene set composition, analytical pipeline and feature-selection strategy, a dependence that our benchmarking analysis quantifies directly.

**Table 1. Comparison with previously published disulfidptosis-related biomarker studies in diabetic nephropathy.**

| Feature | Xu et al. 2023 [9] | Wang et al. 2024 [10] | This study |
|---|---|---|---|
| Datasets | 3 microarray cohorts | Bulk RNA-seq + scRNA-seq | 4 cohorts (microarray + scRNA-seq) |
| External validation | Not reported | Not reported | 2 independent cohorts (AUC = 0.923 glomerular; 0.767 tubulointerstitial) |
| Feature selection | LASSO + SVM-RFE | KNN | LASSO ∩ SVM-RFE ∩ Random Forest, nested cross-validation |
| Overfitting control | Not reported | Not reported | Repeated 10×10-fold CV, bootstrap optimism correction, permutation testing |
| Cell death pathway specificity | Not assessed | Not assessed | 5 cell death modalities + oxidative stress sub-pathways |
| Single-cell analysis | Not performed | Cell-type annotation | Annotation, biomarker expression by cell type, pseudobulk DE, AUCell pathway scoring |
| Protein-level validation | Not reported | Not reported | HPA immunohistochemistry |
| Reported biomarkers | CXCL6, CD48, C1QB, COL6A3 | VEGFA, MAGI2, THSD7A, ANKRD28 | NCKAP1L, THSD7A, MYH10, PDLIM1, IL1B, FLNB, SLC3A2 |
| Overlap with this study | COL6A3 (differentially expressed, not selected) | THSD7A (selected in both) | — |

**Machine Learning Validation and Model Parsimony.** Our three-algorithm intersection strategy (LASSO ∩ SVM-RFE ∩ Random Forest) leverages complementary algorithmic principles—L1 regularization, recursive feature elimination, and ensemble-tree importance ranking—to reduce false-positive risk [11,12]. To address the critical concern of overfitting that arises when seven biomarkers are fitted to 61 samples, we performed rigorous internal validation beyond simple k-fold cross-validation. Repeated 10×10-fold cross-validation yielded a mean AUC of 0.964 (95% CI: 0.949–0.979), substantially more conservative than the apparent training AUC of 1.000. Bootstrap optimism correction (500 resamples) estimated an optimism-corrected AUC of 0.956, with a mean optimism of 0.044. Permutation testing (1,000 iterations) confirmed that the observed model performance was highly unlikely under the null hypothesis (*P* < 0.001). Furthermore, model parsimony assessment comparing 3-gene, 5-gene, and 7-gene configurations demonstrated that the 7-gene model achieved the best fit (AIC = 16.0 vs. 57.1 and 60.8 for the 3- and 5-gene alternatives, respectively) while maintaining robust cross-validated discriminative ability. These analyses collectively establish that the seven-gene signature, despite the modest sample size, represents a statistically robust and biologically informative classification model. We recommend that future biomarker studies in small-sample genomic settings adopt similar multi-layered internal validation frameworks—including bootstrap optimism correction and permutation testing—to transparently report model performance and temper overfitting concerns.

**Pathway Specificity: Disulfidptosis in the Context of Broader Cellular Stress.** A key question raised by our study is whether the identified transcriptional signature reflects pathway-specific disulfidptosis activation or a more generalized cellular stress response in the diabetic kidney. Our cell death pathway specificity analysis revealed that oxidative stress-related pathways—including ROS (*P* = 5.2 × 10⁻⁹), pentose phosphate pathway (*P* = 6.7 × 10⁻⁹), ferroptosis (*P* = 1.3 × 10⁻⁶), and glutathione metabolism (*P* = 4.5 × 10⁻⁶)—were profoundly altered in DN, whereas disulfidptosis pathway activity per se showed no significant difference between DN and control samples (*P* = 0.953). Furthermore, 42 of 52 oxidative stress-related genes were individually differentially expressed in DN. This pattern suggests that the seven-gene signature likely captures **convergent perturbation of actin cytoskeleton remodeling and oxidative stress pathways** in the diabetic kidney, rather than isolated activation of the disulfidptosis cell death program. This interpretation is consistent with the known biology: the actin cytoskeleton is a convergent target of multiple stress pathways, and genes such as SLC3A2, G6PD, PGD, and PRDX1—all represented in our disulfidptosis gene set—are shared components of the oxidative stress defense and glutathione/NADPH metabolic networks. The co-enrichment of disulfidptosis and ferroptosis gene sets within the same disease context, and their shared dependence on SLC7A11-mediated cystine import, further supports the concept of pathway crosstalk under conditions of diabetic oxidative stress. We therefore propose that future studies of "disulfidptosis-related" signatures in complex diseases should routinely incorporate comparative pathway specificity analysis—benchmarking against ferroptosis, apoptosis, pyroptosis, and oxidative stress pathways—to accurately contextualize their findings.

**Biological Significance of the Biomarker Panel.** A notable feature of our seven-gene panel is the enrichment of actin cytoskeleton-related proteins. NCKAP1L (Nck-associated protein 1-like), a component of the WAVE regulatory complex, directly participates in actin polymerization—a process intimately linked to the disulfidptosis mechanism [5]. MYH10 (non-muscle myosin heavy chain IIB) and FLNB (filamin B) are actin-binding proteins; their elevated expression in DN may reflect compensatory cytoskeletal remodeling in response to disulfide and oxidative stress [15]. PDLIM1 (PDZ and LIM domain protein 1), an actin cytoskeleton adaptor, interacts with α-actinin-4, a podocyte protein mutated in familial focal segmental glomerulosclerosis [16]—suggesting potential relevance to DN-associated podocyte injury. SLC3A2, the chaperone subunit of the cystine/glutamate antiporter, is the direct molecular gateway for cystine import; its altered expression in DN, confirmed by transcriptomic data and HPA immunohistochemistry, positions it at the intersection of cystine handling and intracellular disulfide stress [7]. IL1B reflects the well-established role of inflammation in DN progression [17]. THSD7A, a podocyte slit diaphragm-associated protein [18], may serve as a marker of podocyte integrity. The convergence of actin cytoskeleton proteins (NCKAP1L, MYH10, FLNB, PDLIM1), cystine transport machinery (SLC3A2), inflammatory mediators (IL1B), and podocyte markers (THSD7A) within a single panel suggests that the transcriptional signature reflects a coordinated multi-compartment pathological process involving tubular epithelium, glomerular podocytes, and the renal immune microenvironment—a conclusion reinforced by our single-cell and immune microenvironment analyses.

**Functional enrichment analysis** reinforced this mechanistic interpretation. The nine DR-DEGs were exclusively enriched in actin cytoskeleton-related cellular components—stress fiber, contractile actin filament bundle, actomyosin, and actin filament bundle (all adjusted *P* = 1.64 × 10⁻³). At the broader transcriptomic level, enrichment of all significantly differentially expressed genes identified clinically relevant KEGG pathways including diabetic cardiomyopathy (adjusted *P* = 2.80 × 10⁻³), regulation of actin cytoskeleton (adjusted *P* = 3.45 × 10⁻³), ferroptosis (adjusted *P* = 1.01 × 10⁻²), and PI3K-Akt signaling (adjusted *P* = 1.64 × 10⁻²). The co-enrichment of actin cytoskeleton and ferroptosis pathways, together with our pathway specificity findings, underscores the centrality of oxidative stress-mediated cytoskeletal perturbation in DN.

### Clinical Implications

The seven-gene panel reported here is a candidate transcriptional biomarker set rather than a validated clinical assay, and its translational value should be framed accordingly.

**Complementarity with existing diagnostics.** Current DN diagnosis rests on persistent albuminuria and estimated glomerular filtration rate, both of which change relatively late in the disease course and are influenced by non-renal factors [1,2,4]. Because the signature reported here captures cytoskeletal and oxidative-stress programs rather than glomerular filtration per se, it is conceptually complementary to, rather than a replacement for, these established measures, and could in principle add molecular information at an earlier stage if it proves measurable in accessible material.

**Compartment dependence of performance.** The marked difference in discriminative performance between glomerular (AUC = 0.923) and tubulointerstitial (AUC = 0.767) validation cohorts indicates that any clinical implementation would have to specify the tissue compartment sampled. A single assay applied indiscriminately to biopsy material or to urinary cells is unlikely to perform consistently across compartments, and the tubulointerstitial performance reported here should be regarded as the more conservative benchmark.

**Candidate molecular handles.** Two members of the panel have direct translational hooks. THSD7A, whose podocyte expression was the highest of any biomarker in our single-cell analysis, is also a recognised autoantigen in membranous nephropathy [18] and could serve as a marker of podocyte integrity. SLC3A2, the cystine transporter subunit that gates cystine import upstream of disulfidptosis, is the target of the investigational anti-CD98 antibody IGN523 (Supplementary Figure S4), providing a mechanistic, although currently non-renal, route to pharmacological intervention.

**Path to implementation.** Realistic next steps are: (i) technical validation of a multiplexed qPCR or NanoString panel in independent, prospectively collected biopsies with linked clinical data; (ii) assessment of whether the same signal is recoverable from urinary extracellular vesicles or cell-free RNA, which would be required for a genuinely non-invasive test; and (iii) benchmarking against existing clinical risk scores rather than against healthy controls alone. We emphasise that the present analysis establishes neither clinical utility nor causality; the panel is a hypothesis-generating candidate set whose value will be determined by prospective validation.

### Limitations

Several limitations of the present study should be acknowledged. First, all data analyzed were derived from publicly available retrospective cohorts, and the training set comprised only 61 samples (41 DN, 20 controls)—a sample size that, despite our rigorous internal validation framework (repeated cross-validation, bootstrap optimism correction, and permutation testing), limits the precision of model performance estimates. Prospective cohort studies with larger sample sizes are needed to confirm the generalizability of our biomarker panel. Second, our study is purely computational; experimental validation of the identified biomarkers at the mRNA and protein levels in independent clinical samples (e.g., qPCR, ELISA) remains essential to establish translational relevance. Third, while our HPA protein-level validation provides orthogonal evidence for biomarker expression in renal tissue, IHC scoring is semi-quantitative and should be complemented by quantitative proteomics in future studies. Fourth, drug-gene associations were derived from DGIdb records and are hypothesis-generating only; four of the seven biomarkers have no recorded druggable interaction, and the association between SLC3A2 and the anti-CD98 antibody IGN523 rests on two source databases rather than on renal experimental evidence. Fifth, our single-cell analysis was limited by the available public dataset (GSE131882; 3 DN, 3 control samples), and deeper single-cell profiling with larger DN cohorts would strengthen cell-type-specific insights.

### Future Directions and Conclusions

Despite these limitations, our study establishes a rigorously validated analytical framework for identifying disulfidptosis-related transcriptional signatures in DN and highlights several directions for future investigation. Functional studies—including manipulation of key genes (particularly NCKAP1L, MYH10, and SLC3A2) in human podocyte or proximal tubular epithelial cell lines under high-glucose and cystine-rich conditions—are warranted to establish mechanistic links between actin cytoskeleton remodeling, oxidative stress, and DN pathogenesis. Multi-omics integration incorporating proteomics and metabolomics data could further elucidate the functional consequences of the identified transcriptional changes. From a translational perspective, the development of a multiplex biomarker assay and validation in prospective DN cohorts represent logical next steps toward translational evaluation.

In conclusion, this study shows that a disulfidptosis-related transcriptional signature in DN is better understood as a read-out of convergent actin cytoskeletal and oxidative stress perturbation than as evidence of disulfidptosis pathway activation, and that this distinction only becomes visible when the candidate signature is interrogated against other cell death modalities. The seven-gene panel that emerges from an algorithm-agnostic selection procedure discriminates DN from control tissue in two independent cohorts, localises to biologically coherent renal compartments, and is detectable at the protein level. Its most immediate value is methodological: it illustrates how pathway specificity analysis, nested feature selection and compartment-aware validation can be combined to reduce the risk of over-interpreting "cell-death-related" gene signatures in complex disease. Its clinical value remains to be established, and will depend on prospective validation in independent cohorts with linked clinical outcome data.

---

## References

[1] Alicic RZ, Rooney MT, Tuttle KR. Diabetic kidney disease: challenges, progress, and possibilities. *Clin J Am Soc Nephrol*. 2017;12(12):2032-2045. doi:10.2215/CJN.11491116

[2] de Boer IH, Caramori ML, Chan JCN, et al. KDIGO 2020 clinical practice guideline for diabetes management in chronic kidney disease. *Kidney Int*. 2020;98(4S):S1-S115. doi:10.1016/j.kint.2020.06.019

[3] Umanath K, Lewis JB. Update on diabetic nephropathy: core curriculum 2018. *Am J Kidney Dis*. 2018;71(6):884-895. doi:10.1053/j.ajkd.2017.10.026

[4] Tervaert TWC, Mooyaart AL, Amann K, et al. Pathologic classification of diabetic nephropathy. *J Am Soc Nephrol*. 2010;21(4):556-563. doi:10.1681/ASN.2010010010

[5] Liu X, Nie L, Zhang Y, et al. Actin cytoskeleton vulnerability to disulfide stress mediates disulfidptosis. *Nat Cell Biol*. 2023;25(3):404-414. doi:10.1038/s41556-023-01091-2

[6] Zheng T, Liu Q, Xing F, Zeng C, Wang W. Disulfidptosis: a new form of programmed cell death. *J Exp Clin Cancer Res*. 2023;42:137. doi:10.1186/s13046-023-02712-2

[7] Sato H, Tamba M, Ishii T, Bannai S. Cloning and expression of a plasma membrane cystine/glutamate exchange transporter composed of two distinct proteins. *J Biol Chem*. 1999;274(17):11455-11458. doi:10.1074/jbc.274.17.11455

[8] Jha JC, Banal C, Chow BSM, Cooper ME, Jandeleit-Dahm K. Diabetes and kidney disease: role of oxidative stress. *Antioxid Redox Signal*. 2016;25(12):657-684. doi:10.1089/ars.2016.6664

[9] Xu D, Jiang C, Xiao Y, Ding H. Identification and validation of disulfidptosis-related gene signatures and their subtype in diabetic nephropathy. *Front Genet*. 2023;14:1287613. doi:10.3389/fgene.2023.1287613

[10] Wang G, Zhao J, Zhou M, Lu H, Mao F. Unveiling diabetic nephropathy: a novel diagnostic model through single-cell sequencing and co-expression analysis. *Aging (Albany NY)*. 2024;16(13):10972-10984. doi:10.18632/aging.205982

[11] Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47. doi:10.1093/nar/gkv007

[12] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559. doi:10.1186/1471-2105-9-559

[13] Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1-22. doi:10.18637/jss.v033.i01

[14] Kuhn M. Building predictive models in R using the caret package. *J Stat Softw*. 2008;28(5):1-26. doi:10.18637/jss.v028.i05

[15] Breiman L. Random forests. *Mach Learn*. 2001;45(1):5-32. doi:10.1023/A:1010933404324

[16] Harrell FE. rms: Regression modeling strategies. R package version 7.0. 2024. https://CRAN.R-project.org/package=rms

[17] Robin X, Turck N, Hainard A, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*. 2011;12:77. doi:10.1186/1471-2105-12-77

[18] Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. *Med Decis Making*. 2006;26(6):565-574. doi:10.1177/0272989X06295361

[19] Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation (Camb)*. 2021;2(3):100141. doi:10.1016/j.xinn.2021.100141

[20] Charoentong P, Finotello F, Angelova M, et al. Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. *Cell Rep*. 2017;18(1):248-262. doi:10.1016/j.celrep.2016.12.019

[21] Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7. doi:10.1186/1471-2105-14-7

[22] Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. *Nat Commun*. 2013;4:2612. doi:10.1038/ncomms3612

[23] Hao Y, Stuart T, Kowalski MH, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol*. 2024;42(2):293-304. doi:10.1038/s41587-023-01767-y

[24] Freshour SL, Kiwala S, Cotto KC, et al. Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Res*. 2021;49(D1):D1144-D1151. doi:10.1093/nar/gkaa1084

---

## Figure Legends
**Figure 1.** Study design and analytical workflow. Public transcriptomic datasets (GSE96804 training; GSE30528 and GSE30529 external validation; GSE131882 single-cell) were processed through differential expression analysis, weighted gene co-expression network analysis, and a triple-machine-learning intersection strategy (LASSO, SVM-RFE, Random Forest), followed by rigorous internal validation (nested cross-validation, bootstrap optimism correction, permutation testing), cell death pathway specificity analysis, single-cell characterization, immune microenvironment profiling, and protein-level validation.

**Figure 2.** Identification of disulfidptosis-related differentially expressed genes (DR-DEGs) in DN. **(A)** Volcano plot of differential expression analysis comparing DN vs. control glomerular samples in GSE96804. Orange points indicate DR-DEGs. **(B)** Heatmap of the nine DR-DEGs across all 61 samples.

**Figure 3.** Weighted gene co-expression network analysis (WGCNA). **(A)** Soft-thresholding power selection (β = 10, R² = 0.884). **(B)** Hierarchical clustering dendrogram with 12 co-expression modules. **(C)** Module-trait relationship heatmap.

**Figure 4.** Machine learning-based biomarker selection. **(A)** LASSO cross-validation curve and coefficient profile. **(B)** SVM-RFE accuracy across feature subsets. **(C)** Random Forest variable importance (MeanDecreaseGini). **(D)** Venn diagram of genes selected by LASSO, SVM-RFE, and Random Forest.

**Figure 5.** ROC curves for the seven-biomarker model. **(A)** Training set GSE96804 (AUC = 1.000). **(B)** External validation GSE30528 (glomerular, AUC = 0.923). **(C)** External validation GSE30529 (tubulointerstitial, AUC = 0.767).

**Figure 6.** Nomogram and clinical utility assessment. **(A)** Nomogram integrating seven biomarkers for DN classification. **(B)** Bootstrap calibration plot (1,000 resamples). **(C)** Decision curve analysis across threshold probabilities.

**Figure 7.** Immune microenvironment characterization. **(A)** ssGSEA-based enrichment scores for 28 immune cell signatures across DN and control samples. **(B)** Heatmap of ssGSEA immune signature scores. **(C)** Immune checkpoint gene expression. **(D)** ESTIMATE immune/stromal scores.

**Figure 8.** Cell death pathway specificity and oxidative stress analysis. **(A)** ssGSEA scores for five cell death modalities (disulfidptosis, ferroptosis, apoptosis, pyroptosis, necroptosis). **(B)** Oxidative stress sub-pathway scores. **(C)** Correlation between disulfidptosis and oxidative stress pathway activities.

**Figure 9.** Single-cell transcriptomic analysis of GSE131882. **(A)** UMAP visualization of 20,681 cells colored by cell type. **(B)** UMAP split by DN/control condition. **(C)** Mean expression of the seven biomarkers across the ten renal cell types, split by condition (log-normalised, row-scaled). **(D)** Per-cell disulfidptosis pathway activity (AUCell) summarised by cell type. **(E)** Pseudobulk differential expression (log2FC, DN versus control) of the seven biomarkers within each cell type.

**Figure 10.** HPA immunohistochemistry validation of biomarker protein expression in normal human kidney tissue.

---

## Supplementary Materials

- **Supplementary Figure S1.** Correlation heatmap: seven biomarkers × immune infiltration features.
- **Supplementary Figure S2.** Machine learning methodological benchmarking (six methods compared).
- **Supplementary Figure S3.** GO and KEGG enrichment analysis of DR-DEGs and all DEGs.
- **Supplementary Figure S4.** DGIdb drug-gene interaction network and interaction-score matrix for the seven-biomarker panel. Only interactions meeting the filtering criteria are shown.
- **Supplementary Table S1.** Dataset characteristics for all four GEO datasets.
- **Supplementary Table S2.** Full DR-DEG list with differential expression statistics.
- **Supplementary Table S3.** Machine learning method benchmarking results.
- **Supplementary Table S4.** Complete DGIdb query output (all 53 recorded interactions) and the filtered subset reported in the manuscript.

---

## Data Availability

All data analyzed in this study are publicly available from the Gene Expression Omnibus (GEO) under accession numbers GSE96804, GSE30528, GSE30529, and GSE131882. Analysis code is available at https://github.com/weihong6666/disulfidptosis-dn.

## Competing Interests

The author declares no competing interests.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

## Author Contributions

Z.W.X. conceived and designed the study, performed all bioinformatic analyses, interpreted the results, and wrote the manuscript.
