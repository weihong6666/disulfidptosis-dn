# Data

No raw data are stored in this repository. Everything needed to reproduce the analysis is public.

## Transcriptomic datasets

All datasets were downloaded from NCBI GEO.

| Accession | Role | Samples | Tissue | Platform |
| --------------- | ------------------ | --------------------- | ------------------- | -------------------------------- |
| GSE96804 | Training | 41 DN, 20 control | Glomerulus | Affymetrix Human Transcriptome Array 2.0 (GPL17586) |
| GSE30528 | External validation | 9 DN, 13 control | Glomerulus | Affymetrix Human Genome U133A 2.0 (GPL571) |
| GSE30529 | External validation | 10 DN, 12 control | Tubulointerstitium | Affymetrix Human Genome U133A 2.0 (GPL571) |
| GSE131882 | Single cell | 3 DN, 3 control | Kidney | 10x Genomics Chromium (20,681 cells after QC) |

### Processing notes

- **GSE96804**: probe-to-gene mapping parsed the `gene_assignment` column of the GPL17586 annotation file, taking the first annotated symbol per probe. Probes mapping to multiple genes or lacking annotation were dropped; for genes with several probes the probe with the highest mean expression was retained (31,150 genes × 61 samples). Expression was RMA-normalised and log2-transformed. Batch effects were assessed by PCA.
- **GSE30528 / GSE30529**: processed with the same annotation strategy on the GPL571 platform, with group labels taken from the series metadata.
- **GSE131882**: filtered to cells with 200–5,000 detected genes and mitochondrial content below the configured threshold; clustered and annotated against canonical kidney marker sets.

## Disulfidptosis gene set

`disulfidptosis_gene_set.txt` contains the 53 unique disulfidptosis-related genes used as the starting gene set. It was compiled from five sources:

1. the core CRISPR-screen gene set from Liu et al. 2023 (10 genes);
2. actin cytoskeleton genes implicated in disulfide crosslinking (15 genes);
3. WAVE regulatory complex and associated metabolic genes (12 genes);
4. genes nominated by the two prior DN–disulfidptosis studies (8 genes);
5. additional entries from GeneCards and FerrDB V2.

## Drug–gene interactions

Retrieved live from the DGIdb GraphQL API on 4 October 2026. Both the complete raw response (53 records) and the filtered subset reported in the manuscript are provided in `../results/tables/`.
