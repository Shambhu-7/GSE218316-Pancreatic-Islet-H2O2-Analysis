# Shared transcriptional responses to hydrogen peroxide in human pancreatic alpha and beta cells

This repository contains the reproducible analysis and figure-generation
code associated with the manuscript:

**Shared transcriptional responses to hydrogen peroxide in human
pancreatic alpha and beta cells**

The analysis uses publicly available single-cell RNA-sequencing data
from **GEO accession GSE218316** to examine transcriptional responses of
human pancreatic alpha and beta cells following 24-hour hydrogen
peroxide (H₂O₂) exposure.

## Analysis overview

The analysis includes:

-   1,024 pancreatic alpha cells
-   594 pancreatic beta cells
-   Three independent human donors
-   H₂O₂-treated and untreated conditions
-   Donor-level pseudobulk differential-expression analysis using edgeR
-   Separate alpha- and beta-cell differential-expression analyses
-   Joint donor-aware cell-type × condition interaction analysis
-   Hallmark gene-set enrichment analysis
-   Correlation-aware pathway interaction analysis
-   Exploratory global-expression heatmap and single-cell UMAP
    visualization

The separate pseudobulk analyses identified:

-   71 differentially expressed genes in alpha cells
-   29 differentially expressed genes in beta cells
-   9 genes shared between the two cell types

The joint interaction model tested 10,221 genes and detected no
significant cell-type × H₂O₂ interactions at FDR \< 0.05.

## Repository contents

### Analysis

`GSE218316_Reproducible_Analysis.R`

Main reproducible analysis script containing the pseudobulk
differential-expression, interaction, and pathway analyses.

### Figure generation

`GSE218316_Final_Figures.R`

Reproducible script used to generate the final manuscript figures:

-   Figure 1 --- Differential-expression volcano plots
-   Figure 2 --- Shared H₂O₂-responsive gene expression
-   Figure 3 --- Hallmark transcriptional signatures
-   Figure S1 --- DEG overlap
-   Figure S2 --- Pseudobulk multidimensional scaling
-   Figure S3 --- Pseudobulk count depth
-   Figure S4 --- Average pseudobulk counts per cell
-   Figure S5 --- Global expression heatmap of the 1,000 most variable
    genes
-   Figure S6 --- Single-cell UMAP

### Processed data

The processed inputs required to reproduce the analyses are included
under:

`data/GSE218316/`

Files:

-   `counts_annotated.rds`
-   `metadata_checked.rds`
-   `gene_annotation.rds`

### Analysis results

`results/GSE218316_Complete_Reproducible_Results.rds`

contains the saved results from the reproducible analysis.

`GSE218316_Hallmark_GSEA_Verified.rds`

contains the verified Hallmark GSEA results used for figure generation.

## Statistical analysis

Differential expression was performed using donor-level pseudobulk
counts and the edgeR quasi-likelihood framework.

Alpha and beta cells were first analyzed separately using donor-adjusted
models.

Cell-type differences in the transcriptional response to H₂O₂ were
evaluated using a joint model containing donor-by-cell-type and
condition-by-cell-type terms. Statistical inference therefore accounts
for the donor structure of the experiment rather than treating
individual cells as independent biological replicates.

Hallmark pathway enrichment was evaluated using ranked gene-level
statistics. Correlation-aware pathway testing was additionally used when
assessing evidence for differential pathway responses between cell
types.

## Reproducibility

The analysis and figure-generation workflows were validated against the
final manuscript results.

Expected analysis checkpoints include:

-   Alpha cells: 1,024
-   Beta cells: 594
-   Pseudobulk samples: 12
-   Alpha genes tested: 10,036
-   Beta genes tested: 8,493
-   Alpha DEGs (FDR \< 0.05): 71
-   Beta DEGs (FDR \< 0.05): 29
-   Shared DEGs: 9
-   Genes tested in the joint interaction analysis: 10,221
-   Significant cell-type × condition interactions (FDR \< 0.05): 0
-   Minimum interaction FDR: 0.5108

The nine shared H₂O₂-responsive genes are:

`TGM2`, `CRYBA2`, `PLOD2`, `ACTG1`, `CITED2`, `NFKBIA`, `SHISAL2B`,
`S100A10`, and `LOXL4`.


### Within-cell-type CAMERA sensitivity analysis

The final manuscript includes correlation-aware Hallmark
sensitivity analyses performed separately in alpha and beta
cells using donor-adjusted pseudobulk models.

Script:

`scripts/sensitivity/GSE218316_Hallmark_CAMERA_Sensitivity_v4.R`

Run from the repository root:

```bash
Rscript scripts/sensitivity/GSE218316_Hallmark_CAMERA_Sensitivity_v4.R
```

Required R packages: Matrix, edgeR, limma, and msigdbr.

The script reads the three processed RDS inputs in
`data/GSE218316/` and
`GSE218316_Hallmark_GSEA_Verified.rds`.

The manuscript-matched supplementary results are stored in
`results/supplementary/`:

- `Table_S_Hallmark_fgsea_ALL_Combined.csv` — Table S1
- `Table_S_CAMERA_ALL_Hallmark.csv` — Table S2
- `Table_S_FGSEA_vs_CAMERA_Selected9.csv` — Table S3

The CAMERA script regenerates these CSV results in `results/`.
The supplied reference copies are retained in
`results/supplementary/` for comparison.

The R environment record is available at
`reproducibility/CAMERA_Sensitivity_sessionInfo.txt`.

Analysis environment: R 4.5.2, edgeR 4.8.2,
limma 3.66.0, msigdbr 25.1.1,
MSigDB 2025.1.Hs.

These sensitivity analyses are not independent biological
validation.

## Running the analysis

Run the scripts from the repository root.

Main analysis:

`Rscript GSE218316_Reproducible_Analysis.R`

Figure generation:

`Rscript GSE218316_Final_Figures.R`

The figure-generation script writes the final PNG figures to the
`figures/` directory.

## Data availability

The original dataset is publicly available through the NCBI Gene
Expression Omnibus under accession **GSE218316**.

## Citation

The manuscript citation and permanent Zenodo DOI will be added following
repository archiving.
