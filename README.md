# Shared Hydrogen Peroxide Responses in Human Pancreatic Alpha and Beta Cells

## Overview

This repository contains the computational analysis supporting
a secondary analysis of public single-cell RNA-sequencing data
from human pancreatic alpha and beta cells exposed to hydrogen
peroxide (H2O2).

## Data source

GEO accession: GSE218316
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218316

The analysis includes 1,024 alpha cells and 594 beta cells
from three donors, comparing 24-hour H2O2 exposure with
untreated controls.

## Analysis methods

- Donor-level pseudobulk aggregation (12 samples)
- Differential expression using edgeR quasi-likelihood models
- Donor-adjusted models: ~ Donor + Condition
- Interaction model: ~ Donor + Celltype * Condition
- MSigDB Hallmark enrichment using fgseaMultilevel
- Gene ranking: sign(logFC) * sqrt(F)
- Benjamini-Hochberg FDR correction

## Principal findings

- Alpha cells: 71 differentially expressed genes (FDR < 0.05)
- Beta cells: 29 differentially expressed genes (FDR < 0.05)
- Nine shared differentially expressed genes
- No significant cell-type-by-treatment interactions
- Oxidative phosphorylation: alpha NES = -1.433, FDR = 0.0170
- Oxidative phosphorylation: beta NES = +1.876, FDR = 0.0000156

Opposite enrichment directions do not establish a statistically
significant difference between cell types.

## Repository files

- GSE218316_Reproducible_Analysis.R: analysis script
- GSE218316_Hallmark_GSEA_Verified.rds: frozen gene sets and GSEA reference
- results/GSE218316_Complete_Reproducible_Results.rds: validated results

## Software

- R
- edgeR 4.8.2
- fgsea 1.36.2

## Reproducibility

The analysis script has been validated against the original
differential-expression and Hallmark GSEA results.

The script currently requires three processed input files:

- data/GSE218316/counts_annotated.rds
- data/GSE218316/metadata_checked.rds
- data/GSE218316/gene_annotation.rds

These input files are not yet distributed in this repository.
A tested data-access or reconstruction workflow will be added
before the first archived release.

## Citation

Manuscript citation and Zenodo DOI will be added after publication
and repository archiving.
