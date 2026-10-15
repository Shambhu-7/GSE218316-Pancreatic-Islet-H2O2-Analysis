# GSE218316: Reproducible alpha- and beta-cell analysis
# Requires: Matrix, edgeR

library(Matrix)
library(edgeR)

# 1. Import original data
counts <- readRDS("data/GSE218316/counts_annotated.rds")
metadata <- readRDS("data/GSE218316/metadata_checked.rds")
stopifnot(identical(colnames(counts), as.character(metadata$cell_id)))

# 2. Select 24-hour H2O2 and untreated alpha/beta cells
meta <- metadata[metadata$Timepoint == "24h" &
  metadata$Condition %in% c("H2O2", "Untreated") &
  metadata$Celltype %in% c("Alpha", "Beta"), , drop = FALSE]
counts_ab <- counts[, as.character(meta$cell_id), drop = FALSE]
meta$Sample <- paste(meta$Donor, meta$Celltype, meta$Condition, sep = "_")

# 3. Aggregate counts by donor, cell type and condition
sample_order <- as.vector(t(outer(
  paste0("D", 1:3),
  c("Alpha_H2O2", "Alpha_Untreated", "Beta_H2O2", "Beta_Untreated"),
  paste, sep = "_"
)))
groups <- factor(meta$Sample, levels = sample_order)
stopifnot(!anyNA(groups))
group_matrix <- Matrix::sparse.model.matrix(~ 0 + groups)
colnames(group_matrix) <- levels(groups)
pb <- as.matrix(counts_ab %*% group_matrix)
rownames(pb) <- rownames(counts_ab)

# 4. Create pseudobulk sample metadata
sample_info <- data.frame(
  Sample = colnames(pb),
  Donor = sub("_.*", "", colnames(pb)),
  Celltype = sub("^D[0-9]+_(Alpha|Beta)_.*$", "\\1", colnames(pb)),
  Condition = sub("^D[0-9]+_(Alpha|Beta)_", "", colnames(pb)),
  stringsAsFactors = FALSE
)

# 5. Donor-adjusted robust edgeR analysis
run_edger <- function(celltype) {
  idx <- sample_info$Celltype == celltype
  md <- sample_info[idx, , drop = FALSE]
  md$Donor <- factor(md$Donor)
  md$Condition <- factor(md$Condition, levels = c("Untreated", "H2O2"))
  design <- model.matrix(~ Donor + Condition, data = md)
  y <- DGEList(counts = pb[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "ConditionH2O2")
  topTags(test, n = Inf)$table
}

alpha_results <- run_edger("Alpha")
beta_results <- run_edger("Beta")
shared_genes <- intersect(
  rownames(alpha_results)[alpha_results$FDR < 0.05],
  rownames(beta_results)[beta_results$FDR < 0.05]
)

# 6. Report and verify expected results
cat("Alpha DEGs:", sum(alpha_results$FDR < 0.05), "\n")
cat("Beta DEGs:", sum(beta_results$FDR < 0.05), "\n")
cat("Shared DEGs:", length(shared_genes), "\n")
stopifnot(
  ncol(pb) == 12,
  nrow(alpha_results) == 10036,
  nrow(beta_results) == 8493,
  sum(alpha_results$FDR < 0.05) == 71,
  sum(beta_results$FDR < 0.05) == 29,
  length(shared_genes) == 9
)

# 7. Combined interaction analysis
# Original manuscript model: donor + cell type * condition

joint_info <- sample_info
joint_info$Donor <- factor(joint_info$Donor)
joint_info$Celltype <- factor(joint_info$Celltype,
  levels = c("Alpha", "Beta"))
joint_info$Condition <- factor(joint_info$Condition,
  levels = c("Untreated", "H2O2"))

design_joint <- model.matrix(
  ~ Donor + Celltype * Condition,
  data = joint_info
)

y_joint <- DGEList(counts = pb)
keep_joint <- filterByExpr(y_joint, design = design_joint)
y_joint <- y_joint[keep_joint, , keep.lib.sizes = FALSE]
y_joint <- calcNormFactors(y_joint, method = "TMM")
y_joint <- estimateDisp(y_joint, design = design_joint)

fit_joint <- glmQLFit(
  y_joint, design = design_joint, robust = TRUE
)

test_joint <- glmQLFTest(
  fit_joint,
  coef = "CelltypeBeta:ConditionH2O2"
)

interaction_results <- topTags(test_joint, n = Inf)$table

cat("Interaction genes tested:",
    nrow(interaction_results), "\n")
cat("Significant interactions:",
    sum(interaction_results$FDR < 0.05), "\n")
cat("Minimum interaction FDR:",
    min(interaction_results$FDR), "\n")

stopifnot(
  nrow(interaction_results) == 10221,
  sum(interaction_results$FDR < 0.05) == 0,
  isTRUE(all.equal(
    min(interaction_results$FDR),
    0.5108125, tolerance = 1e-6
  ))
)

# ============================================================
# HALLMARK GSEA: ALPHA AND BETA CELLS

# Load gene annotation from repository inputs
annotation_input <- readRDS(
  "data/GSE218316/gene_annotation.rds"
)

# ============================================================

library(fgsea)

# Load frozen Hallmark gene sets
gsea_reference <- readRDS('GSE218316_Hallmark_GSEA_Verified.rds')
hallmark_pathways <- gsea_reference$hallmark_pathways

# Map Ensembl identifiers to gene symbols
make_gsea_ranks <- function(results, annotation) {
  df <- data.frame(
    feature_id = rownames(results),
    results,
    row.names = NULL
  )
  df <- merge(
    df,
    annotation[, c('feature_id', 'gene_symbol')],
    by = 'feature_id'
  )
  df <- df[!is.na(df$gene_symbol) & nzchar(df$gene_symbol), ]
  df$rank_score <- sign(df$logFC) * sqrt(df$F)
  df <- df[order(-df$rank_score), ]
  df <- df[!duplicated(df$gene_symbol), ]
  ranks <- df$rank_score
  names(ranks) <- df$gene_symbol
  ranks
}

alpha_ranked <- make_gsea_ranks(alpha_results, annotation_input)
beta_ranked <- make_gsea_ranks(beta_results, annotation_input)

# Verify the regenerated rankings
stopifnot(
  isTRUE(all.equal(alpha_ranked, gsea_reference$alpha_ranked)),
  isTRUE(all.equal(beta_ranked, gsea_reference$beta_ranked))
)

# Run enrichment analysis
set.seed(123)
alpha_gsea <- fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = alpha_ranked,
  minSize = 15,
  maxSize = 500
)

set.seed(123)
beta_gsea <- fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = beta_ranked,
  minSize = 15,
  maxSize = 500
)

# Verify enrichment results by pathway
check_gsea <- function(observed, expected) {
  idx <- match(expected$pathway, observed$pathway)
  stopifnot(
    !anyNA(idx),
    length(idx) == nrow(observed),
    isTRUE(all.equal(expected$NES, observed$NES[idx], tolerance = 1e-8)),
    isTRUE(all.equal(expected$padj, observed$padj[idx], tolerance = 1e-8)),
    identical(as.integer(expected$size), as.integer(observed$size[idx]))
  )
}

check_gsea(alpha_gsea, gsea_reference$alpha_gsea)
check_gsea(beta_gsea, gsea_reference$beta_gsea)

cat('Hallmark GSEA successfully reproduced for both cell types.\n')
