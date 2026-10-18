# ======================================================================
# GSE218316_Final_Figures.R
#
# Final, standalone, reproducible PNG figure script for:
#   Figure 1
#   Figure 2
#   Figure 3
#   Figure S1
#   Figure S2
#   Figure S3
#   Figure S4
#   Figure S5
#   Figure S6
#
# Study:
#   Shared transcriptional responses to hydrogen peroxide in
#   human pancreatic alpha and beta cells
#
# Run from the repository root.
#
# Required inputs:
#   data/GSE218316/counts_annotated.rds
#   data/GSE218316/metadata_checked.rds
#   data/GSE218316/gene_annotation.rds
#   GSE218316_Hallmark_GSEA_Verified.rds
#
# Statistical inference:
#   Alpha and Beta differential expression:
#       ~ Donor + Condition
#
# Figures S5 and S6 are exploratory global-expression visualizations.
# Differential-expression inference remains donor-level pseudobulk.
#
# PNG only; no PDF output.
# ======================================================================


# ======================================================================
# 0. REQUIRED PACKAGES
# ======================================================================

required_packages <- c(
  "Matrix",
  "edgeR",
  "ggplot2",
  "scales",
  "pheatmap",
  "ggrepel",
  "uwot",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Install these packages before running the script: ",
      paste(missing_packages, collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(scales)
  library(pheatmap)
  library(ggrepel)
  library(uwot)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

set.seed(218316)


# ======================================================================
# 1. LOAD AND VERIFY INPUTS
# ======================================================================

counts <- readRDS(
  "data/GSE218316/counts_annotated.rds"
)

metadata <- readRDS(
  "data/GSE218316/metadata_checked.rds"
)

gene_annotation <- readRDS(
  "data/GSE218316/gene_annotation.rds"
)

gsea_reference <- readRDS(
  "GSE218316_Hallmark_GSEA_Verified.rds"
)

stopifnot(
  ncol(counts) == nrow(metadata),
  identical(
    colnames(counts),
    as.character(metadata$cell_id)
  ),
  !anyDuplicated(metadata$cell_id),
  all(
    c(
      "feature_id",
      "gene_symbol"
    ) %in% colnames(gene_annotation)
  ),
  all(
    c(
      "alpha_gsea",
      "beta_gsea"
    ) %in% names(gsea_reference)
  )
)

gene_annotation$feature_id <- as.character(
  gene_annotation$feature_id
)

gene_annotation$gene_symbol <- as.character(
  gene_annotation$gene_symbol
)


# ======================================================================
# 2. EXACT MANUSCRIPT POPULATION
#    24 h H2O2 vs untreated
#    Alpha/Beta
#    Donors D1-D3
# ======================================================================

keep_cells <- (
  metadata$Timepoint == "24h" &
    metadata$Condition %in% c(
      "H2O2",
      "Untreated"
    ) &
    metadata$Celltype %in% c(
      "Alpha",
      "Beta"
    ) &
    metadata$Donor %in% c(
      "D1",
      "D2",
      "D3"
    )
)

meta <- metadata[
  keep_cells,
  ,
  drop = FALSE
]

counts_ab <- counts[
  ,
  as.character(meta$cell_id),
  drop = FALSE
]

stopifnot(
  nrow(meta) == 1618,
  ncol(counts_ab) == 1618,
  sum(meta$Celltype == "Alpha") == 1024,
  sum(meta$Celltype == "Beta") == 594,
  sum(meta$Condition == "H2O2") == 790,
  sum(meta$Condition == "Untreated") == 828
)

meta$Sample <- paste(
  meta$Donor,
  meta$Celltype,
  meta$Condition,
  sep = "_"
)


# ======================================================================
# 3. CONSTRUCT 12 DONOR-LEVEL PSEUDOBULK SAMPLES
# ======================================================================

sample_order <- as.vector(
  t(
    outer(
      paste0(
        "D",
        1:3
      ),
      c(
        "Alpha_H2O2",
        "Alpha_Untreated",
        "Beta_H2O2",
        "Beta_Untreated"
      ),
      paste,
      sep = "_"
    )
  )
)

groups <- factor(
  meta$Sample,
  levels = sample_order
)

stopifnot(
  !anyNA(groups)
)

group_matrix <- Matrix::sparse.model.matrix(
  ~ 0 + groups
)

colnames(group_matrix) <- levels(groups)

pb <- as.matrix(
  counts_ab %*% group_matrix
)

rownames(pb) <- rownames(counts_ab)

stopifnot(
  ncol(pb) == 12,
  identical(
    colnames(pb),
    sample_order
  )
)


# ======================================================================
# 4. PSEUDOBULK SAMPLE METADATA
# ======================================================================

sample_info <- data.frame(
  Sample = colnames(pb),
  Donor = sub(
    "_.*",
    "",
    colnames(pb)
  ),
  Celltype = sub(
    "^D[0-9]+_(Alpha|Beta)_.*$",
    "\\1",
    colnames(pb)
  ),
  Condition = sub(
    "^D[0-9]+_(Alpha|Beta)_",
    "",
    colnames(pb)
  ),
  stringsAsFactors = FALSE
)

sample_info$Donor <- factor(
  sample_info$Donor,
  levels = c(
    "D1",
    "D2",
    "D3"
  )
)

sample_info$Celltype <- factor(
  sample_info$Celltype,
  levels = c(
    "Alpha",
    "Beta"
  )
)

sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c(
    "Untreated",
    "H2O2"
  )
)

cell_counts <- table(groups)

sample_info$n_cells <- as.integer(
  cell_counts[
    sample_info$Sample
  ]
)

sample_info$library_size <- as.numeric(
  colSums(pb)
)

sample_info$library_millions <-
  sample_info$library_size / 1e6

sample_info$counts_per_cell <-
  sample_info$library_size /
  sample_info$n_cells


# ======================================================================
# 5. DONOR-ADJUSTED edgeR ANALYSES
#    Model: ~ Donor + Condition
# ======================================================================

run_edger <- function(celltype) {

  idx <- sample_info$Celltype == celltype

  md <- droplevels(
    sample_info[
      idx,
      ,
      drop = FALSE
    ]
  )

  md$Donor <- factor(
    md$Donor,
    levels = c(
      "D1",
      "D2",
      "D3"
    )
  )

  md$Condition <- factor(
    md$Condition,
    levels = c(
      "Untreated",
      "H2O2"
    )
  )

  design <- model.matrix(
    ~ Donor + Condition,
    data = md
  )

  y <- DGEList(
    counts = pb[
      ,
      idx,
      drop = FALSE
    ]
  )

  keep <- filterByExpr(
    y,
    design = design
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- calcNormFactors(
    y,
    method = "TMM"
  )

  y <- estimateDisp(
    y,
    design = design
  )

  fit <- glmQLFit(
    y,
    design = design,
    robust = TRUE
  )

  test <- glmQLFTest(
    fit,
    coef = "ConditionH2O2"
  )

  results <- topTags(
    test,
    n = Inf
  )$table

  list(
    results = results,
    y = y,
    metadata = md
  )
}

alpha_fit <- run_edger(
  "Alpha"
)

beta_fit <- run_edger(
  "Beta"
)

alpha_results <- alpha_fit$results
beta_results <- beta_fit$results

alpha_sig <- rownames(
  alpha_results
)[
  alpha_results$FDR < 0.05
]

beta_sig <- rownames(
  beta_results
)[
  beta_results$FDR < 0.05
]

shared_ids <- intersect(
  alpha_sig,
  beta_sig
)

stopifnot(
  nrow(alpha_results) == 10036,
  nrow(beta_results) == 8493,
  length(alpha_sig) == 71,
  length(beta_sig) == 29,
  length(shared_ids) == 9
)


# ======================================================================
# 6. GENE ANNOTATION HELPERS
# ======================================================================

symbol_from_id <- function(ids) {

  ans <- gene_annotation$gene_symbol[
    match(
      ids,
      gene_annotation$feature_id
    )
  ]

  bad <- is.na(ans) |
    ans == ""

  ans[bad] <- ids[bad]

  ans
}

annotate_results <- function(x) {

  out <- data.frame(
    feature_id = rownames(x),
    x,
    row.names = NULL,
    check.names = FALSE
  )

  out$gene_symbol <- symbol_from_id(
    out$feature_id
  )

  out
}

alpha_table <- annotate_results(
  alpha_results
)

beta_table <- annotate_results(
  beta_results
)

expected_shared <- c(
  "TGM2",
  "CRYBA2",
  "PLOD2",
  "ACTG1",
  "CITED2",
  "NFKBIA",
  "SHISAL2B",
  "S100A10",
  "LOXL4"
)

stopifnot(
  setequal(
    symbol_from_id(shared_ids),
    expected_shared
  )
)

cat(
  "\n========================================\n"
)

cat(
  "VERIFIED BEFORE FIGURE GENERATION\n"
)

cat(
  "========================================\n"
)

cat(
  "Alpha cells:",
  sum(meta$Celltype == "Alpha"),
  "\n"
)

cat(
  "Beta cells:",
  sum(meta$Celltype == "Beta"),
  "\n"
)

cat(
  "Pseudobulk samples:",
  ncol(pb),
  "\n"
)

cat(
  "Alpha genes tested:",
  nrow(alpha_results),
  "\n"
)

cat(
  "Beta genes tested:",
  nrow(beta_results),
  "\n"
)

cat(
  "Alpha DEGs:",
  length(alpha_sig),
  "\n"
)

cat(
  "Beta DEGs:",
  length(beta_sig),
  "\n"
)

cat(
  "Shared DEGs:",
  length(shared_ids),
  "\n\n"
)


# ======================================================================
# 7. COMMON TMM-NORMALIZED PSEUDOBULK MATRIX
# ======================================================================

y_all <- DGEList(
  counts = pb
)

y_all <- calcNormFactors(
  y_all,
  method = "TMM"
)

logcpm_all <- cpm(
  y_all,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)


# ======================================================================
# FIGURE 1
# SHARED H2O2-RESPONSIVE GENES ACROSS BIOLOGICAL DONORS
#
# 9 shared genes
# donor-paired delta log2CPM = H2O2 - untreated
# no row scaling
# hierarchical clustering of rows and columns
# ======================================================================

shared_order <- expected_shared

shared_feature_ids <- gene_annotation$feature_id[
  match(
    shared_order,
    gene_annotation$gene_symbol
  )
]

stopifnot(
  length(shared_feature_ids) == 9,
  !anyNA(shared_feature_ids),
  all(
    shared_feature_ids %in%
      rownames(logcpm_all)
  ),
  setequal(
    shared_feature_ids,
    shared_ids
  )
)

heatmap_matrix <- matrix(
  NA_real_,
  nrow = length(shared_order),
  ncol = 6,
  dimnames = list(
    shared_order,
    c(
      "α-D1",
      "α-D2",
      "α-D3",
      "β-D1",
      "β-D2",
      "β-D3"
    )
  )
)

column_index <- 1L

for (ct in c("Alpha", "Beta")) {

  for (d in c("D1", "D2", "D3")) {

    treated_sample <- paste(
      d,
      ct,
      "H2O2",
      sep = "_"
    )

    untreated_sample <- paste(
      d,
      ct,
      "Untreated",
      sep = "_"
    )

    heatmap_matrix[
      ,
      column_index
    ] <- (
      logcpm_all[
        shared_feature_ids,
        treated_sample
      ] -
        logcpm_all[
          shared_feature_ids,
          untreated_sample
        ]
    )

    column_index <- column_index + 1L
  }
}

annotation_col_fig1 <- data.frame(
  Cell_type = factor(
    c(
      rep("Alpha", 3),
      rep("Beta", 3)
    ),
    levels = c(
      "Alpha",
      "Beta"
    )
  ),
  Donor = factor(
    rep(
      c(
        "D1",
        "D2",
        "D3"
      ),
      2
    ),
    levels = c(
      "D1",
      "D2",
      "D3"
    )
  )
)

rownames(annotation_col_fig1) <-
  colnames(heatmap_matrix)

annotation_colors_fig1 <- list(
  Cell_type = c(
    Alpha = "#2878B5",
    Beta = "#D55E00"
  ),
  Donor = c(
    D1 = "#4D4D4D",
    D2 = "#969696",
    D3 = "#D0D0D0"
  )
)

max_abs <- max(
  abs(heatmap_matrix),
  na.rm = TRUE
)

fig1_breaks <- seq(
  -max_abs,
  max_abs,
  length.out = 101
)

fig1_colors <- colorRampPalette(
  c(
    "#2878B5",
    "#A9CFE5",
    "#FFF7D6",
    "#F3C567",
    "#D55E00"
  )
)(100)

pheatmap::pheatmap(
  heatmap_matrix,
  scale = "none",
  color = fig1_colors,
  breaks = fig1_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col_fig1,
  annotation_colors = annotation_colors_fig1,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  border_color = "white",
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellwidth = 48,
  cellheight = 38,
  fontsize = 11,
  fontsize_row = 12,
  fontsize_col = 11,
  angle_col = 45,
  treeheight_row = 45,
  treeheight_col = 45,
  legend_breaks = pretty(
    c(
      -max_abs,
      max_abs
    ),
    n = 5
  ),
  main = "Shared H₂O₂-responsive genes",
  filename =
    "figures/Figure_1_Shared_H2O2_gene_expression.png",
  width = 9,
  height = 7.5
)


# ======================================================================
# FIGURE 2
# DIFFERENTIAL-EXPRESSION VOLCANO PLOTS
#
# Significance definition:
#   FDR < 0.05
#
# No absolute fold-change threshold.
# ======================================================================

make_volcano <- function(tab, celltype) {

  x <- tab

  x$Celltype <- celltype

  x$neglog10FDR <- -log10(
    pmax(
      x$FDR,
      .Machine$double.xmin
    )
  )

  x$Plot_group <- "Not significant"

  x$Plot_group[
    x$FDR < 0.05 &
      x$logFC > 0
  ] <- "Upregulated"

  x$Plot_group[
    x$FDR < 0.05 &
      x$logFC < 0
  ] <- "Downregulated"

  x
}

volcano_df <- rbind(
  make_volcano(
    alpha_table,
    "Alpha"
  ),
  make_volcano(
    beta_table,
    "Beta"
  )
)

volcano_df$Celltype <- factor(
  volcano_df$Celltype,
  levels = c(
    "Alpha",
    "Beta"
  )
)

volcano_df$Plot_group <- factor(
  volcano_df$Plot_group,
  levels = c(
    "Upregulated",
    "Downregulated",
    "Not significant"
  )
)

alpha_label_symbols <- c(
  "TGM2",
  "MSMB",
  "PLCE1",
  "LDHA",
  "NDRG1"
)

beta_label_symbols <- c(
  "C1QL1",
  "DEPP1",
  "IGFBP3",
  "SOX9",
  "MT1M"
)

alpha_label_index <- match(
  alpha_label_symbols,
  alpha_table$gene_symbol
)

beta_label_index <- match(
  beta_label_symbols,
  beta_table$gene_symbol
)

stopifnot(
  !anyNA(alpha_label_index),
  !anyNA(beta_label_index)
)

alpha_labels <- alpha_table[
  alpha_label_index,
  ,
  drop = FALSE
]

alpha_labels$Celltype <- "Alpha"

beta_labels <- beta_table[
  beta_label_index,
  ,
  drop = FALSE
]

beta_labels$Celltype <- "Beta"

label_df <- rbind(
  alpha_labels,
  beta_labels
)

label_df$neglog10FDR <- -log10(
  pmax(
    label_df$FDR,
    .Machine$double.xmin
  )
)

label_df$Celltype <- factor(
  label_df$Celltype,
  levels = c(
    "Alpha",
    "Beta"
  )
)

label_df$Plot_group <- "Not significant"

label_df$Plot_group[
  label_df$FDR < 0.05 &
    label_df$logFC > 0
] <- "Upregulated"

label_df$Plot_group[
  label_df$FDR < 0.05 &
    label_df$logFC < 0
] <- "Downregulated"

label_df$Plot_group <- factor(
  label_df$Plot_group,
  levels = c(
    "Upregulated",
    "Downregulated",
    "Not significant"
  )
)

p2 <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = neglog10FDR,
    colour = Plot_group
  )
) +
  geom_point(
    size = 2.3,
    alpha = 0.85
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 1.0,
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.8,
    colour = "black"
  ) +
  facet_wrap(
    ~ Celltype,
    nrow = 1,
    labeller = as_labeller(
      c(
        Alpha = "Alpha cells",
        Beta = "Beta cells"
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "Upregulated" = "#4FA392",
      "Downregulated" = "#C6A06A",
      "Not significant" = "#D0D2D1"
    )
  ) +
  coord_cartesian(
    xlim = c(
      -4,
      4
    )
  ) +
  labs(
    title = expression(
      bold(
        "Volcano: H"[2] *
          "O"[2] *
          " vs untreated"
      )
    ),
    subtitle =
      "(+ve) log2FC = higher with H2O2; (-ve) log2FC = higher in untreated",
    x = expression(
      log[2] ~
        "Fold Change (H"[2] *
        "O"[2] *
        " / untreated)"
    ),
    y = expression(
      -log[10] ~ "(FDR)"
    ),
    colour = NULL
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 19,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 13,
      hjust = 0,
      margin = margin(
        b = 12
      )
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 16
    ),
    axis.title = element_text(
      size = 15
    ),
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 1
    ),
    axis.ticks = element_line(
      linewidth = 0.8
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 12
    )
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(
      x = logFC,
      y = neglog10FDR,
      label = gene_symbol,
      colour = Plot_group
    ),
    inherit.aes = FALSE,
    size = 4.2,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  )

ggsave(
  filename =
    "figures/Figure_2_Differential_expression_volcano_plots.png",
  plot = p2,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE 3
# SELECTED HALLMARK GSEA SIGNATURES
#
# Nine biologically selected Hallmark sets.
# Bubble size = -log10(FDR)
# Bubble colour = NES
# Filled = FDR < 0.05
# Hollow = FDR >= 0.05
# ======================================================================

alpha_gsea <- as.data.frame(
  gsea_reference$alpha_gsea
)

beta_gsea <- as.data.frame(
  gsea_reference$beta_gsea
)

a_ox <- alpha_gsea[
  alpha_gsea$pathway ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  ,
  drop = FALSE
]

b_ox <- beta_gsea[
  beta_gsea$pathway ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  ,
  drop = FALSE
]

stopifnot(
  nrow(a_ox) == 1,
  nrow(b_ox) == 1,
  isTRUE(
    all.equal(
      a_ox$NES,
      -1.4327130,
      tolerance = 1e-6
    )
  ),
  isTRUE(
    all.equal(
      a_ox$padj,
      0.01700407,
      tolerance = 1e-6
    )
  ),
  isTRUE(
    all.equal(
      b_ox$NES,
      1.876257,
      tolerance = 1e-6
    )
  ),
  isTRUE(
    all.equal(
      b_ox$padj,
      1.558519e-05,
      tolerance = 1e-6
    )
  )
)

selected_pathways <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

pretty_name <- c(
  "HALLMARK_HYPOXIA" =
    "Hypoxia",
  "HALLMARK_GLYCOLYSIS" =
    "Glycolysis",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" =
    "TNF-alpha / NF-kappaB",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" =
    "Reactive oxygen species",
  "HALLMARK_MTORC1_SIGNALING" =
    "mTORC1 signaling",
  "HALLMARK_P53_PATHWAY" =
    "p53 pathway",
  "HALLMARK_APOPTOSIS" =
    "Apoptosis",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" =
    "Unfolded protein response",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" =
    "Oxidative phosphorylation"
)

alpha_pathway_df <- alpha_gsea[
  alpha_gsea$pathway %in%
    selected_pathways,
  ,
  drop = FALSE
]

beta_pathway_df <- beta_gsea[
  beta_gsea$pathway %in%
    selected_pathways,
  ,
  drop = FALSE
]

stopifnot(
  nrow(alpha_pathway_df) == 9,
  nrow(beta_pathway_df) == 9
)

alpha_pathway_df$Celltype <- "Alpha"
beta_pathway_df$Celltype <- "Beta"

pathway_df <- rbind(
  alpha_pathway_df,
  beta_pathway_df
)

pathway_df$Pathway <- unname(
  pretty_name[
    pathway_df$pathway
  ]
)

pathway_df$Pathway <- factor(
  pathway_df$Pathway,
  levels = rev(
    unname(
      pretty_name[
        selected_pathways
      ]
    )
  )
)

pathway_df$Celltype <- factor(
  pathway_df$Celltype,
  levels = c(
    "Alpha",
    "Beta"
  )
)

pathway_df$Significance <- ifelse(
  pathway_df$padj < 0.05,
  "FDR < 0.05",
  "FDR ≥ 0.05"
)

pathway_df$Significance <- factor(
  pathway_df$Significance,
  levels = c(
    "FDR < 0.05",
    "FDR ≥ 0.05"
  )
)

pathway_df$minus_log10_FDR <- pmin(
  -log10(
    pmax(
      pathway_df$padj,
      .Machine$double.xmin
    )
  ),
  20
)

p3 <- ggplot(
  pathway_df,
  aes(
    x = NES,
    y = Pathway,
    colour = NES,
    size = minus_log10_FDR,
    shape = Significance
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey55"
  ) +
  geom_point(
    stroke = 1.2
  ) +
  facet_wrap(
    ~ Celltype,
    nrow = 1,
    labeller = as_labeller(
      c(
        Alpha = "Alpha cells",
        Beta = "Beta cells"
      )
    )
  ) +
  scale_colour_gradient(
    low = "#173A56",
    high = "#4DA8F0",
    name = "NES"
  ) +
  scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "FDR ≥ 0.05" = 1
    ),
    name = "Significance"
  ) +
  scale_size_continuous(
    name = expression(
      -log[10] * "(FDR)"
    ),
    range = c(
      2.5,
      9
    ),
    limits = c(
      0,
      20
    )
  ) +
  labs(
    title =
      "Hallmark GSEA: H₂O₂ vs untreated",
    x =
      "Normalized enrichment score (NES)",
    y = NULL
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.y = element_text(
      size = 11
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

ggsave(
  filename =
    "figures/Figure_3_Hallmark_transcriptional_signatures.png",
  plot = p3,
  width = 11,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE S1
# BORDERLESS VENN DIAGRAM OF DEG OVERLAP
# ======================================================================

alpha_total <- length(
  alpha_sig
)

beta_total <- length(
  beta_sig
)

shared_total <- length(
  shared_ids
)

alpha_only <- length(
  setdiff(
    alpha_sig,
    beta_sig
  )
)

beta_only <- length(
  setdiff(
    beta_sig,
    alpha_sig
  )
)

stopifnot(
  alpha_total == 71,
  beta_total == 29,
  shared_total == 9,
  alpha_only == 62,
  beta_only == 20
)

make_circle <- function(
  cx,
  cy,
  radius,
  n = 500
) {

  theta <- seq(
    0,
    2 * pi,
    length.out = n
  )

  data.frame(
    x = cx +
      radius *
      cos(theta),
    y = cy +
      radius *
      sin(theta)
  )
}

alpha_circle <- make_circle(
  cx = 0,
  cy = 0,
  radius = 1.55
)

beta_circle <- make_circle(
  cx = 1.75,
  cy = 0,
  radius = 1.55
)

pS1 <- ggplot() +
  geom_polygon(
    data = alpha_circle,
    aes(
      x = x,
      y = y
    ),
    fill = "#4FA392",
    alpha = 0.55,
    colour = NA
  ) +
  geom_polygon(
    data = beta_circle,
    aes(
      x = x,
      y = y
    ),
    fill = "#F4A261",
    alpha = 0.55,
    colour = NA
  ) +
  annotate(
    "text",
    x = -0.55,
    y = 0,
    label = alpha_only,
    size = 8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.875,
    y = 0,
    label = shared_total,
    size = 8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.30,
    y = 0,
    label = beta_only,
    size = 8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = -0.45,
    y = 1.82,
    label = "Alpha cells",
    size = 5.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.20,
    y = 1.82,
    label = "Beta cells",
    size = 5.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.875,
    y = 0.42,
    label = "Shared",
    size = 4.4
  ) +
  annotate(
    "text",
    x = -0.45,
    y = -1.82,
    label = "71 DEGs",
    size = 4.8
  ) +
  annotate(
    "text",
    x = 2.20,
    y = -1.82,
    label = "29 DEGs",
    size = 4.8
  ) +
  coord_equal(
    xlim = c(
      -2.0,
      3.75
    ),
    ylim = c(
      -2.25,
      2.25
    ),
    expand = FALSE
  ) +
  labs(
    title =
      "Overlap of H₂O₂-responsive genes",
    subtitle =
      "Differentially expressed genes (FDR < 0.05)"
  ) +
  theme_void(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 12.5,
      hjust = 0.5
    ),
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

ggsave(
  filename =
    "figures/Figure_S1_DEG_overlap.png",
  plot = pS1,
  width = 7.5,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE S2
# PSEUDOBULK MULTIDIMENSIONAL SCALING
# ======================================================================

mds_result <- plotMDS(
  y_all,
  plot = FALSE
)

mds_data <- sample_info

mds_data$MDS1 <-
  mds_result$x

mds_data$MDS2 <-
  mds_result$y

stopifnot(
  nrow(mds_data) == 12
)

pS2 <- ggplot(
  mds_data,
  aes(
    x = MDS1,
    y = MDS2,
    colour = Celltype,
    shape = Condition
  )
) +
  geom_line(
    aes(
      group = interaction(
        Donor,
        Celltype
      )
    ),
    linewidth = 0.9,
    alpha = 0.55
  ) +
  geom_point(
    size = 5
  ) +
  facet_wrap(
    ~ Donor,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      Alpha = "#4FA392",
      Beta = "#F4A261"
    )
  ) +
  scale_shape_manual(
    values = c(
      Untreated = 16,
      H2O2 = 17
    ),
    labels = c(
      Untreated = "Untreated",
      H2O2 = "H₂O₂"
    )
  ) +
  labs(
    title =
      "Pseudobulk multidimensional scaling",
    x =
      "MDS dimension 1",
    y =
      "MDS dimension 2",
    colour =
      "Cell type",
    shape =
      "Treatment"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    panel.grid.major = element_line(
      linewidth = 0.45,
      colour = "grey90"
    ),
    panel.grid.minor = element_line(
      linewidth = 0.25,
      colour = "grey95"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

ggsave(
  filename =
    "figures/Figure_S2_Pseudobulk_MDS.png",
  plot = pS2,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE S3
# PSEUDOBULK COUNT DEPTH
# ======================================================================

s3_data <- sample_info

s3_data$Condition <- factor(
  s3_data$Condition,
  levels = c(
    "Untreated",
    "H2O2"
  )
)

pS3 <- ggplot(
  s3_data,
  aes(
    x = Condition,
    y = library_millions,
    fill = Celltype
  )
) +
  geom_col(
    width = 0.62,
    colour = NA
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        library_millions
      )
    ),
    vjust = -0.45,
    size = 3.8
  ) +
  facet_grid(
    Donor ~ Celltype
  ) +
  scale_fill_manual(
    values = c(
      Alpha = "#4FA392",
      Beta = "#F4A261"
    )
  ) +
  scale_x_discrete(
    labels = c(
      Untreated = "Untreated",
      H2O2 = "H₂O₂"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.16
      )
    )
  ) +
  labs(
    title =
      "Pseudobulk count depth",
    x = NULL,
    y =
      "Pseudobulk library size (million counts)"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    strip.background = element_rect(
      fill = "white",
      colour = "grey75"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      linewidth = 0.45,
      colour = "grey90"
    ),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

ggsave(
  filename =
    "figures/Figure_S3_Pseudobulk_count_depth.png",
  plot = pS3,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE S4
# AVERAGE PSEUDOBULK COUNTS PER CELL
# ======================================================================

pS4 <- ggplot(
  sample_info,
  aes(
    x = Donor,
    y = counts_per_cell,
    colour = Celltype,
    shape = Condition
  )
) +
  geom_point(
    size = 5.2,
    position = position_dodge(
      width = 0.28
    )
  ) +
  scale_colour_manual(
    values = c(
      Alpha = "#4FA392",
      Beta = "#F4A261"
    )
  ) +
  scale_shape_manual(
    values = c(
      Untreated = 16,
      H2O2 = 17
    ),
    labels = c(
      Untreated = "Untreated",
      H2O2 = "H₂O₂"
    )
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  labs(
    title =
      "Average pseudobulk counts per cell",
    x =
      "Biological donor",
    y =
      "Pseudobulk counts per cell",
    colour =
      "Cell type",
    shape =
      "Treatment"
  ) +
  theme_bw(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      linewidth = 0.45,
      colour = "grey90"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

ggsave(
  filename =
    "figures/Figure_S4_Average_counts_per_cell.png",
  plot = pS4,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ======================================================================
# FIGURE S5
# GLOBAL PSEUDOBULK EXPRESSION STRUCTURE
#
# Top 1,000 most variable genes across all 12 pseudobulk samples.
# TMM-normalized log2 CPM.
# Row-wise z-score.
# Unsupervised hierarchical clustering.
# No DEG/FDR selection.
# ======================================================================

logCPM_all <- edgeR::cpm(
  y_all,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

gene_variance <- apply(
  logCPM_all,
  1,
  var
)

gene_variance <- gene_variance[
  is.finite(gene_variance) &
    gene_variance > 0
]

stopifnot(
  length(gene_variance) >= 1000
)

top1000_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[
  1:1000
]

heatmap_matrix_s5 <- logCPM_all[
  top1000_genes,
  ,
  drop = FALSE
]

heatmap_z_s5 <- t(
  scale(
    t(
      heatmap_matrix_s5
    )
  )
)

stopifnot(
  all(
    is.finite(
      heatmap_z_s5
    )
  )
)

sample_labels_s5 <- paste0(
  ifelse(
    sample_info$Celltype == "Alpha",
    "α",
    "β"
  ),
  "-",
  sample_info$Donor,
  "-",
  ifelse(
    sample_info$Condition == "H2O2",
    "H2O2",
    "Unt"
  )
)

colnames(heatmap_z_s5) <-
  sample_labels_s5

annotation_col_s5 <- data.frame(
  Cell_type = sample_info$Celltype,
  Treatment = sample_info$Condition,
  Donor = sample_info$Donor
)

rownames(annotation_col_s5) <-
  sample_labels_s5

annotation_colors_s5 <- list(
  Cell_type = c(
    Alpha = "#4FA392",
    Beta = "#F4A261"
  ),
  Treatment = c(
    Untreated = "#BDBDBD",
    H2O2 = "#7A6FAC"
  ),
  Donor = c(
    D1 = "#CC79A7",
    D2 = "#56B4E9",
    D3 = "#A6B83C"
  )
)

s5_breaks <- seq(
  -3,
  3,
  length.out = 102
)

s5_colors <- colorRampPalette(
  c(
    "#2C7FB8",
    "#F7F7F7",
    "#FDB863"
  )
)(101)

pheatmap::pheatmap(
  heatmap_z_s5,
  scale = "none",
  color = s5_colors,
  breaks = s5_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col_s5,
  annotation_colors = annotation_colors_s5,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 11,
  fontsize_col = 10,
  angle_col = 45,
  main =
    "Top 1,000 variable genes across pseudobulk samples",
  filename =
    "figures/Figure_S5_Top1000_variable_genes_heatmap.png",
  width = 9,
  height = 10
)


# ======================================================================
# FIGURE S6
# SINGLE-CELL UMAP OF THE EXACT 1,618-CELL MANUSCRIPT POPULATION
#
# Raw counts
# -> genes detected in >=10 cells
# -> library-size normalization to 10,000 counts/cell
# -> log1p
# -> top 2,000 variable genes
# -> gene-wise scaling
# -> PCA
# -> UMAP using first 30 PCs
#
# Exploratory visualization only.
# ======================================================================

counts_umap <- as(
  counts_ab,
  "dgCMatrix"
)

umap_meta <- meta

stopifnot(
  ncol(counts_umap) == 1618,
  nrow(umap_meta) == 1618,
  identical(
    colnames(counts_umap),
    as.character(
      umap_meta$cell_id
    )
  )
)

detected_cells <- Matrix::rowSums(
  counts_umap > 0
)

gene_keep_umap <- detected_cells >= 10

counts_filtered_umap <- counts_umap[
  gene_keep_umap,
  ,
  drop = FALSE
]

library_size_umap <- Matrix::colSums(
  counts_filtered_umap
)

stopifnot(
  all(
    library_size_umap > 0
  )
)

norm_counts_umap <- counts_filtered_umap %*%
  Matrix::Diagonal(
    x = 10000 /
      library_size_umap
  )

norm_counts_umap <- as(
  norm_counts_umap,
  "dgCMatrix"
)

norm_counts_umap@x <- log1p(
  norm_counts_umap@x
)

gene_mean_umap <- Matrix::rowMeans(
  norm_counts_umap
)

gene_mean_sq_umap <- Matrix::rowMeans(
  norm_counts_umap ^ 2
)

gene_variance_umap <-
  gene_mean_sq_umap -
  gene_mean_umap ^ 2

gene_variance_umap[
  gene_variance_umap < 0
] <- 0

n_variable_umap <- min(
  2000,
  length(
    gene_variance_umap
  )
)

top_variable_index_umap <- order(
  gene_variance_umap,
  decreasing = TRUE
)[
  seq_len(
    n_variable_umap
  )
]

variable_matrix_umap <- norm_counts_umap[
  top_variable_index_umap,
  ,
  drop = FALSE
]

variable_dense_umap <- as.matrix(
  variable_matrix_umap
)

variable_scaled_umap <- t(
  scale(
    t(
      variable_dense_umap
    )
  )
)

finite_gene_umap <- apply(
  variable_scaled_umap,
  1,
  function(z) {
    all(
      is.finite(z)
    )
  }
)

variable_scaled_umap <- variable_scaled_umap[
  finite_gene_umap,
  ,
  drop = FALSE
]

set.seed(
  218316
)

pca_result_umap <- prcomp(
  t(
    variable_scaled_umap
  ),
  center = FALSE,
  scale. = FALSE,
  rank. = 30
)

stopifnot(
  ncol(
    pca_result_umap$x
  ) >= 30
)

pca_coordinates_umap <- pca_result_umap$x[
  ,
  1:30,
  drop = FALSE
]

set.seed(
  218316
)

umap_coordinates <- uwot::umap(
  pca_coordinates_umap,
  n_neighbors = 30,
  min_dist = 0.30,
  metric = "euclidean",
  n_components = 2,
  init = "spectral",
  n_threads = 1,
  verbose = TRUE
)

umap_df <- data.frame(
  UMAP1 =
    umap_coordinates[
      ,
      1
    ],
  UMAP2 =
    umap_coordinates[
      ,
      2
    ],
  Celltype = factor(
    umap_meta$Celltype,
    levels = c(
      "Alpha",
      "Beta"
    )
  ),
  Condition = factor(
    umap_meta$Condition,
    levels = c(
      "Untreated",
      "H2O2"
    )
  ),
  Donor = factor(
    umap_meta$Donor,
    levels = c(
      "D1",
      "D2",
      "D3"
    )
  )
)

p_umap_celltype <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    colour = Celltype
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.75
  ) +
  scale_colour_manual(
    values = c(
      Alpha = "#4FA392",
      Beta = "#F4A261"
    )
  ) +
  labs(
    title =
      "Cell type",
    colour =
      "Cell type"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    legend.position = "bottom"
  )

p_umap_treatment <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    colour = Condition
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.75
  ) +
  scale_colour_manual(
    values = c(
      Untreated = "#BDBDBD",
      H2O2 = "#7A6FAC"
    ),
    labels = c(
      Untreated = "Untreated",
      H2O2 = "H₂O₂"
    )
  ) +
  labs(
    title =
      "Treatment",
    colour =
      "Treatment"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    legend.position = "bottom"
  )

p_umap_donor <- ggplot(
  umap_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    colour = Donor
  )
) +
  geom_point(
    size = 1.5,
    alpha = 0.75
  ) +
  scale_colour_manual(
    values = c(
      D1 = "#CC79A7",
      D2 = "#56B4E9",
      D3 = "#A6B83C"
    )
  ) +
  labs(
    title =
      "Donor",
    colour =
      "Donor"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    legend.position = "bottom"
  )

pS6 <- (
  p_umap_celltype |
    p_umap_treatment |
    p_umap_donor
) +
  patchwork::plot_annotation(
    title =
      "Single-cell transcriptional structure"
  ) &
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  filename =
    "figures/Figure_S6_Single_cell_UMAP.png",
  plot = pS6,
  width = 13,
  height = 5,
  dpi = 300,
  bg = "white"
)

write.csv(
  cbind(
    umap_meta,
    UMAP1 =
      umap_coordinates[
        ,
        1
      ],
    UMAP2 =
      umap_coordinates[
        ,
        2
      ]
  ),
  file =
    "figures/Figure_S6_UMAP_coordinates.csv",
  row.names = FALSE
)


# ======================================================================
# FINAL REPORT
# ======================================================================

cat(
  "\n========================================\n"
)

cat(
  "ALL FINAL FIGURES GENERATED\n"
)

cat(
  "========================================\n"
)

cat(
  "Figure 1 : figures/Figure_1_Shared_H2O2_gene_expression.png\n"
)

cat(
  "Figure 2 : figures/Figure_2_Differential_expression_volcano_plots.png\n"
)

cat(
  "Figure 3 : figures/Figure_3_Hallmark_transcriptional_signatures.png\n"
)

cat(
  "Figure S1: figures/Figure_S1_DEG_overlap.png\n"
)

cat(
  "Figure S2: figures/Figure_S2_Pseudobulk_MDS.png\n"
)

cat(
  "Figure S3: figures/Figure_S3_Pseudobulk_count_depth.png\n"
)

cat(
  "Figure S4: figures/Figure_S4_Average_counts_per_cell.png\n"
)

cat(
  "Figure S5: figures/Figure_S5_Top1000_variable_genes_heatmap.png\n"
)

cat(
  "Figure S6: figures/Figure_S6_Single_cell_UMAP.png\n"
)

cat(
  "\nUMAP settings:\n"
)

cat(
  "Cells: 1618\n"
)

cat(
  "Variable genes: 2000 maximum\n"
)

cat(
  "PCA dimensions: 30\n"
)

cat(
  "UMAP neighbors: 30\n"
)

cat(
  "UMAP min_dist: 0.30\n"
)

cat(
  "Random seed: 218316\n"
)

cat(
  "\nSession information:\n"
)

print(
  sessionInfo()
)
