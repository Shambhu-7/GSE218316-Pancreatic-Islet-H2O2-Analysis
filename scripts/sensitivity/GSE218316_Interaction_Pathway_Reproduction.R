#!/usr/bin/env Rscript
# Interaction-level Hallmark fgsea + CAMERA for GSE218316.
# Run from the root of GSE218316-Pancreatic-Islet-H2O2-Analysis.
# Does not alter existing manuscript/reference outputs.
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(fgsea); library(msigdbr)
})
stopifnot(file.exists('GSE218316_Reproducible_Analysis.R'))
# The verified main script reconstructs pseudobulk counts and sample annotation.
source('GSE218316_Reproducible_Analysis.R', local=FALSE)
stopifnot(exists('pb'), exists('joint_info'), exists('interaction_results'))
stopifnot(ncol(pb)==12L)

# Retrieve human Hallmark sets and record the exact release from the installed package.
ms <- msigdbr::msigdbr(species='Homo sapiens', collection='H')
stopifnot(all(c('gs_name','gene_symbol') %in% names(ms)))
sets <- split(ms$gene_symbol, ms$gs_name)
sets <- lapply(sets, unique)

# Resolve feature ID -> gene symbol from the original annotation, preserving ID keys.
ann <- readRDS('data/GSE218316/gene_annotation.rds')
if (is.data.frame(ann)) {
  cn <- names(ann)
  id_col <- intersect(c('feature_id','gene_id','ensembl_id','GeneID','id'), cn)
  sym_col <- intersect(c('gene_symbol','symbol','gene_name','GeneSymbol','Symbol'), cn)
  if (!length(sym_col)) stop('Cannot identify gene-symbol column in gene_annotation.rds: ',paste(cn,collapse=', '))
  if (length(id_col)) {
    symmap <- setNames(as.character(ann[[sym_col[1]]]), as.character(ann[[id_col[1]]]))
  } else if (!is.null(rownames(ann))) {
    symmap <- setNames(as.character(ann[[sym_col[1]]]), rownames(ann))
  } else stop('Annotation has no identifiable feature IDs')
} else if (is.atomic(ann) && !is.null(names(ann))) {
  symmap <- as.character(ann)
  names(symmap) <- names(ann)
} else stop('Unsupported annotation structure: ',paste(class(ann),collapse=', '))
if (sum(rownames(pb) %in% names(symmap)) < 5000L)
  stop('Annotation IDs do not match count row IDs; inspect annotation structure before proceeding')

joint_info$Donor <- factor(joint_info$Donor)
joint_info$Celltype <- factor(joint_info$Celltype, levels=c('Alpha','Beta'))
joint_info$Condition <- factor(joint_info$Condition, levels=c('Untreated','H2O2'))
if (anyNA(joint_info[,c('Donor','Celltype','Condition')])) stop('Missing joint-model factor levels')

outdir <- 'results/interaction_pathway_reproduction'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

run_model <- function(name, formula) {
  cat('\n---',name,'---\n')
  design <- model.matrix(formula, data=joint_info)
  stopifnot(qr(design)$rank==ncol(design), nrow(design)==ncol(pb))
  contrast_col <- 'CelltypeBeta:ConditionH2O2'
  if (!contrast_col %in% colnames(design)) stop('Missing interaction coefficient: ',paste(colnames(design),collapse=', '))
  y <- edgeR::DGEList(counts=pb)
  keep <- edgeR::filterByExpr(y,design=design)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y,method='TMM')
  y <- edgeR::estimateDisp(y,design=design)
  fit <- edgeR::glmQLFit(y,design=design,robust=TRUE)
  de <- edgeR::topTags(edgeR::glmQLFTest(fit,coef=contrast_col),n=Inf)$table
  de$feature_id <- rownames(de)
  de$gene_symbol <- unname(symmap[de$feature_id])
  cat('Tested:',nrow(de),' significant:',sum(de$FDR<0.05,na.rm=TRUE),' min FDR:',min(de$FDR,na.rm=TRUE),'\n')
  write.csv(de,file.path(outdir,paste0(name,'_interaction_DE.csv')),row.names=FALSE)

  # Mirror the verified script's signed-ranking and duplicate resolution convention:
  # decreasing signed rank, first feature per gene symbol (not maximum absolute rank).
  rankdf <- de[!is.na(de$gene_symbol) & nzchar(de$gene_symbol) & is.finite(de$F),]
  rankdf$rank_score <- sign(rankdf$logFC)*sqrt(pmax(rankdf$F,0))
  rankdf <- rankdf[order(-rankdf$rank_score,rankdf$feature_id),]
  rankdf <- rankdf[!duplicated(rankdf$gene_symbol),]
  ranks <- setNames(rankdf$rank_score,rankdf$gene_symbol)
  ranks <- sort(ranks,decreasing=TRUE)
  set.seed(218316)
  fg <- fgsea::fgseaMultilevel(pathways=sets,stats=ranks,minSize=15,maxSize=500)
  fg <- as.data.frame(fg)
  if ('leadingEdge' %in% names(fg)) fg$leadingEdge <- vapply(fg$leadingEdge,paste,collapse=';',FUN.VALUE=character(1))
  fg <- fg[order(fg$padj,fg$pval),]
  write.csv(fg,file.path(outdir,paste0(name,'_interaction_fgsea_ALL.csv')),row.names=FALSE)

  # CAMERA competitive gene-set test on voom-transformed joint pseudobulk data.
  # One feature per gene symbol, chosen by average expression independently of effect.
  sy <- unname(symmap[rownames(y)])
  good <- !is.na(sy) & nzchar(sy)
  ycam <- y[good,,keep.lib.sizes=FALSE]
  sy <- sy[good]
  av <- rowMeans(edgeR::cpm(ycam,log=TRUE,prior.count=2))
  ord <- order(sy,-av)
  choose <- ord[!duplicated(sy[ord])]
  ycam <- ycam[choose,,keep.lib.sizes=FALSE]
  sy <- sy[choose]
  rownames(ycam$counts) <- sy
  ycam <- edgeR::calcNormFactors(ycam)
  vv <- limma::voom(ycam,design=design,plot=FALSE)
  ix <- lapply(sets,function(g) which(rownames(vv$E) %in% g))
  ix <- ix[vapply(ix,length,integer(1))>=15L]
  coef_idx <- match(contrast_col,colnames(design))
  cam_est <- limma::camera(vv,index=ix,design=design,contrast=coef_idx,inter.gene.cor=NA,sort=FALSE)
  cam_est$Pathway <- rownames(cam_est)
  cam_est$FDR <- p.adjust(cam_est$PValue,method='BH')
  write.csv(cam_est,file.path(outdir,paste0(name,'_interaction_CAMERA_estimated_ALL.csv')),row.names=FALSE)
  cam_fix <- limma::camera(vv,index=ix,design=design,contrast=coef_idx,inter.gene.cor=0.01,sort=FALSE)
  cam_fix$Pathway <- rownames(cam_fix)
  cam_fix$FDR <- p.adjust(cam_fix$PValue,method='BH')
  write.csv(cam_fix,file.path(outdir,paste0(name,'_interaction_CAMERA_fixed001_ALL.csv')),row.names=FALSE)

  for (method in c('fgsea','CAMERA estimated','CAMERA fixed 0.01')) {
    z <- switch(method,fgsea=fg,'CAMERA estimated'=cam_est,'CAMERA fixed 0.01'=cam_fix)
    namecol <- if (method=='fgsea') 'pathway' else 'Pathway'
    ox <- z[grepl('OXIDATIVE_PHOSPHORYLATION',z[[namecol]],fixed=TRUE),,drop=FALSE]
    cat('\n',method,' OXPHOS:\n',sep='')
    print(ox[,intersect(c(namecol,'NES','pval','padj','Direction','PValue','FDR','Correlation'),names(ox)),drop=FALSE],row.names=FALSE)
  }
  invisible(list(de=de,fgsea=fg,camera_est=cam_est,camera_fixed=cam_fix))
}

primary <- run_model('primary_shared_donor',~ Donor + Celltype * Condition)
sensitivity <- run_model('sensitivity_donor_by_celltype',~ Donor * Celltype + Condition * Celltype)
writeLines(capture.output(sessionInfo()),file.path(outdir,'interaction_pathway_sessionInfo.txt'))
writeLines(c(
  'Models: primary ~ Donor + Celltype * Condition; sensitivity ~ Donor * Celltype + Condition * Celltype',
  paste('msigdbr:',as.character(packageVersion('msigdbr'))),
  paste('fgsea:',as.character(packageVersion('fgsea'))),
  paste('limma:',as.character(packageVersion('limma'))),
  paste('edgeR:',as.character(packageVersion('edgeR'))),
  'fgsea: minSize=15, maxSize=500, default eps; ranks sign(logFC)*sqrt(F), max signed rank per symbol.',
  'CAMERA: voom, one feature per symbol by highest mean logCPM, estimated or fixed correlation 0.01.',
  'IMPORTANT: inspect outputs and compare against manuscript; do not assume previously reported OXPHOS values match.'
),file.path(outdir,'README_methods.txt'))
cat('\nCOMPLETE. Results written to:',normalizePath(outdir),'\n')
