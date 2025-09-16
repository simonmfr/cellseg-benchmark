# adapted from https://www.sc-best-practices.org/conditions/differential_gene_expression.html
library(dplyr)
library(data.table)
library(edgeR)
library(MAST)

edgeR_loop <- function(adata, group_i, edger_methods, test_groups, ref_group="WT",
                       condition_col="genotype", batch_col=NULL) {

  condition <- droplevels(factor(colData(adata)[[condition_col]]))
  if (is.null(condition) || nlevels(condition) < 2) { message(group_i, ": skip (condition invalid or <2 levels)"); return(NULL) }

  present <- intersect(test_groups, levels(condition))
  if (!(ref_group %in% levels(condition)) || !length(present)) { message(group_i, ": skip (ref_group/test_groups missing)"); return(NULL) }

  missing <- setdiff(test_groups, present)
  message(group_i, ": ", paste(present, collapse="/"), " vs ", ref_group,
          if (length(missing)) paste0(" (missing:", paste(missing, collapse=","), ")") else "")

  res <- list()
  for (m in edger_methods) {
    o <- edgeR_fit_model(adata, edger_method=m, condition_col=condition_col, ref_group=ref_group, batch_col=batch_col)
    if (is.null(o)) next
    for (tg in present) {
      tmp <- edgeR_run_test(o$fit, o$design, m, tg, ref_group)
      if (is.null(tmp)) next
      tmp <- dplyr::mutate(as.data.frame(tmp), FC = 2^logFC, edgeR_method = m, test_group = tg, ref = ref_group)
      res[[paste(tg, m, sep="_")]] <- tmp
    }
  }
  if (!length(res)) return(NULL)
  out <- dplyr::bind_rows(res, .id="result_id")
  out$gene <- sub("[.][.][.].*", "", rownames(out)); rownames(out) <- NULL
  dplyr::select(out, gene, FC, dplyr::everything())
  #out$subset <- as.character(group_i)
  #out <- dplyr::relocate(out, subset, .before = 1)
}

edgeR_fit_model <- function(
    adata, 
    edger_method="LRT",
    condition_col="genotype",
    ref_group=NULL,
    batch_col=NULL,
    min_count=2
) {
  condition <- droplevels(factor(colData(adata)[[condition_col]]))
  condition <- stats::relevel(condition, ref = ref_group)
  
  # us batch covariate if exists and has >1 level
  use_batch <- !is.null(batch_col) && batch_col %in% colnames(colData(adata)) &&
               nlevels(droplevels(factor(colData(adata)[[batch_col]]))) > 1
  
  design <- if (use_batch) {
    batch <- droplevels(factor(colData(adata)[[batch_col]]))
    model.matrix(~ condition + batch)
  } else {
    model.matrix(~ condition)
  }
  
  message("  Design: ~ ", condition_col, if (use_batch) paste0(" + ", batch_col) else "")
  
  rdof <- nrow(design) - qr(design)$rank
  if (rdof <= 0) {
    message("  Skip: no residual df (likely too few samples or condition confounded with batch)")
    return(NULL)
  }
    
  y <- edgeR::DGEList(assay(adata, "X"), group = condition)
  keep <- edgeR::filterByExpr(y, min.count = min_count, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # replace library sizes with volume
  vol <- colData(adata)[["volume_sum"]]
  if (is.null(vol) || any(vol <= 0)) stop("Invalid volume in colData(sce), expects 'volume_sum'.")
  y$samples$lib.size <- vol
  y <- edgeR::calcNormFactors(y)
    
  fit <- if (edger_method == "QL") {
    edgeR::glmQLFit(edgeR::estimateDisp(y, design), design)
  } else if (edger_method == "LRT") {
    edgeR::glmFit(edgeR::estimateGLMRobustDisp(y, design), design)
  }
  list(fit=fit, design=design)
}

edgeR_run_test <- function(fit, design, edger_method, test_group, ref_group){
  
  c_test <- paste0("condition", make.names(test_group))
  c_ref  <- paste0("condition", make.names(ref_group))

  if (c_test %in% colnames(design) && !(c_ref %in% colnames(design))) {
    # intercept design
    tt <- if (edger_method == "QL") edgeR::glmQLFTest(fit, coef = c_test)
          else                      edgeR::glmLRT(fit,  coef = c_test)

  } else if (all(c(c_test, c_ref) %in% colnames(design))) {
    # no-intercept design
    con <- limma::makeContrasts(contrasts = sprintf("%s-%s", c_test, c_ref), levels = design)
    tt  <- if (edger_method == "QL") edgeR::glmQLFTest(fit, contrast = con)
           else                      edgeR::glmLRT(fit,  contrast = con)

  } else {
    stop("Cannot form contrast. Needed ", c_test,
         if (c_ref %in% colnames(design)) "" else paste0(" and ", c_ref),
         ". Available columns: ", paste(colnames(design), collapse=", "))
  }

  de <- edgeR::topTags(tt, n = Inf)$table
  de$test <- paste0(test_group, "vs", ref_group)
  de
}

mastre_run <- function(adata,
                        group_i,
                        test_groups,
                        ref_group,
                        condition_col,
                        sample_col,
                        batch_col) {
  
  sca <- MAST::SceToSingleCellAssay(adata)

  colData(sca)$Condition <- factor(colData(sca)[[condition_col]])
  colData(sca)$Condition <- stats::relevel(colData(sca)$Condition, ref = ref_group)
  colData(sca)$Sample    <- factor(colData(sca)[[sample_col]])
  colData(sca)$Batch     <- factor(colData(sca)[[batch_col]])

  message("Fitting hurdle model: ~ ", condition_col, " + ", batch_col, " + (1 | ", sample_col, ")...")
  fit <- MAST::zlm(~ Condition + Batch + (1 | Sample),
                   sca,
                   method = "glmer",
                   ebayes = FALSE,
                   strictConvergence = FALSE)

  message("Collecting contrasts...")
  res_list <- lapply(test_groups, function(g) {
    contrast_i <- paste0("Condition", g)
    summaryCond <- MAST::summary(fit, doLRT = contrast_i)
    dt <- summaryCond$datatable

    p_main <- dt[component == "H" & contrast == contrast_i,
                 .(gene = primerid, PValue = `Pr(>Chisq)`)]
    logfc  <- dt[component == "logFC" & contrast == contrast_i,
                 .(gene = primerid, log2FC = coef, wald = z)]

    de <- merge(p_main, logfc, by = "gene", all.x = TRUE)
    de[, FDR := p.adjust(PValue, method = "BH")]
    setorder(de, PValue)

    de[, `:=`(
      subset_group = group_i,
      method       = "MAST-RE",
      test         = paste0(g, "vs", ref_group),
      test_group   = g,
      ref          = ref_group,
      FC           = 2^log2FC
    )]
    de
  })

  if (length(res_list) == 0) return(NULL)
  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  setDF(res)
  return(res)
}
