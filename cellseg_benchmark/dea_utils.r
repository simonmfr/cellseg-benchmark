suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(edgeR)   
  library(variancePartition)
  library(limma)
  library(BiocParallel)
  library(MAST)
})

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

dream_fit_model <- function(
  adata,
  condition_col = "genotype",
  ref_group     = "WT",
  batch_col     = NULL,
  threads       = 4,
  min_count     = 2
){
  df <- as.data.frame(colData(adata))

  condition <- droplevels(factor(df[[condition_col]]))
  if (is.null(condition) || nlevels(condition) < 2) stop("condition invalid or <2 levels")
  if (!(ref_group %in% levels(condition)))         stop("ref_group not in condition levels")
  condition <- stats::relevel(condition, ref = ref_group)
  df$condition <- condition

  use_batch <- !is.null(batch_col) &&
               batch_col %in% colnames(df) &&
               nlevels(droplevels(factor(df[[batch_col]]))) > 1
  if (use_batch) {
    df$batch <- droplevels(factor(df[[batch_col]]))
    form <- ~ condition + batch
  } else {
    form <- ~ condition
  }

  design <- model.matrix(form, df)
  rk   <- qr(design)$rank
  rdof <- nrow(design) - rk

  summarize_counts <- function(df){
    out <- c()
    if ("condition" %in% colnames(df)) {
      tc <- table(df$condition)
      out <- c(out, paste0("conditions [", length(tc), "]: ",
                           paste(sprintf("%s (n=%d)", names(tc), tc), collapse=", ")))
    }
    if ("batch" %in% colnames(df)) {
      tb <- table(df$batch)
      out <- c(out, paste0("batches [", length(tb), "]: ",
                           paste(sprintf("batch %s (n=%d)", names(tb), tb), collapse=", ")))
    }
    paste(out, collapse=" | ")
  }

  if (rk < ncol(design)) {
    ne <- limma::nonEstimable(design)
    msg_counts <- summarize_counts(df)
    message("  Skip: design not full rank. Non-estimable coef: ",
            paste(ne, collapse=", "), " | ", msg_counts)
    return(NULL)
  }

  if (rdof <= 0) {
    msg_counts <- summarize_counts(df)
    message("  Skip: no residual df. | ", msg_counts)
    return(NULL)
  }
      
  y <- edgeR::DGEList(assay(adata, "X"), group = df$condition)
  keep <- edgeR::filterByExpr(y, min.count = min_count, design = model.matrix(form, df))
  y <- y[keep, , keep.lib.sizes = FALSE]

  vol <- colData(adata)[["volume_sum"]]
  if (is.null(vol) || any(vol <= 0)) stop("Invalid 'volume_sum' in colData(adata)")
  y$samples$lib.size <- vol
  y <- edgeR::calcNormFactors(y)

  bp   <- BiocParallel::SnowParam(as.numeric(threads), "SOCK", progressbar = TRUE)
  vobj <- variancePartition::voomWithDreamWeights(y, form, df, BPPARAM = bp)
  fit  <- variancePartition::dream(vobj, form, df)
  fit  <- limma::eBayes(fit)

  list(fit = fit, ref_group = ref_group)
}

dream_extract <- function(fit_obj, test_groups){
  fit       <- fit_obj$fit
  ref_group <- fit_obj$ref_group

  res_list <- lapply(test_groups, function(tg){
    coef_name <- paste0("condition", make.names(tg))
    if (!(coef_name %in% colnames(coef(fit)))) return(NULL)

    tt <- variancePartition::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
    if (!nrow(tt)) return(NULL)

    data.table::setnames(tt, old = c("P.Value", "adj.P.Val"),
                             new = c("PValue",  "FDR"),
                             skip_absent = TRUE)
      
    tt <- as.data.table(tt, keep.rownames = "gene") |>
      dplyr::mutate(
        FC         = 2^logFC,
        method     = "DREAM",
        test_group = tg,
        ref        = ref_group,
        test       = paste0(tg, "vs", ref_group)
      )

    data.table::setcolorder(tt, c("gene","FC", setdiff(names(tt), c("gene","FC"))))
    tt
  })

  data.table::rbindlist(Filter(Negate(is.null), res_list))
}

dream_loop <- function(
  adata,
  group_i,
  test_groups,
  ref_group     = "WT",
  condition_col = "genotype",
  batch_col     = NULL,
  threads       = 4,
  min_count     = 2
){
  condition <- droplevels(factor(colData(adata)[[condition_col]]))
  if (is.null(condition) || nlevels(condition) < 2) {
    message(group_i, ": skip (condition invalid or <2 levels)"); return(NULL)
  }
  present <- intersect(test_groups, levels(condition))
  if (!(ref_group %in% levels(condition)) || !length(present)) {
    message(group_i, ": skip (ref_group/test_groups missing)"); return(NULL)
  }
  missing <- setdiff(test_groups, present)
  message(group_i, ": ", paste(present, collapse="/"), " vs ", ref_group,
          if (length(missing)) paste0(" (missing:", paste(missing, collapse=","), ")") else "")

  o <- tryCatch(
    dream_fit_model(adata, condition_col, ref_group, batch_col,
                    threads = threads, min_count = min_count),
    error = function(e){ message("  Skip: ", e$message); return(NULL) }
  )
  if (is.null(o)) return(NULL)

  dream_extract(o, present)
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
