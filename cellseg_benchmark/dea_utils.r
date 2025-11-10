suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(edgeR)   
  library(variancePartition)
  library(limma)
  library(BiocParallel)
  library(MAST)
  library(SummarizedExperiment)
  library(Matrix)
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
      met <- paste0("edgeR_", m)
      tmp <- dplyr::mutate(as.data.frame(tmp), FC = 2^logFC, method = met, test_group = tg, ref = ref_group)
      res[[paste(tg, met, sep="_")]] <- tmp
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

mast_run_cached <- function(adata,
                            group_i,
                            test_groups,
                            ref_group,
                            condition_col,
                            sample_col,
                            batch_col = NULL,
                            cache_dir = ".",
                            overwrite = FALSE,
                            clear_cache = TRUE,
                            sample_random_effects = TRUE,
                            max_cells_per_condition = 10000,
                            brain_region_subset = NULL) {

  mode_tag <- if (isTRUE(sample_random_effects)) "RE" else "FE"

  ## interpret cap
  unlimited_cap <- is.null(max_cells_per_condition) || !is.finite(max_cells_per_condition)
  cap_val <- if (unlimited_cap) "ALL" else as.integer(max_cells_per_condition)

  ## --- concise tag (use user input, not expanded list) ---
  br_tag <- if (!is.null(brain_region_subset) && nzchar(brain_region_subset))
    paste0("_Subset-", gsub("\\s+", "", brain_region_subset)) else ""

  ## paths (cap- & subset-aware)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(cache_dir, paste0("DE_", group_i, "_", mode_tag, "_cap", cap_val, br_tag))
  res_path <- paste0(prefix, "_res.rds")
  contrast_path <- function(g) paste0(prefix, "_", g, ".rds")
  log_file <- paste0(prefix, "_R.log")

  ## fast return if cached
  if (file.exists(res_path) && !overwrite) {
    message("-> Using cached final result: ", res_path)
    out <- readRDS(res_path)
    return(out)
  }

  ## logging (only for fresh run)
  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output"); sink(log_con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

  t_all_start <- Sys.time()

  ## 1) Prep + filter
  t1_start <- Sys.time()
  message("-> Step 1: preparing SCA")
  sca <- MAST::SceToSingleCellAssay(adata)
  cd <- SummarizedExperiment::colData(sca)
  cd[[condition_col]] <- factor(cd[[condition_col]])
  cd[[sample_col]]    <- factor(cd[[sample_col]])
  if (!ref_group %in% levels(cd[[condition_col]]))
    stop("ref_group not found in condition_col")
  cd[[condition_col]] <- stats::relevel(cd[[condition_col]], ref = ref_group)
  SummarizedExperiment::colData(sca) <- cd

  message("   total cells: ", ncol(sca),
          " | conditions: ", length(unique(cd[[condition_col]])))
  message("   max_cells_per_condition: ",
          if (unlimited_cap) "unlimited" else cap_val)

  ## subsample per condition if needed (skip if unlimited)
  if (!unlimited_cap && ncol(sca) > cap_val) {
    set.seed(1L)
    cond_now <- SummarizedExperiment::colData(sca)[[condition_col]]
    idx_keep <- unlist(lapply(
      split(seq_len(ncol(sca)), cond_now),
      function(ix)
        if (length(ix) > cap_val)
          sample(ix, cap_val)
        else ix
    ), use.names = FALSE)
    new_n <- length(unique(idx_keep))
    if (new_n < ncol(sca)) {
      sca <- sca[, sort(unique(idx_keep))]
      message("   subsampled to ", new_n, " cells total")
    }
  }

  ## recompute after possible subsampling
  cd  <- SummarizedExperiment::colData(sca)
  grp <- droplevels(factor(cd[[condition_col]]))

  ## gene filtering
  Y  <- SummarizedExperiment::assay(sca)
  nz <- as(Y != 0, "dMatrix")
  det_overall <- Matrix::rowMeans(nz)
  det_by_grp <- lapply(levels(grp), function(lvl) {
    idx <- which(!is.na(grp) & grp == lvl)
    Matrix::rowMeans(nz[, idx, drop = FALSE])
  })
  det_any <- do.call(pmax, det_by_grp)

  message("   filtering genes: keep if detected in ≥10% of cells in any condition OR ≥20% overall")
  keep <- (det_any >= 0.10) | (det_overall >= 0.20)
  sca  <- sca[keep, ]
  message("   kept: ", sum(keep), " / ", nrow(Y), " genes")

  t1_end <- Sys.time()
  message(sprintf("   [Step 1 done in %.2f min]", as.numeric(difftime(t1_end, t1_start, units = "mins"))))

  ## 2) Model
  t2_start <- Sys.time()
  cd <- SummarizedExperiment::colData(sca)
  use_batch <- !is.null(batch_col) &&
    batch_col %in% colnames(cd) &&
    nlevels(droplevels(factor(cd[[batch_col]]))) > 1

  if (isTRUE(sample_random_effects)) {
    model_str <- if (use_batch) {
      paste0("~ ", condition_col, " + ", batch_col, " + (1 | ", sample_col, ")")
    } else {
      paste0("~ ", condition_col, " + (1 | ", sample_col, ")")
    }
    zlm_method <- "glmer"
  } else {
    model_str <- if (use_batch) {
      paste0("~ ", condition_col, " + ", batch_col)
    } else {
      paste0("~ ", condition_col)
    }
    zlm_method <- "bayesglm"
  }

  message("-> Step 2: fitting model (", model_str, "; method = ", zlm_method, ")")
  form <- as.formula(model_str)
  fit <- MAST::zlm(form, sca, method = zlm_method,
                   ebayes = FALSE, strictConvergence = FALSE, parallel = TRUE)
  t2_end <- Sys.time()
  message(sprintf("   [Step 2 done in %.2f min]", as.numeric(difftime(t2_end, t2_start, units = "mins"))))

  ## 3) Contrasts
  t3_start <- Sys.time()
  message("-> Step 3: contrasts")
  groups_to_test <- test_groups[vapply(test_groups, \(g)
    !file.exists(contrast_path(g)) || overwrite, logical(1))]

  res_list <- lapply(setNames(test_groups, test_groups), function(g) {
    if (!g %in% groups_to_test) {
      message("   using cached: ", g)
      return(readRDS(contrast_path(g)))
    }
    NULL
  })

  if (length(groups_to_test) > 0) {
    message("   computing: ", paste(groups_to_test, collapse = ", "))
    avail_coefs <- if ("coefD" %in% slotNames(fit)) colnames(fit@coefD) else character()
    cond_coefs  <- grep(paste0("^", condition_col), avail_coefs, value = TRUE)
    cond_groups <- sub(paste0("^", condition_col), "", cond_coefs)
    coefs_to_test <- cond_coefs[cond_groups %in% groups_to_test]

    if (length(coefs_to_test) > 0) {
      dt_all <- MAST::summary(fit, doLRT = coefs_to_test, parallel = FALSE)$datatable
      p_all   <- dt_all[component == "H",     .(gene = primerid, contrast, PValue = `Pr(>Chisq)`)]
      lfc_all <- dt_all[component == "logFC", .(gene = primerid, contrast, log2FC = coef, wald = z)]
      merged  <- merge(p_all, lfc_all, by = c("gene", "contrast"), all.x = TRUE)

      for (g in groups_to_test) {
        ci <- paste0(condition_col, g)
        de <- merged[contrast == ci]
        if (!nrow(de)) next
        de[, FDR := p.adjust(PValue, "BH")]
        data.table::setorder(de, PValue)
        de[, `:=`(
          subset_group = group_i,
          method = paste0("MAST-", mode_tag),
          test = paste0(g, "vs", ref_group),
          test_group = g,
          ref = ref_group,
          FC = 2^log2FC,
          model = model_str
        )]
        saveRDS(as.data.frame(de), contrast_path(g))
        res_list[[g]] <- de
      }
    } else {
      message("   no valid contrasts for ", condition_col)
    }
  }

  ## 4) Combine + save
  res <- data.table::rbindlist(res_list, use.names = TRUE)
  if ("contrast" %in% names(res)) res[, contrast := NULL]
  if (!"model" %in% names(res)) res[, model := model_str]
  saveRDS(as.data.frame(res), res_path)
  message("-> Done. Saved: ", res_path)

  ## 5) Cleanup
  unlink(vapply(test_groups, contrast_path, character(1)), force = TRUE)
  if (clear_cache) {
    unlink(res_path, force = TRUE)
    message("-> Final res cache cleared (clear_cache = TRUE)")
  } else {
    message("-> Final res cache kept (clear_cache = FALSE)")
  }

  t3_end <- Sys.time()
  message(sprintf("   [Step 3 done in %.2f min]", as.numeric(difftime(t3_end, t3_start, units = "mins"))))
  t_all_end <- Sys.time()
  message(sprintf("[Total runtime %.2f min]", as.numeric(difftime(t_all_end, t_all_start, units = "mins"))))

  return(res)
}