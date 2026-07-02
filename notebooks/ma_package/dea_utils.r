suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(variancePartition)
  library(limma)
  library(BiocParallel)
  library(MAST)
  library(SummarizedExperiment)
  library(Matrix)
})

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

  ## concise tag (use user input, not expanded list)
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

  ## factor / droplevels hygiene
  cd[[condition_col]] <- droplevels(factor(cd[[condition_col]]))
  cd[[sample_col]]    <- droplevels(factor(cd[[sample_col]]))
  if (!is.null(batch_col) && batch_col %in% colnames(cd)) {
    cd[[batch_col]] <- droplevels(factor(cd[[batch_col]]))
  }

  ## relevel condition and write back
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

  cd  <- SummarizedExperiment::colData(sca)
  grp <- droplevels(factor(cd[[condition_col]]))

  n_cells <- ncol(sca)
  n_sample_levels <- nlevels(cd[[sample_col]])

  ## 2a) global size check
  min_cells = 250L # required for GLMM stability
  if (n_cells < min_cells) {
    message("   Skipping group '", group_i,
            "': too few cells for a mixed-effects model (",
            n_cells, " cells across conditions).")
    return(NULL)
  }

  ## 2b) drop tiny samples, then re-check
  tab_sample <- table(cd[[sample_col]])
  small_samples <- names(tab_sample[tab_sample < 10L])

  if (length(small_samples) > 0L) {
    message("   Removing ", length(small_samples),
            " samples with <10 cells: ",
            paste(small_samples, collapse = ", "))

    keep_idx <- !(cd[[sample_col]] %in% small_samples)
    sca <- sca[, keep_idx, drop = FALSE]

    cd  <- SummarizedExperiment::colData(sca)
    grp <- droplevels(factor(cd[[condition_col]]))
    n_cells <- ncol(sca)

    if (n_cells < min_cells) {
      message("   Skipping group '", group_i,
              "': too few cells after removing samples with few cells (",
              n_cells, " cells).")
      return(NULL)
    }
  }

  if (n_sample_levels <= 1L || n_sample_levels >= n_cells) {
    message("   Skipping group '", group_i,
            "': random effect '", sample_col,
            "' not estimable (", n_sample_levels,
            " levels, ", n_cells, " cells).")
    return(NULL)
  }

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

  fit <- tryCatch(
    {
      MAST::zlm(form, sca, method = zlm_method,
                ebayes = FALSE, strictConvergence = FALSE, parallel = TRUE)
    },
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("Downdated VtV is not positive definite", msg, fixed = TRUE)) {
        message("   Skipping group '", group_i,
                "': mixed-effects model did not converge (non–positive definite random-effects covariance).")
        return(NULL)
      }
      stop(e)
    }
  )

  if (is.null(fit)) {
    return(NULL)
  }

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