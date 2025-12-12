# Formula:
#   score ~ surgery + tissue_source + (1 | patient)
# 
#   - score is the GBMDeonvoluteR score for a given cell type
#   - surgery tests the score change (primary → recurrent)
#   - tissue_source adjusts for known confounding factors in Nomura et al
#   - (1 | patient) handles paired samples by including patient as a random effect

fit_celltype_model <- function(long_df,
                                   cell_type_value,
                                   cell_type_col = "cell_type",
                                   model_formula = score ~ surgery + tissue_source + (1 | patient),
                                   coef_name = "surgeryRecurrent") {
 
    dat <- long_df[long_df[[cell_type_col]] == cell_type_value, , drop = FALSE]
    n <- nrow(dat)
    
    # Base result row with NAs filled in; we always return one row per cell type
    result_row <- data.frame(
        cell_type = cell_type_value,
        term = coef_name,
        estimate = NA_real_,
        se = NA_real_,
        df = NA_real_,
        t = NA_real_,
        p = NA_real_,
        n = n,
        model_status = NA_character_,
        stringsAsFactors = FALSE
    )
    
    # Not enough rows to fit a model
    if (n < 3) {
        result_row$model_status = "insufficient_rows"
        return(list(result = result_row))
    }
    
    # Try fitting model
    fit <- try(lmerTest::lmer(model_formula, data = dat), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
        result_row$model_status = "fit_error"
        return(list(result = result_row, error = as.character(fit)))
    }
    
    coef_tab <- summary(fit)$coef
    
    if (!coef_name %in% rownames(coef_tab)) {
        result_row$model_status = "coef_not_found"
        return(list(result = result_row, fit = fit))
    }
    
    this <- coef_tab[coef_name, , drop = FALSE]
    
    result_row$estimate = this[,"Estimate"]
    result_row$se       = this[,"Std. Error"]
    result_row$df       = this[,"df"]
    result_row$t        = this[,"t value"]
    result_row$p        = this[,"Pr(>|t|)"]
    result_row$model_status = "ok"
    
    return(list(result = result_row, fit = fit))
}


run_celltype_mixed_models = function(long_df,
                                     cell_type_col = "cell_type",
                                     celltypes = NULL,
                                     model_formula = score ~ surgery + tissue_source + (1 | patient),
                                     coef_name = "surgeryRecurrent",
                                     p_adjust_method = "BH",
                                     p_threshold = 0.05
                                     ) {
    
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' is required. Please install it with install.packages('dplyr').")
    }
    
    if (!cell_type_col %in% colnames(long_df)) {
        stop("Column '", cell_type_col, "' not found in long_df.")
    }
    
    # Do NOT silently coerce to factor
    if (!is.factor(long_df[[cell_type_col]])) {
        stop("Column '", cell_type_col, "' must be a factor. ",
             "Please convert it (e.g. long_df$", cell_type_col, " = factor(long_df$", cell_type_col, ")) before calling.")
    }
    
    data_celltypes = levels(long_df[[cell_type_col]])
    if (is.null(celltypes)) {
        celltypes = data_celltypes
    } else {
        # Check that all unique values present in the data are covered by supplied vector
        if (!all(celltypes %in% data_celltypes)) {
            missing_vals = setdiff(data_celltypes, celltypes)
            stop(
                "The supplied 'celltypes' vector does not include all levels of '", cell_type_col, "'. Missing: ",
                paste(missing_vals, collapse = ", ")
            )
        }
    }
    
    fits_list <-  lapply(celltypes, function(ct) {
        fit_celltype_model(
            long_df = long_df,
            cell_type_value = ct,
            cell_type_col = cell_type_col,
            model_formula = model_formula,
            coef_name = coef_name
        )
    })
    
    # Combine result rows (one per cell type) with bind_rows
    result_rows <- dplyr::bind_rows(lapply(fits_list, `[[`, "result"))
    
    if (nrow(result_rows) == 0) {
        warning("No results were produced; check your inputs.")
        return(list(
            results = result_rows,
            by_cell_type = fits_list
        ))
    }
    
    # Compute p-adjust, log2FC etc. only for successfully fitted models
    ok_idx <- which(result_rows$model_status == "ok")
    
    if (length(ok_idx) > 0) {
        result_rows$padj <- NA_real_
        result_rows$log2FC <- NA_real_
        result_rows$significant <- NA
        result_rows$neglog10_padj <- NA_real_
        
        result_rows$padj[ok_idx] <- 
            p.adjust(result_rows$p[ok_idx], method = p_adjust_method)
        
        result_rows$log2FC[ok_idx] <- 
            result_rows$estimate[ok_idx] / log(2)
        
        result_rows$significant[ok_idx] <- 
            result_rows$padj[ok_idx] < p_threshold
        
        result_rows$neglog10_padj[ok_idx] <- 
            -log10(result_rows$padj[ok_idx])
    } else warning("No models were successfully fitted (no 'ok' model_status).")

    out = list(
        results = result_rows,
        by_cell_type = fits_list
    )
    
    names(out$by_cell_type) <- out$results$cell_type
    
    return(out)
}