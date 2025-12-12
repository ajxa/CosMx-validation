# Author: Shoaib Ajaib
# Date: 19/11/2025
# Usage: This script contains the code/logic to compare the GBMDeconvoluteR scores
#        obtained for each of the datasets and to quantify the relative cell type
#        proportions obtained for each sample pair (primary -> recurrence).
# PACKAGES ---------------------------------------------------------------------
library(IMCfuncs)
library(readr)
library(dplyr)
library(stringr)
library(magrittr)
library(purrr)
library(ggplot2)
library(openxlsx)
library(patchwork)
library(viridis)
library(EnhancedVolcano)
library(tidyr)
# IO ---------------------------------------------------------------------------
io <- list(
    inputs = list(
        processed_scores = "Data/public_single_cell/Processed_GBMDeconvoluteR_scores_2025-11-19T18-03-21.rds",
        markers_and_labs = "Data/public_single_cell/GBMDeconv_cell_markers_and_labels_cleaned_2025-11-19T11-18-40.rds"
    ),
    output = list(
        comp_scores = "Outputs/single_cell_validation/comp_scores"
    ),
    plots = list(),
    args = list(),
    data = list()
)

# nd(io$output$comp_scores, recursive = TRUE, add_timestamp = FALSE)
io$output$out_temp <- nd(path = io$output$comp_scores)

# LOAD PROCESSED SCORES AND MARKERS --------------------------------------------
scores <- readRDS(io$inputs$processed_scores)
markers <- readRDS(io$inputs$markers_and_labs)

markers$labels$cell2group <- markers$markers %>%
    select(cell_label, cell_type) %>%
    distinct()

markers$labels$cell2group <- setNames(
    object = markers$labels$cell2group$cell_type,
    nm = markers$labels$cell2group$cell_label
)

io$args$colors$surgery <- plot_colours$surgery
names(io$args$colors$surgery) <- c("Primary", "Recurrent")

io$args$colors$cell_types <- setNames(
    object = c("#1AE4B6FF", "#30123BFF", "#FABA39FF", "#7A0403FF"),
    nm = names(markers$labels$main)
)

io$args$colors$diverging <- c("#446faa", "#FFFFFF", "#BB4444")

io$args$colors$positive <- c("#D1E6EF", "#ABC4DE", "#9099CA", "#8566B1", "#762D81", "#540046")


io$args$colors$cell_label_types <-  
    rep(io$args$colors$cell_types, lengths(markers$labels$main))

names(io$args$colors$cell_label_types) <- 
    unlist(markers$labels$main, use.names = FALSE) 


io$args$comps <- expand.grid(
    df = "scores",
    tbl = names(scores),
    response_var = "cell_types",
    comp_var = "surgery",
    group_var = "cell_types",
    paired_test = c(TRUE, FALSE),
    stat = "wilcox",
    p_adjust_method = "fdr",
    unnest_y_levels = list(unlist(markers$labels$main, use.names = FALSE)),
    plot_colours = list(io$args$colors$surgery),
    stringsAsFactors = FALSE
)

io$args$comps$plot_names <-  ifelse(io$args$comps$paired_test, 
    yes = paste0(io$args$comps$tbl, "_paired"),
    no =  paste0(io$args$comps$tbl, "_unpaired")
    )

io$args$comps$p_signif_col <- ifelse(io$args$comps$paired_test, "p_adj", "p")

# WILCOXON COMPARISON STATISTICS -----------------------------------------------
get_comp_stats <- function(df,
                           tbl,
                           response_var,
                           comp_var,
                           group_var,
                           unnest_y_levels,
                           stat,
                           paired_test,
                           p_adjust_method,
                           p_signif_col,
                           plot_colours,
                           plot_names,
                           ...) {
    if (is.null(tbl)) {
        df_to_use <- get(df, envir = globalenv())
    } else {
        df_to_use <- get(df, envir = globalenv())[[tbl]]
    }
    
    stats <- compare_groups(
        df = df_to_use,
        df_name = plot_names,
        response_var = response_var,
        nested_response_var = TRUE,
        unnest_levels = unnest_y_levels,
        comp_var = comp_var,
        group_var = group_var,
        facetted_plot_dims = TRUE,
        stat_y_pos_multiplier = 3,
        stat = stat,
        paired_test = paired_test,
        p_adjust_method = p_adjust_method,
        p_signif_col = p_signif_col
    )
    
    
    comp_boxplot <- comp_boxplot(
        df = df_to_use,
        x = comp_var,
        y = response_var,
        facet_by = response_var,
        unnest_y = TRUE,
        unnest_y_levels = unnest_y_levels,
        paired = paired_test,
        fill = comp_var,
        ylab = "GBMDeconvoluteR Score",
        color_palette = plot_colours,
        title = plot_names
    )
    
    
    return(
        list(stats = stats, boxplot = comp_boxplot)
    )
}

io$data$comp_stats <- pmap(io$args$comps, get_comp_stats)
names(io$data$comp_stats) <- io$args$comps$plot_names

# ADD STATS TO COMPARISONS AND SAVE BOXPLOTS -----------------------------------
add_stats_to_plot <- function(plot_stats, box_plot) {
    out <- box_plot +
        
        # Add statistics to the boxplot
        ggpubr::stat_pvalue_manual(plot_stats, label = "p_signif", tip.length = 0, size = 7) +
        
        # Add 10% spaces between the p-value labels and the plot border
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.1, 0.1))
            # ,labels = scales::percent_format()
        ) +
        # Add custom theme
        IMCfuncs::facetted_comp_bxp_theme(hide_x_labs = T)
    
    return(out)
}

io$data$comp_stats <- map(io$data$comp_stats, ~ {
    .x$boxplot <- add_stats_to_plot(plot_stats = .x$stats, box_plot = .x$boxplot)
    return(.x)
})

io$output$out_boxplots <- file.path(io$output$out_temp, "boxplots")
ndirs(io$output$out_boxplots)

save_boxplots <- function(
        .x, .y, plot_name_prefix = "", out_dir = io$output$out_boxplots
        ){
    svglite::svglite(
        filename = IMCfuncs::nf(paste0(plot_name_prefix, .y, ".svg"), out_dir),
        width = 20,
        height = 25
    )
    
    print(.x$boxplot)
    
    dev.off()
}

iwalk(io$data$comp_stats, save_boxplots)

rm(get_comp_stats, add_stats_to_plot, save_boxplots)

# FORMAT COMPARISON STATS TABLES -----------------------------------------------
add_Log2FC <- function(comparison, tbl, comp_data_list = scores) {
    comparison$stats <- IMCfuncs::calculate_Log2FC(
        score_df = comp_data_list[[tbl]],
        from = "Primary",
        to = "Recurrent",
        scores = "cell_types",
        add_to_stats_df = TRUE,
        stats_df = comparison$stats
    )
    
    return(comparison)
}

io$data$comp_stats <- map2(
    .x = io$data$comp_stats, .y = io$args$comps[["tbl"]], ~add_Log2FC(comparison = .x, tbl = .y)
    )

format_stats_tbl <- function(comparisons, unique_cell_types=NULL) {

    if(is.null(unique_cell_types)) stop("Please provide unique_cell_types vector")
    
    out <- comparisons$stats %>%
        select(cell_types, dataset, p, p_adj, p_signif, log2FC) %>%
        mutate(across(dataset, ~ str_remove_all(., "\\s+\\(.*\\)") %>% tolower())) %>%
        # mutate(across(c(p, p_adj), ~sprintf("%.2e", .))) %>%
        # mutate(across(log2FC, ~sprintf("%.3f", .))) %>%
        mutate(change = ifelse(log2FC > 0, "up", "down")) %>%
        mutate(across(cell_types, ~ factor(.x, levels = unique_cell_types))) %>%
        arrange(cell_types)
    
    out <- out %>%
        tidyr::pivot_wider(
            names_from = "dataset",
            values_from = c(p, p_adj, p_signif, log2FC, change),
            names_glue = "{dataset}_{.value}"
        )
    
    return(out)
}

io$data$cleaned_stats <- map(io$data$comp_stats, ~{
    format_stats_tbl(
        comparisons = .x, 
        unique_cell_types = unlist(markers$labels$main, use.names = FALSE)
        )
})

# SAVE STATISTICS TABLES AS EXCEL WORKBOOK -------------------------------------
io$output$out_stats <- file.path(io$output$out_temp, "stats")
ndirs(io$output$out_stats)

create_formatted_workbook <- function(workbook_tbls, outname, outpath) {
    wb <- openxlsx::createWorkbook()
    
    # Workbook Styling ----
    openxlsx::modifyBaseFont(wb,
                             fontSize = 14,
                             fontName = "Arial",
                             fontColour = "black"
    )
    
    
    headerstyle <- openxlsx::createStyle(
        fontColour = "#ffffff",
        fgFill = "#003087",
        halign = "left",
        valign = "top",
        textDecoration = "bold"
    )
    
    
    normal_text_style <- openxlsx::createStyle(halign = "left", valign = "top")
    
    # Add worksheet tabs ----
    
    add_df_to_workbook <- function(sheet_df, tab_name, workbook) {
        # Add the sheet to the workbook
        openxlsx::addWorksheet(
            wb,
            sheetName = tab_name,
            gridLines = FALSE,
            zoom = 100
        )
        
        openxlsx::writeData(
            wb,
            sheet = tab_name,
            x = sheet_df,
            startCol = 1,
            startRow = 1,
            withFilter = FALSE,
            borders = "all",
            borderColour = "darkgrey",
            borderStyle = "hair",
            colNames = TRUE,
            rowNames = FALSE,
            headerStyle = headerstyle
        )
        
        # Format the table
        openxlsx::setColWidths(wb, tab_name, cols = 1:ncol(sheet_df), widths = nchar(names(sheet_df)) * 1)
        
        openxlsx::addStyle(
            wb,
            sheet = tab_name,
            style = normal_text_style,
            rows = 1:(nrow(sheet_df) * 1.1),
            cols = ncol(sheet_df),
            gridExpand = TRUE,
            stack = TRUE
        )
    }
    
    purrr::imap(workbook_tbls, ~ add_df_to_workbook(
        sheet_df = .x,
        tab_name = .y,
        workbook = wb
    ))
    
    # Save the created workbook ----
    
    out_filename <- nf(glue::glue("{outname}.xlsx"), outpath) 
    
    openxlsx::saveWorkbook(
        wb = wb,
        file = out_filename,
        overwrite = TRUE
    )
    
    cli::cli_alert_success("Saved formatted workbook: {out_filename}")
    
}

create_formatted_workbook(
    workbook_tbls = io$data$cleaned_stats,
    outname = "GBMDeconvoluteR_score_comparisons",
    outpath = io$output$out_stats
)

rm(add_Log2FC, create_formatted_workbook, format_stats_tbl)

# MIXED-MODEL USING ALL DATA ---------------------------------------------------
convert_to_long_df <- 
    function(
        df, 
        celltype_col = "cell_types",
        unique_sample_col = "Mixture",
        clean=TRUE,
        cleaned_celltype_order = unlist(markers$labels$main, use.names = FALSE)
        ){
    
    if(!celltype_col %in% colnames(df)) stop("celltype_col not found in df supplied")
    
    unique_cell_types <- unique(names(df[[celltype_col]][[1]]))
    
    if(!all(unique_cell_types %in% cleaned_celltype_order)){
        stop("Not all cleaned_celltype_order found in celltype_col of df supplied")
    }
    
    pivot_long_cols <- colnames(df)[colnames(df) != celltype_col]
    
    if(length(pivot_long_cols) == 0) stop("No columns to pivot longer found in df supplied apart from celltype_col")
    
    long_df <- df %>%
        tidyr::unnest_wider(col = celltype_col, names_repair = "unique") %>%
        tidyr::pivot_longer(
            cols = !matches(paste0("^(", paste0(pivot_long_cols, collapse = "|"), ")$")),
            names_to = "cell_type",
            values_to = "score"
        )
    
    # sense check (should be TRUE):
    predicted_nrow <- nrow(long_df) == (length(unique(df[[unique_sample_col]])) * length(unique_cell_types ))
    
    if(!predicted_nrow) {
        stop("Number of rows in long_df does not match expected number of rows")
    }else {
        cli::cli_alert_success("Successfully converted to long dataframe with {nrow(long_df)} rows")
    }
    
    if(clean){
     
    long_df$surgery <- factor(long_df$surgery, levels = c("Primary", "Recurrent"))
    long_df$tissue_source <- factor(long_df$tissue_source)
    long_df$patient <- factor(long_df$patient)
    long_df$score <- as.numeric(long_df$score)
    
    long_df$cell_type <- factor(
        x = long_df$cell_type, 
        levels = cleaned_celltype_order,
        ordered = TRUE
        )
    
    return(long_df)
        
    }else return(long_df)
    
} 

scores$all_long <- convert_to_long_df(scores$all)
scores$nomura_long <- convert_to_long_df(scores$nomura)

library(lme4)
library(lmerTest)

source("Process/functions/mixed_models.R")

io$data$mixed_model_res <- list()

io$data$mixed_model_res$all <- run_celltype_mixed_models(
    long_df = scores$all_long,
    celltypes = as.character(unique(scores$all_long$cell_type)),
    p_threshold = 0.05, 
    p_adjust_method = "fdr" 
    ) 

io$data$mixed_model_res$nomura <- run_celltype_mixed_models(
    long_df = scores$nomura_long,
    celltypes = as.character(unique(scores$nomura_long$cell_type)),
    p_threshold = 0.05, 
    p_adjust_method = "fdr" 
) 

rm(fit_celltype_model, run_celltype_mixed_models)

# ENHANCED VOLCANO PLOTS -------------------------------------------------------
source("Process/functions/enhanced_volcano_plots.R")

save_svg_plot <- function(plot, filename, prefix = NULL, outdir, width=15, height=15){
    
   if(is.null(prefix)){
       
       out_filename <- IMCfuncs::nf(
           filename = paste0(filename, ".svg"),
           filepath =  outdir
           )
       
   }else{

       out_filename <- IMCfuncs::nf(
           filename = paste(prefix,paste0(filename, ".svg"), sep = "_"), 
           filepath = outdir
       )
       
   }
    
    svglite::svglite(
        filename = out_filename,
        width = width,
        height = height
    )
    print(plot)
    dev.off()
    
    cli::cli_alert_success("\nSaved plot: {out_filename}")
}

io$output$out_volcano <- file.path(io$output$out_temp, "volcano_plots")
ndirs(io$output$out_volcano)

# Wilcoxon test volcano plots:

io$plots$wilcox_prep <- imap(io$data$comp_stats, ~{
    
    out <- prep_volcano_data(
        stats_tbl = .x$stats,
        current_cells = io$args$colors$cell_label_types,
        method_info = list(
            method = unique(.x$stats[["method"]]),
            `paired test` = unique(.x$stats[["paired_test"]]),
            `p adjust method` = unique(.x$stats[["p_adj_method"]])
            ), 
        plot_title = str_to_title(str_replace(.y, "_[a-z]+$", ""))
        )
    
    return(out)
    })

io$plots$wilcox <- map(io$plots$wilcox_prep, plot_volcano)

iwalk(io$plots$wilcox, save_svg_plot, prefix = "Wilcoxon", outdir = io$output$out_volcano)

# Mixed-model volcano plots:
io$plots$mixed_model_prep <- list()

io$plots$mixed_model_prep$all <- prep_volcano_data(
    stats_tbl = io$data$mixed_model_res$all$results,
    current_cells = io$args$colors$cell_label_types,
    cell_type_col = "cell_type",
    p_value_col = "padj",
    plot_title = "All Datasets (Wang & Nomura)",
    method_info = list(
        method = "Mixed-effects model",
        formula = "score ~ surgery + tissue_source + (1 | patient)",
        `confounding batch covariate` = "tissue_source",
        `p adjust method` = "fdr",
        `n patients (pairs)` = length(unique(scores$all_long$patient))
    )
)

io$plots$mixed_model_prep$nomura <- prep_volcano_data(
    stats_tbl = io$data$mixed_model_res$nomura$results,
    current_cells = io$args$colors$cell_label_types,
    cell_type_col = "cell_type",
    p_value_col = "padj",
    plot_title = "Nomura et al. Dataset",
    method_info = list(
        method = "Mixed-effects model",
        formula = "score ~ surgery + tissue_source + (1 | patient)",
        `confounding batch covariate` = "tissue_source",
        `p adjust method` = "fdr",
        `n patients (pairs)` = length(unique(scores$nomura_long$patient))
    )
)

io$plots$mixed_model <- map(io$plots$mixed_model_prep, plot_volcano)

iwalk(
    io$plots$mixed_model, 
    .f = save_svg_plot,  
    prefix = "Mixed-model", 
    outdir = io$output$out_volcano
    )
      
      
openxlsx::write.xlsx(
    x = lapply(io$data$mixed_model_res, `[[`, "results") ,
    file = IMCfuncs::nf("mixed_model_results.xlsx", io$output$out_stats),
    overwrite = TRUE
)

rm(prep_volcano_data, plot_volcano, set_axis_lims)

# DELTA SCORE CORRELATIONS -----------------------------------------------------
scores$all_delta_long <-
    scores$all_long %>%
    group_by(patient, cell_type) %>%
    summarise(
        tissue_source = first(tissue_source),
        Primary   = score[surgery == "Primary"],
        Recurrent = score[surgery == "Recurrent"],
        .groups = "drop"
    ) %>%
    mutate(delta_score = Recurrent - Primary)

scores$all_delta_wide <-
    scores$all_delta_long %>%
    select(patient, cell_type, delta_score) %>%
    pivot_wider(names_from = cell_type, values_from = delta_score)

scores$all_delta_mat <- as.matrix(scores$all_delta_wide[, -1])
rownames(scores$all_delta_mat) <- scores$all_delta_wide$patient

# Adjusting the cor matrix to account for the tissue_source batch:
scores$all_batch_df <- scores$all_delta_long %>% distinct(patient, tissue_source)

scores$all_batch_vec <- 
    scores$all_batch_df$tissue_source[match(rownames(scores$all_delta_mat), scores$all_batch_df$patient)]

# this will return the Δ scores that have had batch (tissue_source) variation removed
io$data$delta_score_batch_resid_mat <- 
    apply(scores$all_delta_mat, 2, function(y) {
    fit <- lm(y ~ scores$all_batch_vec)
    resid(fit)
    })

source("Process/functions/correlations_galore.R")

io$data$cor_resid_mat <- cor(
    x = io$data$delta_score_batch_resid_mat, 
    use = "pairwise.complete.obs", 
    method = "spearman"
    )

io$data$cor_resid_pval_mat <- cor_mtest(
    mat = io$data$cor_resid_mat, method = "spearman", exact = FALSE
    )

io$data$cor_resid_adj_pval_mat <- adjust_p_matrix(
    p_mat = io$data$cor_resid_pval_mat,
    method = "BH", 
    diag_value = 1
    )

io$data$corr_plot_data <- clean_corr_plot_data(
    corr_matrix = io$data$cor_resid_mat,
    p_matrix = io$data$cor_resid_adj_pval_mat,
    remove_diagonal = FALSE, 
    p_sig_level = 0.05,
    minimum_cor = 0.5
)

# PLOT CORRELATION MATRIX ------------------------------------------------------
io$plots$labs <- list()

io$plots$labs$plot_subtitle <- glue::glue(
    "Heatmap of Spearman correlation coefficients using Δ(rec-prim) GBMDeconvoluteR cell-type scores.",
    "\nThe Δ scores were corrected (linear model) to exclude the 'tissue_source' batch effect.",
    "\nColoured tiles represent the strength and direction of correlation coefficients (abs(min_r)>=0.5).",
    "\nOnly statistically significant comparisons are shown (fdr, adjusterd p < 0.05).",
    "\nParied patients: Nomura et.al (n=52), Wang et.al (n=28)",
    "\n\n\n"
)

# To show the diagonals:
io$data$corr_plot_data$plot_value[which(io$data$corr_plot_data$Var1 == io$data$corr_plot_data$Var2)] = 1

p <- plot_corr(
    corr_plot_data = io$data$corr_plot_data,
    title = "Nomura & Wang Pseudobulks: Δ GBMDeconvoluteR Score Correlations",
    subtitle = io$plots$labs$plot_subtitle,
    diverging_colours = io$args$colors$diverging,
    pos_colours = io$args$colors$positive,
    show_labels = FALSE
)

out <- add_highlight_regions(
    baseplot = p,
    corr_matrix = io$data$cor_resid_mat,
    groups = unname(markers$labels$cell2group[levels(io$data$corr_plot_data$Var1)]),
    highlight_colors = io$args$colors$cell_types
)


save_svg_plot(
    plot = out,
    filename = "Delta_score_heatmap",
    outdir = io$output$out_temp,
    width = 17,
    height = 17
)

# SAVE OUTPUTS -----------------------------------------------------------------
io$data$inputs <- scores
io$data$markers <- markers

saveRDS(
    object = io, 
    file = nf("run_data.rds", io$output$out_temp)
    )

# END --------------------------------------------------------------------------