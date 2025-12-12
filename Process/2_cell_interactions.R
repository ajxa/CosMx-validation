# Author: Shoaib Ajaib
# Date: 15/10/2025
# Usage: this script is designed to perform spatial analysis on the data 
#        generated in the previous steps. It also comprises code to look at the
#        cell-cell interactions using different spatial graph constructions and
#        visualise the results. 
# USAGE ------------------------------------------------------------------------
# This script is designed to perform spatial analysis on the data generated in
# the previous steps. The following types of analysis will be performed:
#
# 1. Create spatial interaction graphs
# 3. Determine and test cellular interactions

# OPTIONS ----------------------------------------------------------------------
# options(scipen = 999)

# PACKAGES ---------------------------------------------------------------------
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(SpatialExperiment)
library(IMCfuncs)
library(imcRtools)

# I/O --------------------------------------------------------------------------
attempt <- 3 # change accordingly

io <- list(
    inputs = list(
        input_dir = "Data/spatialexperiment",
        metadata = "Data/metadata"
    ),
    outputs = list(
        out_dir = file.path("outputs/spatial_analysis", paste0("attempt_",attempt))
    ),
    plots = list()
)

# create out directory
ndirs(io$outputs)

# create time-stamped output directory
io$outputs$temp_out <- nd(path = io$outputs$out_dir)

# obtain the most-recent data from a time-stamped directory
find_file <- function(dir_path,
                      file_pattern,
                      dir_pattern = "^[T0-9-]+$",
                      dir_select = c("recent", "oldest"),
                      overide_using_current_dir = FALSE) {
    if (missing(dir_path)) stop("No directory path provided")
    if (missing(file_pattern)) stop("No file pattern provided")
    
    if( overide_using_current_dir){
        dirs_found <- list.files(dir_path, pattern = file_pattern, full.names = TRUE)
        return(dirs_found)
    }
    
    dirs_found <- list.files(dir_path, pattern = dir_pattern)
    
    if (length(dirs_found) == 0) stop("No directories found")
    
    dir_selected <- switch(match.arg(dir_select, several.ok = FALSE),
                           recent = dirs_found[order(dirs_found, decreasing = TRUE)][[1]],
                           oldest = dirs_found[order(dirs_found, decreasing = FALSE)][[1]]
    )
    
    file_found <- list.files(
        path = file.path(dir_path, dir_selected),
        pattern = file_pattern,
        ignore.case = TRUE,
        recursive = FALSE,
        full.names = TRUE
    )
    
    if (length(file_found) == 0) {
        stop("No files found")
    } else if (length(file_found) > 1) {
        stop("Multiple files found")
    } else {
        return(file_found)
    }
}

io$inputs$data <- find_file(
    dir_path = io$inputs$input_dir, 
    overide_using_current_dir = TRUE,
    # file_pattern = "GBMDeconv",
    file_pattern = "Ensemble",
    
)

rm(find_file)

# LOAD DATA --------------------------------------------------------------------
spe <- readRDS(io$inputs$data)

# CLEAN UP LABELS AND SUBSET CELLS ---------------------------------------------
spe$patient_surgery_fov <-  paste(spe$patient_surgery, spe$fov,sep = "_")

length(unique(spe$patient_surgery_fov))
length(unique(spe$patient_surgery))

# Remove cell labels which are not used in the paper2 analysis
spe <- spe[,!spe$remove_cell]
spe
# LOOK AT THE DISTRIBUTION OF CELL TYPE LABELS ---------------------------------
cell_label_distributions <- set_names(
    x = c("patient", "patient_surgery", "patient_surgery_fov"),
    nm  = c("patient", "patient_surgery", "patient_surgery_fov")
) %>%
    map(~{
        as.data.frame(table(spe[[.x]], spe$cell_label)) %>%
            tidyr::pivot_wider(
                names_from = Var2,
                values_from = Freq
            ) %>%
            dplyr::rename(
                sample_id = Var1
            )
        
        
    })

openxlsx::write.xlsx(
    cell_label_distributions,
    nf("Cell_label_counts.xlsx", filepath = io$outputs$temp_out)
)

# PIXEL EXPANSION-BASED SPATIAL GRAPHS -----------------------------------------
spe <- imcRtools::buildSpatialGraph(
        object = spe, 
        img_id = "patient_surgery_fov", # Build graphs for each FOV
        type = "expansion", 
        threshold = 27.5, # This is in microns and based on Greenwald et.al 2021
        coords = c("x_centroid_um", "y_centroid_um"),
        name = "expansion_graph"
        )

spe_expansion_interact <- testInteractions(
    object = spe,
    group_by = "patient_surgery",
    label = "cell_label_fct",
    colPairName = "expansion_graph",
    method = "histocat",
    iter = 1000,
    p_threshold = 0.05,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
)

spe_expansion_ia_df <- spe_expansion_interact %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        count_method = "histocat",
        patient = stringr::str_extract(group_by, "^(?i)[a-z0-9]+")
    )

# DELAUNAY TRIANGULATION-BASED SPATIAL GRAPHS ----------------------------------
spe <- imcRtools::buildSpatialGraph(
    object = spe, 
    img_id = "patient_surgery_fov", # Build graphs for each FOV
    type = "delaunay", 
    threshold = 50, # This is in microns
    coords = c("x_centroid_um", "y_centroid_um"),
    name = "delaunay_50"
)

spe_delaunay_interact <- testInteractions(
    object = spe,
    group_by = "patient_surgery",
    label = "cell_label_fct",
    colPairName = "delaunay_50",
    method = "histocat",
    iter = 1000,
    p_threshold = 0.05,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
)

spe_delaunay_ia_df <- spe_delaunay_interact %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
        count_method = "histocat",
        patient = stringr::str_extract(group_by, "^(?i)[a-z0-9]+")
    )

# SAVE INTERACTION RESULTS  ----------------------------------------------------
openxlsx::write.xlsx(
    x = list(
            expansion_27_5 = spe_expansion_ia_df,
            delaunay_50 = spe_delaunay_ia_df
            ),
    nf("Spatial_interaction_results.xlsx", filepath = io$outputs$temp_out)
)

# PLOTTING CELL INTERACTIONS FUNCTIONS -----------------------------------------
clean_ia_data <- function(ia_df,
                          filter_by = NA,
                          filter_val = NA,
                          filter_min_patients = TRUE,
                          min_patients = 3,
                          count_method = c("histocat", "classic", "patch"),
                          graph_type = "expansion (threshold = 27.5)",
                          group_by = "patient_surgery",
                          p_threshold = 0.01,
                          label_order = levels(lab_spe$manual_gating)) {
    if (length(filter_val) > 1) stop("Only one filter value can be selected")
    
    filter_on <- if (is.na(filter_by)) "none" else filter_by
    filter_on_vals <- if (is.na(filter_val)) "none" else filter_val
    count_filter <- match.arg(count_method, several.ok = FALSE)
    
    if (filter_on != "none") {
        plot_data <- ia_df %>%
            dplyr::filter(!!sym(filter_on) %in% filter_on_vals) %>%
            dplyr::filter(count_method == count_filter)
        
        if (nrow(plot_data) == 0) stop("No data found for the given filter value")
    } else {
        plot_data <- ia_df %>%
            dplyr::filter(count_method == count_filter)
    }
    
    plot_data <- plot_data %>%
        dplyr::group_by(from_label, to_label) %>%
        dplyr::summarize(
            sum_sigval = sum(sigval, na.rm = TRUE),
            pct_sigval = abs(sum_sigval) / n() * 100,
            min_pval = min(p, na.rm = TRUE),
            interact_type = ifelse(sum_sigval == 0, NA, ifelse(sum_sigval > 0, "Interacting", "Avoiding")),
            .groups = "drop"
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            filter_on = filter_on,
            filter_vals = paste0(filter_on_vals, collapse = ","),
            count_method = count_filter,
            min_patient_filtered = filter_min_patients,
            min_patients = ifelse(filter_min_patients, min_patients, NA),
            cell_state = ifelse(filter_on %in% "state_level", filter_on_vals, NA)
        ) %>%
        dplyr::mutate(
            across(interact_type, ~ factor(., levels = c("Interacting", "Avoiding"))),
            across(c("from_label", "to_label"), ~ factor(., levels = label_order))
        )
    
    if (filter_min_patients) plot_data$interact_type[abs(plot_data$sum_sigval) < min_patients] <- NA
    
    # Create a list of labels for the plot:
    title_suffix <- "All Samples"
    
    label_list <- list(
        title = glue::glue("Cell-Cell Interactions - {title_suffix}"),
        subtitle = glue::glue(
            "graph: {graph_type}",
            "count method: {unique(ia_df$count_method)}",
            "minimum {group_by}: {min_patients}",
            "p threshold: < {p_threshold}",
            .sep = "\n"
        ),
        caption = "\n*only significant interactions are displayed"
    )
    
    return(
        list(
            plot_df = plot_data,
            plot_labs = label_list
        )
    )
}

highlight_colours <-  spe@metadata$colors$cell_groups[
    c("Immune", "Cancer", "Normal", "Vasculature")
]


cell_groups <- as.data.frame(colData(spe)) %>% select(cell_label, cell_type) %>% distinct()
cell_groups <- table(cell_groups$cell_type)[names(highlight_colours)]
cell_groups <- rep(names(cell_groups), cell_groups)

add_highlight_regions <- function(baseplot,
                                  groups = cell_groups,
                                  highlight_colors = "black",
                                  highlight_width = 3,
                                  legend_key_size = 6) {
    unique_groups <- unique(groups)
    
    highlight_regions <- data.frame(
        group = factor(unique_groups, levels = names(highlight_colors)),
        xmin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        xmax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        ymin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        ymax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        color = highlight_colors[unique(groups)]
    )
    
    if (length(highlight_colors) == 1) {
        baseplot +
            geom_rect(
                data = highlight_regions,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                color = highlight_colors, fill = NA,
                linewidth = highlight_width,
                inherit.aes = FALSE
            )
    } else {
        baseplot +
            geom_rect(
                data = highlight_regions,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
                fill = NA, linewidth = highlight_width, inherit.aes = FALSE
            ) +
            scale_color_manual(
                name = "",
                values = highlight_colors
            ) +
            guides(
                color = guide_legend(override.aes = list(size = legend_key_size)),
            )
    }
}

plot_ia <- function(ia_clean_list,
                    point_size = 18,
                    highlight_box_colors = highlight_colours,
                    plot_title = "CosMx cell-cell interactions") {
    
    p <- ia_clean_list$plot_df %>%
        ggplot(
            aes(x = from_label, y = to_label)
        ) +
        geom_tile(fill = "#fffdfa", color = "grey50", linewidth = 0.1, linetype = 2) +
        geom_point(
            data = subset(ia_clean_list$plot_df, !is.na(interact_type)),
            aes(fill = interact_type, shape = interact_type),
            color = "black", size = point_size, stroke = 0.5, na.rm = TRUE
        ) +
        scale_shape_manual(
            name = "",
            values = c(
                "Interacting" = 21,
                "Avoiding" = 22
            )
        ) +
        scale_fill_manual(
            name = "",
            values = c(
                "Interacting" = "darkgreen",
                "Avoiding" = "darkred"
            ),
            na.value = "white"
        ) +
        ggplot2::labs(
            title =  plot_title,
            subtitle = ia_clean_list$plot_labs$subtitle,
            caption = ia_clean_list$plot_labs$caption
        ) +
        ggplot2::xlab("from cell type ...") +
        ggplot2::ylab("to cell type ...") +
        theme_minimal(base_size = 16) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(size = 16, angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(size = 16, face = "bold"),
            axis.title = element_text(size = 20, face = "bold", colour = "grey50"),
            plot.title = element_text(size = 24, face = "bold"),
            plot.subtitle = element_text(size = 18, face = "italic"),
            plot.caption = element_text(size = 14, face = "italic")
        ) 
    
    p <- add_highlight_regions(baseplot = p, highlight_colors = highlight_box_colors)
    
    return(p)
}

# PLOTTING ALL EXPANSION GRAPH CELL INTERACTIONS -------------------------------
ia_expansion_cleaned <- clean_ia_data(
    ia_df = spe_expansion_ia_df,
    min_patients = 2,
    group_by = "FOV's", 
    p_threshold = 0.05,
    label_order = unlist(spe@metadata$labels$main,use.names = FALSE)
                            )

io$plots$ia_expansion <- plot_ia(
    ia_clean_list = ia_expansion_cleaned, 
    highlight_box_colors = highlight_colours,
    plot_title = "All cell label interactions"
    )

svglite::svglite(
    filename = nf("expansion_min2_interacts.svg", io$outputs$temp_out),
    width = 20,
    height = 20,
)
print(io$plots$ia_expansion)
dev.off()

# PLOTTING ALL DELAUNAY TRIANGULATION INTERACTIONS -----------------------------
ia_delaunay_cleaned <- clean_ia_data(
    ia_df = spe_delaunay_ia_df,
    min_patients = 2,
    group_by = "FOV's",
    graph_type = "Delaunay (max_dist=50)",  
    p_threshold = 0.05,
    label_order = unlist(spe@metadata$labels$main,use.names = FALSE)
)

io$plots$ia_delaunay <- plot_ia(
    ia_clean_list = ia_delaunay_cleaned, 
    highlight_box_colors = highlight_colours,
    plot_title = "All cell label interactions"
)

svglite::svglite(
    filename = nf("delaunay_min2_interacts.svg", io$outputs$temp_out),
    width = 20,
    height = 20,
)
print(io$plots$ia_delaunay)
dev.off()

# END --------------------------------------------------------------------------