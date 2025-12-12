# Author: Shoaib Ajaib
# Date: 15/10/2025
# Usage: the script looks at the marker coverage for the GBMDeconvoluteR markers
#        across different SpatialExperiment objects created during the data 
#        processing. 

# USAGE ------------------------------------------------------------------------
# This script contains code to look at the marker coverage for the 
# GBMDeconvoluteR cell type markers and also clean up some of the cell object
# metadata.

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
io <- list(
    inputs = list(
        input_dir = "Data",
        data_dir = "Data/spatialexperiment",
        metadata = "Data/metadata"
    ),
    outputs = list(
        out_dir = "outputs/spatial_analysis/marker_comps"
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

io$inputs$spe_objects <- list(
    attempt_1_spe = "all_cosmx_spe_exclude_states", 
    attempt_2_spe = "GBMDeconv_[0-9-]+", 
    attempt_3_spe = "Ensemble"
    )

io$inputs$spe_objects <-  purrr::imap(io$inputs$spe_objects, ~{
    found <- find_file(
        file_pattern = .x, dir_path = io$inputs$data_dir, 
        overide_using_current_dir = TRUE
              )
    cli::cli_alert_info("{.y}:\t found {length(found)} file.")
    
    return(found)
    
    })

io$inputs$markers <- find_file(
    dir_path = io$inputs$metadata, 
    overide_using_current_dir = TRUE,
    file_pattern = "GBMDeconv_cell_markers_and_labels"
)

rm(find_file)

# LOAD DATA --------------------------------------------------------------------
markers <- readRDS(io$inputs$markers)

clean_cell_labels <- function(cell_labels) {
    cell_labels <- str_replace_all(cell_labels, "s$", "")
    cell_labels <- str_replace_all(cell_labels, "_", " ")
    cell_labels <- str_replace_all(cell_labels, "TAM", "Macrophage")
    return(cell_labels)
}

markers$cleaned <- lapply(markers$labels$main, clean_cell_labels)

markers$cleaned$Cancer <-  c("MES","NPC", "MES3", "PN", "MES1", "MES2","NPC1", "NPC2", "AC", "OPC")
markers$cleaned$Vasculature <-  c("Endothelial", "Mural cell", "Pericyte")
markers$cleaned$Normal <-  c("OPC normal", "Stem cell", "Ependymal", "Radial glial", "Astrocyte", "Oligodendrocyte", "Neuron")

spe_objects <- lapply(io$inputs$spe_objects, readRDS)

# The labels from the first attempt were named differently.
spe_objects$attempt_1_spe$cell_label_fct <-  factor(
    x = spe_objects$attempt_1_spe$renamed_label,
    levels = unlist(markers$cleaned, use.names = FALSE),
    ordered = TRUE
)

spe_objects$attempt_1_spe$cell_type <- spe_objects$attempt_1_spe$cell_category

table(spe_objects$attempt_1_spe$cell_label_fct)
table(spe_objects$attempt_1_spe$cell_type)

spe_objects$attempt_2_spe <-  spe_objects$attempt_2_spe[,!spe_objects$attempt_2_spe$remove_cell]

spe_objects$attempt_2_spe$cell_label_fct <-  factor(
    x = spe_objects$attempt_2_spe$cell_label,
    levels = unlist(markers$cleaned, use.names = FALSE),
    ordered = TRUE
)

spe_objects$attempt_2_spe$cell_type <- spe_objects$attempt_2_spe$cell_category

table(spe_objects$attempt_2_spe$cell_label_fct)
table(spe_objects$attempt_2_spe$cell_type)

spe_objects$attempt_3_spe$cell_label_fct <-  factor(
    x = spe_objects$attempt_3_spe$cell_label,
    levels = unlist(markers$cleaned, use.names = FALSE),
    ordered = TRUE
)

table(spe_objects$attempt_3_spe$cell_label_fct)
table(spe_objects$attempt_3_spe$cell_type)

cell_label_types_counts <- purrr::imap(unname(spe_objects), ~{
    
    out <- as.data.frame(
        table(.x$cell_label_fct)
    ) 
    
    colnames(out) <- c("cell_label", "n")
    
    out$attempt <- paste("Attempt",.y)
    
    return(out)
    
    }) %>% bind_rows()

markers$lab2cat <- setNames(
    object = unlist(markers$cleaned, use.names = FALSE),
    nm = rep(names(markers$cleaned), lengths(markers$cleaned))
)

markers$cat2lab <- setNames(
    object = rep(names(markers$cleaned), lengths(markers$cleaned)),
    nm = unlist(markers$cleaned, use.names = FALSE)
)

cell_label_types_counts$cell_type <- unname(markers$cat2lab[cell_label_types_counts$cell_label])

cell_label_types_counts$cell_type <- factor(
    cell_label_types_counts$cell_type,
    levels = names(markers$cleaned),
    ordered = TRUE
    )



svglite::svglite(
    filename = nf("cell_types_across_attempts.svg", io$outputs$out_dir),
    width = 20,
    height = 20,
)

cell_label_types_counts %>%
    ggplot(aes(x = cell_type, y = n, fill = cell_type)) +
    geom_col(position = "dodge") +
    labs(
        x = "Cell type",
        y = "Count (n)",
    ) +
    coord_flip() +
    facet_wrap(~attempt) +
    scale_fill_manual(values = spe_objects$attempt_3_spe@metadata$colors$cell_groups) +
    IMCfuncs::facetted_cell_prop_theme()

dev.off()


svglite::svglite(
    filename = nf("cell_labels_across_attempts.svg", io$outputs$out_dir),
    width = 20,
    height = 20,
)

cell_label_types_counts %>%
    ggplot(aes(x = cell_label, y = n, fill = cell_type)) +
    geom_col(position = "dodge") +
    labs(
        x = "Cell Label",
        y = "Count (n)",
    ) +
    coord_flip() +
    facet_wrap(~attempt) +
    scale_fill_manual(values = spe_objects$attempt_3_spe@metadata$colors$cell_groups) +
    IMCfuncs::facetted_cell_prop_theme()

dev.off()



# LOOK AT THE DISTRIBUTION OF CELL TYPE LABELS ---------------------------------

sort(table(markers$markers$cell_label),decreasing = TRUE)

all_6k_markers <- rownames(spe)

check_marker_coverage <- function(
        ref_markers,
        ref_cell_label_col = "cell_label",
        ref_marker_col = "marker",
        ref_type_col = "cell_type",
        lookup_markers,
        min_conserved_percent = 50){
    
    out <- list(
        marker_level = NULL,
        cell_population_level = NULL,
        type_level = NULL
    )
    
    out$marker_level <- ref_markers %>%  
        dplyr::mutate(missing_from_lookup = !(!!sym(ref_marker_col) %in% lookup_markers))
    
    
    out$cell_population_level <- out$marker_level %>%
        dplyr::group_by(!!sym(ref_cell_label_col), !!sym(ref_type_col)) %>%
        dplyr::summarise(
            `total_markers` = dplyr::n(),
            `present` = length(which(!missing_from_lookup)),
            `missing` = length(which(missing_from_lookup)), .groups = "keep"
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            present_percent = round((present / total_markers) * 100, 2),
            above_threshold = present_percent > 50
        ) %>%
        dplyr::arrange(dplyr::desc(present_percent)) %>%
        dplyr::select(
            cell_label = !!sym(ref_cell_label_col),
            cell_population =!!sym(ref_type_col), 
            total_markers, present, missing,
            present_percent, 
            above_threshold
        )
    
    out$type_level <-  out$cell_population_level %>%
        dplyr::group_by(cell_population) %>%
        dplyr::summarise(coverage = sum(missing)/sum(total_markers)*100) %>%
        dplyr::ungroup()
    
    return(out)
    
}

marker_coverage <- check_marker_coverage(
    ref_markers = markers$markers, 
    lookup_markers = all_6k_markers
    )


openxlsx::write.xlsx(
    marker_coverage,
    nf("Marker_coverage.xlsx", filepath = io$outputs$temp_out)
)

# PLOT MARKER COVERAGE ---------------------------------------------------------

marker_coverage$cell_population_level$cell_label <- factor(
    marker_coverage$cell_population_level$cell_label,
    levels = marker_coverage$cell_population_level$cell_label[
        order(marker_coverage$cell_population_level$present_percent)]
)


coverage_plot <- marker_coverage$cell_population_level %>%
    ggplot(aes(x = cell_label, y = present_percent)) +
    # Add a vertical line at 50%
    geom_hline(yintercept = 50, linetype = "dashed", color = "darkred") +
    geom_bar(stat = "identity", fill= "slateblue", alpha = 0.7) +
    coord_flip() +
    labs(
        x = "Cell Type Label",
        y = "Marker Coverage (%)",
        title = "6K CosMx GBMDeconvoluteR Marker Coverage"
    ) +
    IMCfuncs::facetted_cell_prop_theme()


svglite::svglite(
    filename = nf("Marker_coverage.svg", io$outputs$temp_out),
    width = 10,
    height = 10,
)
coverage_plot
dev.off()

# END --------------------------------------------------------------------------