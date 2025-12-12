# Author: Shoaib Ajaib
# Date: 06/11/2025
# Usage: This script visualises cell objects from SpatialExperiment data and 
#        corresponding cell masks. It generates images highlighting 
#        specific cell types within selected fields of view (FOVs).

# PACKAGES ---------------------------------------------------------------------
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(SpatialExperiment)
library(cytomapper)
library(IMCfuncs)

# I/O --------------------------------------------------------------------------
io <- list(
    inputs = list(
        input_dir = "Data/spatialexperiment",
        metadata = "Data/metadata",
        functions = list.files("Process/functions",full.names = TRUE),
        cell_masks = "Data/CellLabels",
        fov_mapping = "Data/metadata/Run105_GBM_fov_positions_file.csv"
    ),
    outputs = list(
        out_dir = "outputs/spatial_analysis/object_images"
    ),
    plots = list()
)

for (i in seq_along(io$inputs$functions)) source(io$inputs$functions[[i]], verbose = FALSE)

# Check to see if there are any processed masks
io$inputs$processed_masks <- list.files(io$inputs$cell_masks, pattern = "rds$")

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
    file_pattern = "Ensemble"
)

rm(find_file, i)

# LOAD DATA --------------------------------------------------------------------
spe <- readRDS(io$inputs$data)
spe$fov <- paste0("F", sprintf("%04d", spe$fov))

spe@metadata$labels$lab2cat <- setNames(
    object = rep(names(spe@metadata$labels$main), lengths(spe@metadata$labels$main)),
    nm = unlist(spe@metadata$labels$main, use.names = FALSE)
)

spe@metadata$colors$cells  <- splice_a_viridis(
    labels = levels(spe$cell_label_fct),
    label_to_category = spe@metadata$labels$lab2cat
    )

plot_palette_swatches(spe@metadata$colors$cells, title = "Cell Labels")

if(length(io$inputs$processed_masks) == 1){
    masks <- readRDS(file.path(io$inputs$cell_masks, io$inputs$processed_masks))
} else{
    
    fov_mapping <- readr::read_csv(
        file = io$inputs$fov_mapping,
        col_select = 1:5, skip = 1, 
        col_names = c("fov", "x_pos", "y_pos", "z_pos", "patient_surgery"), 
        show_col_types = FALSE
    ) %>%
        mutate(across(fov, ~sprintf("F%04d", .x)))
    
    masks <- cytomapper::loadImages(params$masks_path, as.is = TRUE)
    
    # Map between FOVs in the masks and those in the SPE object
    masks@elementMetadata <- S4Vectors::DataFrame(fov_mapping)
    
    stopifnot(all(masks@elementMetadata$fov %in% unique(colData(spe)$fov)))
    
    saveRDS(masks, file = nf("processed_cell_masks.rds", io$inputs$cell_masks))
}

masks@metadata$mask_length_px <- check_mask_dims(masks)
masks@metadata$mask_length_um <- masks@metadata$mask_length_px * 0.12
masks@metadata$scale_bar_lenth_um <- round(masks@metadata$mask_length_um * 0.2)

masks@metadata$scale_bar_info <- list(
    length = masks@metadata$mask_length_px * 0.2,
    label = paste0(masks@metadata$scale_bar_lenth_um, "µm"),
    position = "bottomright",
    colour = "darkred",
    margin = c(50, 50)
)

# PLOT CELL OBJECTS INITIAL RUN ------------------------------------------------
fovs <- as.data.frame(masks@elementMetadata)

# Walton64_P    c("F0022","F0029")
# Walton40_R    c("F0049", "F0056")

plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0022","F0029"),
    labels_to_highlight = c("Plasma B", "Oligodendrocyte"),
    image_titles = c("Walton64P_fov22", "Walton64P_fov29"),
    scale_bar_info = masks@metadata$scale_bar_info
        )
  
plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0049", "F0056"),
    labels_to_highlight = c("Plasma B", "Oligodendrocyte"),
    image_titles = c("Walton40R_fov49", "Walton40R_fov56"),
    scale_bar_info = masks@metadata$scale_bar_info
)

# PLOT CELL OBJECTS RUN 2 ------------------------------------------------------
fovs <- as.data.frame(masks@elementMetadata)

# Imp17_R       c("F0045")
# Imp6_R        c("F0018")
# Wal73_R       c("F0052")
# Wal73_P       c("F0007")
# Walton64_P    c("F0022","F0029")     


plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0045"),
    labels_to_highlight = c("B cell", "OPC normal"),
    image_titles = c("Imp17R_fov45"),
    scale_bar_info = masks@metadata$scale_bar_info
)

plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0018"),
    labels_to_highlight = c("B cell", "OPC normal"),
    image_titles = c("Imp6R_fov18"),
    scale_bar_info = masks@metadata$scale_bar_info
)

plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0052"),
    labels_to_highlight = c("B cell", "OPC normal"),
    image_titles = c("Wal73R_fov52"),
    scale_bar_info = masks@metadata$scale_bar_info
)


#  B cell -> OPC normal & OPC normal -> plasma B

plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0022","F0029"),
    labels_to_highlight = c("Plasma B", "B cell", "OPC normal"),
    image_titles = c("Walton64P_fov22", "Walton64P_fov29"),
    scale_bar_info = masks@metadata$scale_bar_info
)


plot_cell_objects(
    spe = spe, masks = masks, 
    mask_filt = c("F0007"),
    labels_to_highlight = c("Plasma B", "OPC normal"),
    image_titles = c("Walton73P_fov7"),
    scale_bar_info = masks@metadata$scale_bar_info
)

# END --------------------------------------------------------------------------