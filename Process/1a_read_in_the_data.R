# Author: Shoaib Ajaib
# Date: 15/10/2025
# Usage: this script takes the Seurat class object(s) and converts them into 
#        SpatialExperiment objects which inherits from the SingleCellExperiment
#        class and is designed to represent spatially resolved transcriptomics data

# PACKAGES ---------------------------------------------------------------------
library(imcRtools)
library(SpatialExperiment)
library(Seurat)
library(cytomapper)
library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(viridis)
library(purrr)
library(IMCfuncs)


# I/O --------------------------------------------------------------------------
io <- list(
    
    inputs = list(
        functions = "process/functions/splice_a_viridis.R",
        data = "data/seurat",
        metadata = "data/metadata",
        to_load = list()
        ),
    
    output = list(
        raw_data_out = "data/spatialexperiment"
        )
    
)

source(io$inputs$functions)

if(!dir.exists(io$output$raw_data_out)) dir.create(io$output$raw_data_out)

io$inputs$to_load$seurat_object <- list.files(
    io$inputs$data, 
    # pattern = "^GBMDeconv", 
    pattern = "Ensemble", 
    full.names = TRUE
    )

io$inputs$to_load$cosmx_6k_metadata <- list.files(
    io$inputs$metadata, pattern = "Run105_GBM_metadata", full.names = TRUE
    )

# io$inputs$to_load$cosmx_6k_fovs <- list.files(
#     io$inputs$metadata, pattern = "Run105_GBM_fov",full.names = TRUE
#     )

io$inputs$to_load$cell_markers <- list.files(
    io$inputs$metadata, pattern = "gbmdeconv_cell_markers.+\\.rds$",
    ignore.case =  TRUE, full.names = TRUE
    )

# READ IN THE DATA -------------------------------------------------------------
if(all(lengths(io$inputs$to_load) == 1)){
    
    seurat_obj <- readRDS(io$inputs$to_load$seurat_object)
    metadata_6k <- readr::read_csv(io$inputs$to_load$cosmx_6k_metadata, show_col_types = FALSE)
    cell_markers <- readRDS(io$inputs$to_load$cell_markers)
    # fovs_6k <- readr::read_csv(io$inputs$to_load$cosmx_6k_fovs, show_col_types = FALSE)
 
} else stop("Some input files are missing or there are multiple files found.")
    
# CLEAN UP AND EXTRACT THE SAMLES/CELLS ----------------------------------------
seurat_metadata <-  seurat_obj@meta.data

patients <- list(
    all = seurat_metadata %>% 
        group_by(Patient_ID, P_R) %>%
        distinct(Patient_ID) %>%
        pull()
        )

patients$paired <-  patients$all[duplicated(patients$all)]
patients$singletons <- setdiff(patients$all, patients$paired)


print_fun <- function() {
    cli::cli_alert_info(c("Patients:"))
    cli::cli_ul("Paired = {length(patients$paired)} ({length (patients$paired)*2} total)")
    cli::cli_ul("Single = {length(patients$singletons)}")
    
    cli::cat_line()
    
    cli::cli_alert_info(c("Feilds of View:"))
    cli::cli_ul("Total = {length(unique(metadata_6k$fov))}")
   
}

print_fun()

# CLEAN UP SEURAT OBJECT METADATA ----------------------------------------------
clean_cell_labels <- function(cell_labels) {
    cell_labels <- str_replace_all(cell_labels, "s$", "")
    cell_labels <- str_replace_all(cell_labels, "_", " ")
    cell_labels <- str_replace_all(cell_labels, "TAM", "Macrophage")
    return(cell_labels)
}

cell_markers$labels$fine <-  lapply(cell_markers$labels$fine, \(x) return(clean_cell_labels(x)))
cell_markers$labels$main <-  lapply(cell_markers$labels$main, \(x) return(clean_cell_labels(x)))

cell_markers$filt_markers <- unlist(cell_markers$labels$main, use.names = FALSE)    

cell_markers$labels$lab2cat <- setNames(
    object = rep(names(cell_markers$labels$main), lengths(cell_markers$labels$main)),
    nm = unlist(cell_markers$labels$main, use.names = FALSE)
    )

cell_markers$labels$lab2cat
    
seurat_metadata$cell_type <- unname(cell_markers$labels$lab2cat[seurat_metadata$Resolved_Label])

cleaned_metadata <- seurat_metadata %>%
    dplyr::select(
        patient = Patient_ID,
        surgery = P_R,
        fov,
        cell_ID,
        object_id = id,
        Area,
        width_um = Width,   
        height_um = Height, 
        # cell_label = leiden_cell,
        # cell_type = leiden_type,
        cell_label = Resolved_Label,
        cell_type
    ) %>%
    mutate(
        across(surgery, ~ifelse(.x == "P", "Prim", "Rec")),
        across(cell_ID, ~{as.numeric(str_split(.x, "_",simplify = T)[,4])}),
        across(width_um, ~.x * 0.12),
        across(height_um, ~.x * 0.12)
        # across(cell_label, as.character),
        # across(cell_label , ~clean_cell_labels(.x))
    ) %>%
    mutate(
        patient_surgery = paste(patient, surgery, sep = "_"),
        paired_patient = ifelse(patient %in% patients$paired, TRUE, FALSE), 
        .after = surgery
    )

cleaned_metadata <- metadata_6k %>%
    mutate(
        x_centroid_um = CenterX_local_px * 0.12, 
        y_centroid_um = CenterY_local_px * 0.12 
    ) %>%
    select(
        fov, 
        cell_ID, 
        x_centroid_um, 
        y_centroid_um 
    ) %>%
    left_join(cleaned_metadata, ., by = c("fov", "cell_ID"))

filt_exprs_matrix <- as.matrix(seurat_obj@assays$Nanostring$data)


if(all(colnames(filt_exprs_matrix) == cleaned_metadata$object_id)) {
        cli::cli_alert_success(("Very good, carry on."))
    } else cli::cli_alert_danger(("Something has gone terribly wrong."))


# CLEAN UP LABELS --------------------------------------------------------------
cleaned_metadata$remove_cell <- !cleaned_metadata$cell_label %in% cell_markers$filt_markers

# cleaned_metadata$cell_category <- ifelse(
#     cleaned_metadata$cell_label %in% cell_markers$labels$main$Immune, "Immune",
#     ifelse(
#         cleaned_metadata$cell_label %in% cell_markers$labels$main$Cancer, "Cancer",
#         ifelse(
#             cleaned_metadata$cell_label %in% cell_markers$labels$main$Normal, "Normal",
#             ifelse(
#                 cleaned_metadata$cell_label %in% cell_markers$labels$main$Vasculature, "Vasculature",
#                 NA
#                  
#             )
#         )
#     )
# )

sort(table(cleaned_metadata$cell_type))
sort(table(cleaned_metadata$cell_label))
table(cleaned_metadata$remove_cell)
# sort(table(cleaned_metadata$cell_label[which(!cleaned_metadata$remove_cell)]))

# COLOURS -----------------------------------------------------------
io$inputs$colors <- list(
    dataset_pheno = c(
        Prim = "#2166AC", Rec = "#B2182B", up = "#004529", down = "#FF7F00"
    ),
    cell_groups = c(
        Immune = "#1AE4B6FF", 
        Cancer = "#30123BFF", 
        Normal = "#FABA39FF", 
        Vasculature = "#7A0403FF"
    ),
    diverging = c(
        "#446faa", "#FFFFFF", "#BB4444"
    ),
    positive = c(
        "#D1E6EF", "#ABC4DE", "#9099CA", "#8566B1", "#762D81", "#540046"
    ),
    imc_visual = c(
        "#0000FF", "#00FF00", "#FF0000", "#FF00FF", "#00FFFF", "#FFFF00", "#FFFFFF", "#FFA500"
    )
)


io$inputs$colors$cells <- splice_a_viridis(
    labels = cell_markers$filt_markers,
    label_to_category = cell_markers$labels$lab2cat
    ) 

io$inputs$colors$cells$fills <- c(io$inputs$colors$cells$fills, "Unknown" = "#BDBDBD")
io$inputs$colors$cells$outlines <- c(io$inputs$colors$cells$outlines, "Unknown" = "black")

plot_palette_swatches(io$inputs$colors$cells, title = "Cell Labels")


# CREATE SPATIAL EXPERIMENT OBJECTS ---------------------------------------------
spe <- SpatialExperiment(
    assay = list(data = filt_exprs_matrix),
    colData = cleaned_metadata,
    spatialCoordsNames = c("x_centroid_um", "y_centroid_um")
)

spe@metadata <- list(
    colors = io$inputs$colors,
    labels = cell_markers$labels
    )

spe$cell_label_fct <- factor(
    x = spe$cell_label,
    levels = cell_markers$filt_markers,
    ordered = TRUE
    )

# Remove cells with missing/unknown cell labels
spe <- spe[,!spe$remove_cell]

# SAVE OBJECTS -----------------------------------------------------------------
saveRDS(spe, nf("GBMDeconv_Ensemble.rds", filepath = io$output$raw_data_out)) 

# END --------------------------------------------------------------------------