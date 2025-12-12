# Author: Shoaib Ajaib
# Date: 13/11/2025
# Usage: This script uses the single-cell pseudo bulk profiles for the publicly 
#        available sources to perform cell deconvolution using the 
#        GBMDeconvoluteR marker-panel in order to validate the prevalence/associations 
#        that have already been done with the paired, bulk-RNAseq.
# PACKAGES ---------------------------------------------------------------------
library(readr)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(janitor)
library(IMCfuncs)
library(ggplot2)

# I/O --------------------------------------------------------------------------
io <- list(
    inputs = list(
        combined_pb = "Data/public_single_cell/combined_cpm_pseudobulks.rds",
        combined_pb_metadata = "Data/public_single_cell/combined_pseudobulks_metadata.rds",
        cleaned_markers = "Data/public_single_cell/GBMDeconv_cell_markers_and_labels_cleaned_2025-11-19T11-18-40.rds",
        gbmdeconv_markers = "Data/metadata/GBMDeconv_cell_markers_and_labels_2025-11-03T11-01-29.rds",
        ti_markers = "Data/public_single_cell/Wang_et_al_2017_GBM_TI_markers.rds"
    ),
    outputs = list(
        scores = "Data/public_single_cell",
        plots = "Data/public_single_cell/QC"
    ),
    plots = list()
)

source("Process/functions/GBMDeconvoluteR.R")

# LOAD DATA --------------------------------------------------------------------
pb_data <- readRDS(io$inputs$combined_pb)
pb_metadata <- readRDS(io$inputs$combined_pb_metadata)
markers <- readRDS(io$inputs$cleaned_markers)
markers$TI_genes <- readRDS(io$inputs$ti_markers)

# VISUALISE BATCH EFFECTS ------------------------------------------------------
# Just the Nomura data

nomura_meta <- pb_metadata %>% filter(dataset == "nomura")
nomura_data <- pb_data[, nomura_meta$id]
nomura_cpm_log <- log2(nomura_data + 1)

wang_meta <- pb_metadata %>% filter(dataset == "wang")
wang_data <- pb_data[, wang_meta$id]


pca <- prcomp(t(nomura_cpm_log), scale.=FALSE)

pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    tissue_source = nomura_meta$tissue_source,
    lab = nomura_meta$lab,
    patient = nomura_meta$patient,
    surgery = nomura_meta$surgery
)

svglite::svglite(
    filename = nf("nomura_pseudobulk_PCA_by_tissue_source.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = tissue_source)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()


svglite::svglite(
    filename = nf("nomura_pseudobulk_PCA_by_lab_source.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = lab)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()


svglite::svglite(
    filename = nf("nomura_pseudobulk_PCA_by_surgery.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = surgery)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()

pb_cpm_log <- log2(pb_data + 1)

pca <- prcomp(t(pb_cpm_log), scale.=FALSE)

pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    tissue_source = pb_metadata$tissue_source,
    lab = pb_metadata$lab,
    patient = pb_metadata$patient,
    surgery = pb_metadata$surgery
)

svglite::svglite(
    filename = nf("pseudobulk_PCA_by_tissue_source.svg", io$outputs$plots),
    width = 20, height = 20
    )

ggplot(pca_df, aes(x = PC1, y = PC2, colour = tissue_source)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()


svglite::svglite(
    filename = nf("pseudobulk_PCA_by_surgery_source.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = surgery)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()


combat_logcpm <- sva::ComBat(
    dat = pb_cpm_log,
    batch = pb_metadata$tissue_source
)

pca_after <- prcomp(t(combat_logcpm))

pca_df_after <- data.frame(
    PC1 = pca_after$x[, 1],
    PC2 = pca_after$x[, 2],
    tissue_source = pb_metadata$tissue_source,
    lab = pb_metadata$lab,
    patient = pb_metadata$patient,
    surgery = pb_metadata$surgery
)

svglite::svglite(
    filename = nf("ComBat_pseudobulk_PCA_by_tissue_source.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df_after, aes(x = PC1, y = PC2, colour = tissue_source)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca_after)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca_after)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()


svglite::svglite(
    filename = nf("pseudobulk_PCA_by_surgery_source.svg", io$outputs$plots),
    width = 20, height = 20
)

ggplot(pca_df_after, aes(x = PC1, y = PC2, colour = surgery)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    viridis::scale_color_viridis(option = "H", discrete = TRUE) +
    labs(
        title = "PCA of Combined Pseudobulk Expression",
        x = paste0("PC1 (", round(summary(pca_after)$importance[2,1] * 100, 1), "% var)"),
        y = paste0("PC2 (", round(summary(pca_after)$importance[2,2] * 100, 1), "% var)")
    ) +
    IMCfuncs::facetted_cell_prop_theme() 

dev.off()

saveRDS(combat_logcpm, nf("ComBat_corrected_pseudobulks.rds", io$outputs$scores))

# RUN GBMDECONVOLUTER ----------------------------------------------------------

out <- list(
    all = NULL,
    nomura = NULL,
    wang = NULL
    ) 

out$all <- run_GBMDeconvoluteR(
    data = pb_data,
    markers = markers$markers,
    TI_markers = markers$TI_genes
)

out$nomura <- run_GBMDeconvoluteR(
    data = nomura_data,
    markers = markers$markers,
    TI_markers = markers$TI_genes
)

out$wang <- run_GBMDeconvoluteR(
    data = wang_data,
    markers = markers$markers,
    TI_markers = markers$TI_genes
)

saveRDS(out, nf("Raw_GBMDeconvoluteR_scores.rds", io$outputs$scores))


cleaned <- map(out, transform_scores)

cleaned_meta <- list(
    all = pb_metadata,
    nomura = nomura_meta,
    wang = wang_meta
)


cleaned <- map2(cleaned, cleaned_meta, ~{
    
    out <- left_join(
        x = .x, y = .y, by = join_by("Mixture" == "id")
        ) %>%
        dplyr::relocate(cell_types, .after = last_col())
    
    return(out)
    
})
saveRDS(cleaned, nf("Processed_GBMDeconvoluteR_scores.rds", io$output$scores))


# END --------------------------------------------------------------------------