# Author: Shoaib Ajaib
# Date: 01/12/2025
# Usage: This script creates a IDHwt GBM specific single-cell reference using 
#        data from two studies: 
#        - Neftel et al., 2019 (GSE131928)
#        - GBMap extended reference data 
#       The reference includes both malignant cell states, immune cells, normal
#       brain cells and also vasculature. This can then be used to annotate other
#       single-cell datasets.
# 
# PACKAGES ---------------------------------------------------------------------
library(readr)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(janitor)
library(IMCfuncs)
library(ggplot2)
library(Matrix)

# I/O --------------------------------------------------------------------------
io <- list(
    inputs = list(
        wang_sc_data = "Data/public_single_cell/wang_et_al_2022/wang_et_al_all_single_cells.rds",
        nomura_sc_data_1 = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_12.rds",
        nomura_sc_data_2 = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_3.rds",
        cleaned_markers = "Data/public_single_cell/GBMDeconv_cell_markers_and_labels_cleaned_2025-11-19T11-18-40.rds"
    ),
    neftel = list(
        expression = "Data/public_single_cell/cell_labelling/Neftel_2019/expression/IDHwtGBM.processed.SS2.logTPM.txt",
        metadata = "Data/public_single_cell/cell_labelling/Neftel_2019/metadata/IDHwt.GBM.Metadata.SS2.txt"
    ),
    gbmap = list(
        expression = "Data/public_single_cell/cell_labelling/GBMap/data/cell_reference/gbmap_ref_matrix.mtx",
        genes = "Data/public_single_cell/cell_labelling/GBMap/data/cell_reference/gbmap_ref_genes.csv",
        metadata = "Data/public_single_cell/cell_labelling/GBMap/data/cell_reference/gbmap_ref_metadata.csv"
        ),
    
    outputs = list(
        labelled = "Data/public_single_cell/cell_labelling"
        ),
    
    plots = list()
)

# CREATE MALIGNANT CELL REFERENCE ----------------------------------------------
reference <- readr::read_table(io$neftel$expression, show_col_types = FALSE)

metadata <- list()

metadata$ss_cells <- readr::read_table(
    file = io$neftel$metadata,
    col_names = TRUE,
    show_col_types = FALSE
        )

metadata$malignant_cell_scores <- metadata$ss_cells %>%
    dplyr::filter(GBMType == "Adult") %>%
    dplyr::filter(CellAssignment == "Malignant") %>%
    select(cell_id = NAME,
           MES1= MESlike1, 
           MES2= MESlike2,
           AC = AClike,
           OPC = OPClike,
           NPC1 = NPClike1,
           NPC2 = NPClike2
    ) %>%
    tidyr::drop_na() %>%
    tibble::column_to_rownames(var = "cell_id")
    

label_malignant_cells <- function(score_df, four_states = FALSE) {
    
    if(four_states){
        state_scores <- tinyscalop::as_four_state_gbm(score_df)
    } else state_scores <- score_df
    
    # 1. Identify Primary State
    # -------------------------------------
    max_idx <- max.col(state_scores, ties.method = "first")
    states  <- colnames(state_scores)
    
    primary_labels <- states[max_idx]
    primary_scores <- state_scores[cbind(1:nrow(state_scores), match(primary_labels, colnames(state_scores)))]
    
    # 2. Calculate Population Thresholds (Criterion 2)
    # -------------------------------------
    # "Higher than that of 10% of the cells that map to this meta-module"
    # We calculate the 10th percentile score for each state based on cells 
    # that assigned that state as their Primary.
    
    thresholds <- sapply(states, function(s) {
        # Extract scores of cells where 's' is the primary state
        p_scores <- primary_scores[primary_labels == s]
        
        # If no cells map to this state primarily, set threshold to Infinity (impossible to be hybrid)
        if(length(p_scores) == 0) return(Inf)
        
        return(quantile(p_scores, probs = 0.1))
    })
    names(thresholds) <- states
    
    message("State thresholds (10th percentile of primary):")
    print(thresholds)
    
    out_df <- state_scores
    out_df$primary_label <- primary_labels
    out_df$primary_score <- primary_scores
    
    
    return(
        list(
            thresholds = thresholds,
            scores = out_df
    )
    )
}

find_hybrids <- function(
        scored_df, 
        thresholds, 
        min_score = 1.0, 
        min_diff = 0.3
        ) {
    
    # 1. Setup: Isolate the numeric score matrix
    # ------------------------------------------
    score_cols <- names(thresholds)
    
    # Check if we are in "6-State Mode" (Subtypes present)
    has_subtypes <- all(c("MES1", "MES2", "NPC1", "NPC2") %in% score_cols)
    
    # Create a matrix for calculation
    mat <- as.matrix(scored_df[, score_cols])
    
    # 2. Identify Ranks (Primary, Secondary, Tertiary)
    # ------------------------------------------
    
    # A. Primary (Rank 1)
    idx_1 <- max.col(mat, ties.method = "first")
    primary_labels <- colnames(mat)[idx_1]
    
    # B. Secondary (Rank 2)
    # Mask Primary to find Secondary
    mat_temp <- mat
    mat_temp[cbind(1:nrow(mat), idx_1)] <- -Inf
    idx_2 <- max.col(mat_temp, ties.method = "first")
    sec_scores <- mat[cbind(1:nrow(mat), idx_2)]
    sec_labels <- colnames(mat)[idx_2]
    
    # C. Tertiary (Rank 3)
    # Mask Secondary to find Tertiary
    mat_temp[cbind(1:nrow(mat), idx_2)] <- -Inf
    idx_3 <- max.col(mat_temp, ties.method = "first")
    tert_scores <- mat[cbind(1:nrow(mat), idx_3)]
    
    # 3. Apply Exclusion Logic (MES1/MES2 & NPC1/NPC2)
    # ------------------------------------------
    # Default assumption: They are different lineages
    is_different_lineage <- rep(TRUE, nrow(mat))
    
    if(has_subtypes) {
        # Define Lineage Map
        # We strip the numbers to get the base lineage (MES1 -> MES, NPC2 -> NPC)
        # AC and OPC remain unchanged
        get_lineage <- function(x) gsub("[12]", "", x)
        
        prim_lineage <- get_lineage(primary_labels)
        sec_lineage  <- get_lineage(sec_labels)
        
        # Check: Is the lineage the same? (e.g. MES1 vs MES2)
        # If Primary is MES1 and Secondary is MES2, is_different_lineage = FALSE
        is_different_lineage <- prim_lineage != sec_lineage
    }
    
    # 4. Apply Hybrid Criteria (Vectorized)
    # ------------------------------------------
    res <- scored_df
    
    # Add Metadata columns
    res$primary_label   <- primary_labels
    res$secondary_label <- sec_labels
    res$secondary_score <- sec_scores
    res$tertiary_score  <- tert_scores
    
    # Crit 0: Must be different lineages (Only relevant for 6-state input)
    res$crit0_diff_lineage <- is_different_lineage
    
    # Crit 1: Secondary Score > min_score (e.g. 1.0)
    res$crit1_high_score <- res$secondary_score > min_score
    
    # Crit 2: Secondary Score > 10th percentile of that state's population
    # We look up the specific threshold for the *secondary* label
    res$crit2_population <- res$secondary_score > thresholds[res$secondary_label]
    
    # Crit 3: Separation > 0.3
    res$crit3_separation <- (res$secondary_score - res$tertiary_score) > min_diff
    
    # 5. Final Hybrid Determination
    # ------------------------------------------
    # A cell is hybrid ONLY if ALL criteria are TRUE
    res$is_hybrid <- res$crit0_diff_lineage & 
        res$crit1_high_score & 
        res$crit2_population & 
        res$crit3_separation
    
    # Create Final Label Column
    # If hybrid: "Primary-Secondary", Else: "Primary"
    # Note: For 6-state input, this keeps the subtype numbers (e.g. "MES1-AC")
    res$final_state <- ifelse(res$is_hybrid,
                              paste(res$primary_label, res$secondary_label, sep = "-"),
                              res$primary_label)
    
    return(res)
}

metadata$labelled_states <- label_malignant_cells(
    score_df = metadata$malignant_cell_scores, 
    four_states = FALSE
    )

metadata$labelled_states <-  find_hybrids(
    scored_df = metadata$labelled_states$scores,
    thresholds =  metadata$labelled_states$thresholds
    )

metadata$labelled_cells <-  metadata$labelled_states %>%
    dplyr::filter(!is_hybrid) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    select(cell_id,
           MES1:NPC2,
           state = final_state
           )
table(metadata$labelled_cells$state)    

malignant_ref <- reference %>%
    select(
        gene = GENE, 
        all_of(metadata$labelled_cells$cell_id)
        ) 

metadata$malignant_labels_vec <- metadata$labelled_cells$state
names(metadata$malignant_labels_vec) <- metadata$labelled_cells$cell_id

metadata$malignant_labels_vec

to_save <- list(
    expression = malignant_ref,
    labels = metadata$malignant_labels_vec
)

# Sense check
all(names(to_save$labels) %in% colnames(to_save$expression))

saveRDS(
    to_save, 
    nf("Neftel2019_SC_state_ref.rds", io$outputs$labelled)
    )

rm(
    malignant_ref, markers, reference, metadata, to_save,
    find_hybrids, label_malignant_cells
    )


# NORMAL CELL REFERNCE FROM GBMAP ----------------------------------------------
reference <- list()

reference$gbmap <- list(data = readMM(io$gbmap$expression))
reference$gbmap$genes <- read.csv(io$gbmap$genes, header=F)$V1
reference$gbmap$metadata <- read.csv(io$gbmap$metadata, row.names=1)

rownames(reference$gbmap$data) <- reference$gbmap$genes
colnames(reference$gbmap$data) <- rownames(reference$gbmap$metadata)

reference$gbmap$labels <- reference$gbmap$metadata$svm_label
names(reference$gbmap$labels) <- rownames(reference$gbmap$metadata)

# Convert the raw GBMap data to CPM & put on a log2 scale
total_counts <- Matrix::colSums(reference$gbmap$data)
reference$gbmap$cpm_data <- Matrix::t(Matrix::t(reference$gbmap$data) / total_counts) * 1e6
reference$gbmap$log2_cpm_data  <- log2(reference$gbmap$cpm_data + 1)

# CLEAN UP THE MALIGNANT CELLS REFERENCE ---------------------------------------
reference$neftel <- readRDS(
    file = list.files(io$outputs$labelled, pattern = "SC_state_ref", full.names = TRUE)
    )

reference$neftel$tpm_data <- reference$neftel$expression %>%
    tibble::column_to_rownames(var = "gene") %>%
    as.matrix()

max_val <- max(reference$neftel$tpm_data)

if (max_val > 50) {
    message("Data appears to be Linear TPM. Log-transforming...")
    reference$neftel$log2_tpm_data <- log2(reference$neftel$tpm_data + 1)
} else {
    message("Data appears to be already Log-transformed. Using as is.")
    reference$neftel$log2_tpm_data <- reference$neftel$tpm_data
}

rm(max_val, total_counts)

# HARMONISE AND DOWNSAMPLE REFERENCE DATA --------------------------------------
target_n <- min(c(table(reference$gbmap$labels), table(reference$neftel$labels)))

set.seed(123)

common_genes <- intersect(
    x = rownames(reference$gbmap$log2_cpm_data),
    y =  rownames(reference$neftel$log2_tpm_data)
    )

print(paste("Common genes found:", length(common_genes)))

# Subset matrices to common genes
reference$gbmap$ref_mat  <- reference$gbmap$log2_cpm_data[common_genes, ]
reference$neftel$ref_mat <- reference$neftel$log2_tpm_data[common_genes, ]

downsample_to_exact_n <- function(mat, labels, n) {
    
    kept_cells <- c()
    unique_types <- unique(labels)
    
    for (type in unique_types) {
        cells_of_type <- names(labels)[labels == type]
 
        if (length(cells_of_type) < n) {
            warning(paste("Type", type, "has fewer than", n, "cells. Keeping all", length(cells_of_type)))
            kept_cells <- c(kept_cells, cells_of_type)
        } else {
            kept <- sample(cells_of_type, n)
            kept_cells <- c(kept_cells, kept)
        }
    }
    
    # Return list
    return(list(
        matrix = mat[, kept_cells],
        labels = labels[kept_cells]
    ))
}

reference$gbmap$ref_mat <- downsample_to_exact_n(
    mat = reference$gbmap$ref_mat, 
    labels = reference$gbmap$labels, 
    n = target_n
    )

reference$neftel$ref_mat <- downsample_to_exact_n(
    mat = reference$neftel$ref_mat, 
    labels = reference$neftel$labels, 
    n = target_n
)

# FINAL REFERENCE DATASET MERGE ------------------------------------------------
# sense check
all(
    c(table(reference$gbmap$ref_mat$labels), 
      table(reference$neftel$ref_mat$labels)
      ) == target_n
)

GBM_sc_ref_matrix <- cbind(
    reference$neftel$ref_mat$matrix, reference$gbmap$ref_mat$matrix
    )

GBM_sc_ref_labels <- c(
    reference$neftel$ref_mat$labels, reference$gbmap$ref_mat$labels
    )

# 3. Save Final Reference
GBM_cell_sig_reference <- list(
    matrix = GBM_sc_ref_matrix,
    labels = GBM_sc_ref_labels
)

saveRDS(
    object = GBM_cell_sig_reference, 
    file = nf(filepath = io$outputs$labelled,
              filename =  "GBM_cell_signature_reference.rds"
              )
        )


# END --------------------------------------------------------------------------