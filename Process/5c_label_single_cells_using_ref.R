# Author: Shoaib Ajaib
# Date: 03/12/2025
# Usage: This script takes the single cell reference and uses it annotate
#        single cells from Wang et al. 2022 and Nomura et al. 2025 with malignant
#        cell states (Neftel et al. 2019) and normal cell types (GBMap).
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
library(LiblineaR)
library(Seurat)

# I/O --------------------------------------------------------------------------
io <- list(
    inputs = list(
        sc_data = list(
            wang = "Data/public_single_cell/wang_et_al_2022/wang_et_al_all_single_cells.rds",
            nomura_1 = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_12.rds",
            nomura_2 = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_3.rds"
        ),
        sc_reference = "Data/public_single_cell/cell_labelling/GBM_cell_signature_reference_2025-12-03T15-04-54.rds",
        sc_svm_model = "Data/public_single_cell/cell_labelling/GBM_cell_label_svm_model_2025-12-03T17-59-26.rds",
        cleaned_markers = "Data/public_single_cell/GBMDeconv_cell_markers_and_labels_cleaned_2025-11-19T11-18-40.rds",
        svm_3k_model = "Data/public_single_cell/cell_labelling/GBM_cell_label_3k_svm_model_2025-12-03T23-36-12.rds",
        svm_3k_genes = "Data/public_single_cell/cell_labelling/svm_3k_model_top_genes_2025-12-03T23-37-18.rds"
    ),

    outputs = list(
        labelled = "Data/public_single_cell/cell_labelling"
    ),
    
    plots = list()
)

io$outputs$chunked <- nd("chunked_data", path = io$outputs$labelled, add_timestamp = FALSE)

# LOAD DATA --------------------------------------------------------------------
sc_ref <- readRDS(io$inputs$sc_reference)
markers <- readRDS(io$inputs$cleaned_markers)$markers

sc_data <- readRDS(io$inputs$sc_data$wang)

markers$in_wang <- markers$marker %in% rownames(sc_data)

markers <- markers %>% 
    filter(in_wang) %>%
    split(., as.factor(.$cell_label)) %>%
    lapply(\(x) x[["marker"]])

keep <- grep("opc normal", names(markers), ignore.case = T, value = T, invert = T)
markers <- markers[keep]

# MARKER-BASED SCORING ---------------------------------------------------------
sobj <- CreateSeuratObject(counts = sc_data)
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", verbose = FALSE)

sobj <- AddModuleScore(sobj, features = markers, name = "Score")

# 5. Extract Scores
# The metadata now contains columns: Score1, Score2, ..., ScoreN
# corresponding to your marker_list order
score_cols <- paste0("Score", 1:length(markers))
scores <- sobj@meta.data[, score_cols]
colnames(scores) <- names(markers)

# 6. Assign Label (Max Score)
# We find the column with the highest value for each row
max_scores <- apply(scores, 1, max)
final_labels <- colnames(scores)[max.col(scores, ties.method = "first")]

# 7. Rejection (Optional)
# If the highest score is very low (e.g. < 0), the cell expresses NONE of the markers well.
final_labels[max_scores < 0.1] <- "Unassigned"

# 8. Return Metadata
results <- data.frame(
    cell_id = colnames(sc_data),
    marker_label = final_labels,
    max_score = max_scores
)

cell_types <- readRDS(io$inputs$cleaned_markers)
cell_types <- c(unlist(cell_types$labels$main,use.names = FALSE), "Unassigned")


p <- results %>%
    mutate(
        across(marker_label, ~factor(., levels = cell_types))
    ) %>%
    ggplot(aes(x = marker_label)) +
    geom_bar(fill = "slateblue") +
    theme_minimal() +
    xlab("Predicted Cell Label") +
    ylab("Cell Count") +
    ggtitle("Just scoring markers - Wang et al (158k cells)") +
    IMCfuncs::facetted_cell_prop_theme()

svglite::svglite(
    nf("wang_marker_scoring_with_opc.svg", io$outputs$labelled),
    width = 15, height = 10
)
print(p)
dev.off()



# TRAIN THE SVM MODEL USING ALL GENES ------------------------------------------
ref_mat <- t(as.matrix(sc_ref$matrix))
ref_labels <- as.factor(sc_ref$labels)

# kernel = "linear" and probability = TRUE are needed for the Rejection method
svm_model <- svm(
    x = ref_mat,
    y = ref_labels,
    kernel = "linear",
    probability = TRUE,
    cost = 1,
    scale = TRUE # Scales features to mean=0 sd=1 (Recommended)
)


saveRDS(
    object = svm_model,
    file = nf("GBM_cell_label_svm_model.rds", io$outputs$labelled)
)

# TRAIN THE SVM MODEL USING 3000 VARIABLE GENES --------------------------------
gene_vars <- apply(sc_ref$matrix, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:3000]

ref_subset <- t(as.matrix(sc_ref$matrix[top_genes, ])) # Dense matrix is now tiny (~200MB)

train_mat <- t(as.matrix(sc_ref$matrix[top_genes, ]))
train_labels <- as.factor(sc_ref$labels)

svm_model_fast <- LiblineaR(
    data = train_mat, 
    target = train_labels, 
    type = 7, 
    bias = TRUE, 
    verbose = FALSE
)

saveRDS(
    object = svm_model_fast,
    file = nf("GBM_cell_label_3k_svm_model.rds", io$outputs$labelled)
)


saveRDS(top_genes, file = nf("svm_3k_model_top_genes.rds", io$outputs$labelled))

# CHUNKING THE SINGLE CELL DATA FOR CELL LABEL PREDICTION ----------------------
source("Process/functions/inferCNA_and_svm_cell_label_predict.R")
ref_genes <- rownames(sc_ref$matrix)

# chunks <- prepare_chunks(
#     file_path = io$inputs$sc_data$nomura_2,
#     dataset_name = "nomura_2",
#     ref_genes = ref_genes, 
#     output_dir = io$outputs$chunked
#     )

# PREDICTING CELL LABELS USING REF ---------------------------------------------
model <- readRDS(io$inputs$svm_3k_model)
genes <- readRDS(io$inputs$svm_3k_genes)
all_chunk_files <- list.files(io$outputs$chunked, full.names = TRUE)

all_results <- list()

for (f in all_chunk_files) {
    res <- predict_chunk_liblinear(f, model, genes)
    all_results[[basename(f)]] <- res
    gc()
}

labelled_cells <- bind_rows(all_results)

write.csv(
    labelled_cells, 
    nf("Wang_Nomura_labelled_cells.csv", io$outputs$labelled), 
    row.names = FALSE)

# SIGNTURE REFERENCE QC  -------------------------------------------------------
model <- readRDS(io$inputs$svm_3k_model)
genes <- readRDS(io$inputs$svm_3k_genes)

markers <- markers %>% 
    split(., as.factor(.$cell_label)) %>%
    lapply(\(x) x[["marker"]])

labelled_cells <- readr::read_csv("Data/public_single_cell/cell_labelling/Wang_Nomura_labelled_cells_2025-12-03T23-59-55.csv")

# plot the labelled cell counts

markers$labels$main$Normal <- c("OPC normal", "Radial Glial",  "Astrocyte", "Oligodendrocyte", "Neuron") 

cell_types <- c(unlist(markers$labels$main,use.names = FALSE), "Unassigned")

p <- labelled_cells %>%
    mutate(
        across(svm_label, ~factor(., levels = cell_types))
        ) %>%
    ggplot(aes(x = svm_label)) +
    geom_bar(fill = "slateblue") +
    theme_minimal() +
    xlab("Predicted Cell Label") +
    ylab("Cell Count") +
    ggtitle("Distribution of Predicted Cell Labels in Single Cell Data") +
    IMCfuncs::facetted_cell_prop_theme()

svglite::svglite(
    nf("Predicted_cell_label_distribution_wang_nomura.svg", io$outputs$labelled),
    width = 15, height = 10
)
print(p)
dev.off()



# END --------------------------------------------------------------------------