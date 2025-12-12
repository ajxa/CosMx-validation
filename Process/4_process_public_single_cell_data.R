# Author: Shoaib Ajaib
# Date: 10/11/2025
# Usage: This script is intended to clean up the single-cell data obtained from
#        publicly available sources and prepare it for downstream analysis. 
#        Specifically, these data will be used to complete cell deconvolution
#        using the GBMDeconvoluteR marker-panel in order to validate the 
#        prevalence/associations that have already been done with the paired, 
#        bulk-RNAseq that I have already analysed.
# Public Datasets: 
#       1. Wang et al. Nat Cancer 3, 1534–1552 (2022). https://doi.org/10.1038/s43018-022-00475-x
#       2. Nomura et al. Nat Genet 57, 1155–1167 (2025). https://doi.org/10.1038/s41588-025-02167-5
# PACKAGES ---------------------------------------------------------------------
library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(purrr)
library(IMCfuncs)

# I/O --------------------------------------------------------------------------
io <- list(
    inputs = list(
        wang_etal = list(
            data_dir = "Data/public_single_cell/wang_et_al_2022/GSE174554_RAW",
            metadata_file = "Data/public_single_cell/wang_et_al_2022/metadata.xlsx",
            cleaned_metadata = "Data/public_single_cell/wang_et_al_2022/wang_et_al_metadata_filtered.xlsx",
            cleaned_sc = "Data/public_single_cell/wang_et_al_2022/wang_et_al_all_single_cells.rds",
            pseudo_bulks = "Data/public_single_cell/wang_et_al_2022/wang_et_al_pseudobulk_cpm.rds"
        ),
        nomura_etal = list(
            data_dir = "Data/public_single_cell/nomura_et_al_2025/GSE274546_RAW",
            unzipped_raw = "Data/public_single_cell/nomura_et_al_2025/RAW",
            metadata_file = "Data/public_single_cell/nomura_et_al_2025/metadata.xlsx",
            cleaned_metadata = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_metadata_filtered.xlsx",
            cleaned_sc = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_all_single_cells.rds",
            pseudo_bulks = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_pseudobulk_cpm.rds"
        )
    ),
    outputs = list(
        data = list(
            wang_et_al = "Data/public_single_cell/wang_et_al_2022",
            nomura_et_al = "Data/public_single_cell/nomura_et_al_2025"
        )
    ),
    plots = list()
)

# for (i in seq_along(io$inputs$functions)) source(io$inputs$functions[[i]], verbose = FALSE)

# create out directory
ndirs(io$outputs$data)

# create time-stamped output directory
# io$outputs$temp_out <- nd(path = io$outputs$out_dir)

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

rm(find_file, i)

# NOUMRA ET AL. CLEAN METADATA -------------------------------------------------

if(!file.exists(io$inputs$nomura_etal$cleaned_metadata)){
    
    files <-  list.files(io$inputs$nomura_etal$data_dir)
    
    metadata <-  openxlsx::read.xlsx(io$inputs$nomura_etal$metadata_file)
    
    metadata <- metadata %>%
        select(id = Sample.ID,
               tumour_id = Tumor.ID,
               surgery = Primary.or.recurrent,
               age = Age.at.initial.diagnosis,
               sex = Gender,
               tumour_loc = Tumor.location,
               idh_status = IDH.mutation.status,
               PFS = `Surgical.interval.(month)`
        ) %>%
        mutate(across(surgery, ~str_replace(.x, "\\s+$", ""))) %>%
        mutate(patient_id = str_extract(id, "^P\\d+"), .after = id) %>%
        filter(surgery %in% c("Primary", "1st Recurrent")) %>%
        mutate(
            across(surgery, ~str_replace(.x, "1st\\s+", "")),
            across(idh_status, ~str_replace(.x, "wt$", " wildtype")),
            across(age, as.numeric),
            across(PFS, as.numeric)
        ) 
    
    singletons <- names(which(table(metadata$patient_id) == 1))
    if(length(singletons) > 0){
        metadata <- metadata %>% filter(!patient_id %in% singletons)
    }
    
    meta_split <- split(metadata, metadata$patient_id) |>
        lapply(\(x){
            
            to_update <- list(
                age = unique(na.omit(x$age)),
                sex = unique(na.omit(x$sex)),
                PFS = unique(na.omit(x$PFS))
            )
            
            if(all(lengths(to_update) == 1)){
                for(i in names(to_update)){
                    x[[i]] <- to_update[[i]]
                }
                return(x)
            } else{
                stop("Multiple values found for the same patient!")
            }
            
        })
    
    metadata_cleaned <- bind_rows(meta_split)
    
    rm(singletons, meta_split)
    
    file_ids <- unique(str_split(files, "_",simplify = T)[,2])
    
    x <- metadata_cleaned$id %in% file_ids
    
    if(all(x)){
        cli::cli_alert_success("All metadata IDs found in files")
    } else if(any(!x)) cli::cli_alert_warning("{length(which(!x))} metadata IDs not found in files")
    
    openxlsx::write.xlsx(
        x = metadata_cleaned,
        file = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_metadata_filtered.xlsx"
    )
    
    
    
}

if(file.exists(io$inputs$nomura_etal$cleaned_metadata)){
    
    metadata <- openxlsx::read.xlsx(io$inputs$nomura_etal$cleaned_metadata)
}

# NOMURA ET AL. DATA -----------------------------------------------------------
file_paths <-  list.files(io$inputs$nomura_etal$data_dir, full.names = TRUE)
file_names <-  list.files(io$inputs$nomura_etal$data_dir)
file_ids <- str_extract(file_names, "P\\d+T\\d")

file_paths <- file_paths[file_ids %in% metadata$id]
file_ids <- file_ids[file_ids %in% metadata$id]

for(i in seq_along(file_paths)){
    
    if(!file.exists(new_paths[i])){
        R.utils::gunzip(
            filename = file_paths[i], 
            destname = new_paths[i], 
            remove = FALSE
        )
    }
    
}

rm(file_names, i)

read_sparse_rds <- function(path) {
    
    mat_sparse <- Matrix::Matrix(data = readRDS(path), sparse = TRUE)
    
    # Defensive checks
    if (!inherits(mat_sparse, "dgCMatrix")) {
        stop("Object from ", path, " could not be coerced to a sparse dgCMatrix.")
    }
    if (is.null(rownames(mat_sparse))) {
        stop("Sparse matrix from ", path, " has no rownames (genes).")
    }
    
    return(mat_sparse)
}

tictoc::tic("Reading in sparse matrices")
mat_list <- lapply(new_paths, read_sparse_rds)
tictoc::toc()    

# sense checks
lapply(mat_list, \(x){
    if(length(rownames(x)) != length(unique(rownames(x)))) return(FALSE) else return(TRUE)
}) |> unlist() |> all()

unique_genes <- lapply(mat_list, \(x) return(unique(rownames(x))))

bar <- purrr::reduce(unique_genes,\(acc, nxt) if(all(acc == nxt)) return(unique(c(acc, nxt))) else stop("Gene names differ between matrices!"))
all(bar == unique_genes[[1]])    
rm(bar, unique_genes)

# Change cell names to include file ID prefix
mat_list <-  purrr::map2(mat_list, file_ids, ~{
    colnames(.x) <- paste0(.y, "_", seq_len(ncol(.x)))
    return(.x)
})

sapply(mat_list, ncol) |> sum()

tictoc::tic("Building first-half of the UMI matrix")
combined_matrix <-  purrr::reduce(mat_list[1:52], Matrix::cbind2)
tictoc::toc()

mat_list <- mat_list[-(1:52)]

obj_size_gb <- function(x) {
    sprintf("%.3f GB", as.numeric(object.size(x)) / (1024^3))
}

dim(combined_matrix)
obj_size_gb(combined_matrix)

saveRDS(
    combined_matrix, 
    file = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_1.rds"
    )

rm(combined_matrix)
gc()

tictoc::tic("Building second-half of the UMI matrix")
combined_matrix <-  purrr::reduce(mat_list[1:26], Matrix::cbind2)
tictoc::toc()

mat_list <- mat_list[-(1:26)]

dim(combined_matrix)
obj_size_gb(combined_matrix)

saveRDS(
    combined_matrix, 
    file = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_2.rds"
)

rm(combined_matrix)
gc()

tictoc::tic("Building third-half of the UMI matrix")
combined_matrix <-  purrr::reduce(mat_list, Matrix::cbind2)
tictoc::toc()

dim(combined_matrix)
obj_size_gb(combined_matrix)

saveRDS(
    combined_matrix, 
    file = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_3.rds"
)

rm(combined_matrix, mat_list, file_ids, file_paths, new_paths)
gc()

# Combine the firsr two parts into a single matrix
part1 <- readRDS("Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_1.rds")
part2 <- readRDS("Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_2.rds")

combined_matrix <- Matrix::cbind2(part1, part2)

rm(part1, part2)
gc()

nnz = length(combined_matrix@x) + length(part3@x)    # non-zeros in both
ncols = ncol(combined_matrix) + ncol(part3)
bytes_est = nnz * (8 + 4) + (ncols + 1) * 4   # x + i + p
gb_est = bytes_est / 1024^3
sprintf("Final matrix ≈ %.2f GB (ignoring overhead)", gb_est)

saveRDS(
    combined_matrix, 
    file = "Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_12.rds"
)

rm(combined_matrix, part3)
gc()

# FINAL NOMURA ET AL. COMBINED MATRIX ------------------------------------------
part12 <- readRDS("Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_12.rds")
part3 <- readRDS("Data/public_single_cell/nomura_et_al_2025/nomura_et_al_umi_matrix_3.rds")


# WANG ET AL. DATA -------------------------------------------------------------
# snRNA-seq on a cohort of patient-matched primary-recurrent paired specimens (n=86)
# 36 patients profiled paired (primary untreated tumor and matched first recurrence)

# Get rid of the unnecessary files that are not needed for this analysis:

files <-  list.files(io$inputs$wang_etal$data_dir)
remove_files = c("snATAC", "3days", "Control", "Treated", "Sham", "IR",
                 "Transcriptomics.tar.gz", "Proteomics.zip")

rem =  stringr::str_split(files, "_",simplify = T)[,3] %in% remove_files
bad = file.path(data_dir, files[rem])

exists = file.exists(bad)
if (any(!exists)) message("These paths do not exist:\n", paste(bad[!exists], collapse = "\n"))

okay_to_delete = FALSE

if (okay_to_delete) {
    ok = unlink(bad[exists], recursive = FALSE, force = TRUE)
    if (any(ok != 0)) message("Failed to delete:\n", paste(bad[exists][ok != 0], collapse = "\n"))
}

# WANG ET AL. CLEAN METADATA ---------------------------------------------------
files <-  list.files(io$inputs$wang_etal$data_dir)
metadata <-  openxlsx::read.xlsx(io$inputs$wang_etal$metadata_file)

metadata <- metadata %>%
    dplyr::filter(!is.na(`snRNA-seq`)) %>%
    select(id = ID,
           patient = `Pair#`, 
           surgery = Stage,
           age = Age,
           sex = Sex,
           tumour_loc = Tumor.site,
           idh_status = IDH,
           PFS = Elapsed.time.to.recurrence,
           OS = Overall.survival
           ) %>%
    mutate(across(patient, as.numeric)) %>%
    filter(!is.na(patient)) %>%
    filter(idh_status == "IDH wildtype") %>%
    arrange(patient)


file_ids = unique(str_split(files, "_",simplify = T)[,2])
x = metadata$id %in% file_ids

found <- metadata[x,]
not_found <- metadata[!x,]

found_pairs <- found %>%
    group_by(patient) %>%
    summarise(n = n()) %>%
    filter(n == 2) %>%
    pull(patient)

found <- found %>% filter(patient %in% found_pairs)
    
openxlsx::write.xlsx(
    x = found,
    file = "Data/public_single_cell/wang_et_al_2022/wang_et_al_metadata_filtered.xlsx"
    )

metadata <- found
ids <-  unique(metadata$id)


rm(list = grep("^metadata$|^io$|^files$", ls(), invert = TRUE, value = TRUE))

# READ IN SINGLE-CELL DATA -----------------------------------------------------
read_one_triplet <- function(matrix_file) {
    # infer the other two files from the matrix file
    prefix = sub("_matrix\\.mtx\\.gz$", "", matrix_file)
    f_features = paste0(prefix, "_features.tsv.gz")
    f_barcodes = paste0(prefix, "_barcodes.tsv.gz")
    
    # read
    feats = read.table(f_features, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    genes = if (ncol(feats) >= 2) feats[[2]] else feats[[1]]
    bc = read.table(f_barcodes, sep = "\t", header = FALSE, stringsAsFactors = FALSE)[[1]]
    m = Matrix::readMM(matrix_file)
    
    # name rows/cols
    rownames(m) = genes
    colnames(m) = bc
    m
}

# Expand a matrix to include a common set of genes (rows), filling missing with 0
expand_to_genes <- function(m, all_genes) {
    if (identical(rownames(m), all_genes)) return(m)
    out = Matrix::Matrix(0, nrow = length(all_genes), ncol = ncol(m),
                 sparse = TRUE, dimnames = list(all_genes, colnames(m)))
    # place rows by match
    out[match(rownames(m), all_genes), ] = m
    out
}

all_mtx <-  list.files(io$inputs$wang_etal$data_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)

# Extract the sample id sitting between the first and second underscore
extract_id <-  function(path) {
    b = basename(path)
    sub("^GSM[0-9]+_([^_]+).*", "\\1", b)
}

mtx_id <- vapply(all_mtx, extract_id, character(1))
by_id <- split(all_mtx, mtx_id)

# only IDs present in your metadata:
by_id <- by_id[names(by_id) %in% ids]

# read, combine batches within id, and rename cells
mats_by_id <- lapply(names(by_id), function(id) {
    # read every batch (each matrix file is one batch)
    mats = lapply(by_id[[id]], read_one_triplet)
    # combine columns for this id
    m_id = do.call(cbind, mats)
    # rename cells as "<id>_<cellindex>"
    colnames(m_id) = sprintf("%s_%d", id, seq_len(ncol(m_id)))
    m_id
})
names(mats_by_id) <- names(by_id)

# merge all ids to a single genes x cells matrix
all_genes <- sort(unique(unlist(lapply(mats_by_id, rownames))))
mats_by_id_expanded <- lapply(mats_by_id, expand_to_genes, all_genes = all_genes)
counts_all <- do.call(cbind, mats_by_id_expanded)

# Result:
# counts_all is a sparse dgCMatrix with rows = genes, cols = cells
# column names look like: SF10099_1, SF10099_2, …, SF11873_1, …
counts_all

nnz = length(counts_all@x)
cat("Non-zero entries:", nnz, "\n")
cat("Fraction non-zero:", nnz / prod(dim(counts_all)), "\n")

saveRDS(counts_all, file = file.path(io$outputs$wang_et_al, "wang_et_al_all_single_cells.rds"))

rm(list = grep("^metadata$|^counts_all$", ls(), invert = TRUE, value = TRUE))

# PSEUDOBULK THE SINGLE-CELL DATA ----------------------------------------------

# 1) Drop genes with zero counts across all cells
keep_gene <- Matrix::rowSums(counts_all) > 0

cli::cli_alert_info("genes to remove:\t{length(which(!keep_gene))}")
cli::cli_alert_info("genes to keep:\t{length(which(keep_gene))}")

counts_all <-  counts_all[keep_gene, , drop = FALSE]

dim(counts_all)

# 2) Map each cell to its sample id (prefix before the trailing _<cellindex>)
cell_to_id <- sub("_[0-9]+$", "", colnames(counts_all))

# 3) Build a cells → id sparse design matrix with columns in the metadata order
id_levels <- unique(metadata$id)                     
grp <- factor(cell_to_id, levels = id_levels)
G <- Matrix::sparse.model.matrix(~ 0 + grp) 
colnames(G) <- id_levels

# 4) Pseudobulk: sum UMIs per gene within each id
pb_sparse <- counts_all %*% G

# 5) drop genes that are zero across all pseudobulks
keep_gene2 <- Matrix::rowSums(pb_sparse) > 0
pb_sparse <- pb_sparse[keep_gene2, , drop = FALSE]

pb_counts <- as.matrix(pb_sparse)
pb_counts[1:5,1:5]

# 6) Library-size normalisation (CPM; linear scale, no log)
lib <-  colSums(pb_counts)
stopifnot(all(lib > 0))
pb_cpm <- sweep(pb_counts, 2, lib, "/") * 1e6

# 7) Ensure columns are ordered by meta$id
pb_cpm <- pb_cpm[, metadata$id, drop = FALSE]

all(colnames(pb_cpm) == metadata$id)
pb_cpm[1:5, 1:5]

saveRDS(
    pb_cpm, 
    "Data/public_single_cell/wang_et_al_2022/wang_et_al_pseudobulk_cpm.rds"
    )

# END --------------------------------------------------------------------------