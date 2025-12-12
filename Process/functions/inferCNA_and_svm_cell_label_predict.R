prepare_chunks <- function(
        file_path,
        dataset_name = NULL,
        ref_genes, 
        chunk_size = 20000, 
        output_dir = io$outputs$chunked
        ) {
    if (is.null(dataset_name)){
        dataset_name <- gsub("\\.rds$", "", basename(file_path))
    }
    
    message(paste(">>> Loading:", dataset_name))
    
    obj <- readRDS(file_path)

    common_genes <- intersect(rownames(obj), ref_genes)
    common_genes <- common_genes[!grepl("^MT-", common_genes)]
    
    message(paste("    Genes matching reference:", length(common_genes)))
    
    raw_counts <- obj[common_genes, ]
    
    all_cells <- colnames(raw_counts)
    n_chunks <- ceiling(length(all_cells) / chunk_size)
    chunks <- split(all_cells, sample(1:n_chunks, length(all_cells), replace = TRUE))
    
    created_files <- c()
    
    for (i in 1:length(chunks)) {
        
        chunk_cells <- chunks[[i]]
        chunk_mat <- raw_counts[, chunk_cells, drop = FALSE]
        
        fname <- file.path(output_dir, paste0(dataset_name, "_chunk_", i, ".rds"))
        
        saveRDS(chunk_mat, fname)
        created_files <- c(created_files, fname)
        
       message(paste("    Saved chunk", i, "of", n_chunks))
    }
    
    rm(raw_counts, obj)
    gc()
    
    return(created_files)
}


predict_chunk_liblinear <- function(file_path, model, ref_genes) {
    
    dataset_name <- gsub("\\.rds$", "", basename(file_path))
    message(paste0(">>> Predicting: ", dataset_name))
    
    raw_counts <- readRDS(file_path)
    
    lib_sizes <- Matrix::colSums(raw_counts)
    lib_sizes[lib_sizes == 0] <- 1
    chunk_cpm <- raw_counts %*% Matrix::Diagonal(x = 1e6 / lib_sizes)
    chunk_log <- log2(chunk_cpm + 1)
    
    current_genes <- rownames(chunk_log)
    common_genes <- intersect(current_genes, ref_genes)
    missing_genes <- setdiff(ref_genes, current_genes)
    
    subset_mat <- chunk_log[common_genes, , drop = FALSE]
    

    input_mat <- as.matrix(subset_mat)
    
    if(length(missing_genes) > 0) {
        zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(input_mat))
        rownames(zero_mat) <- missing_genes
        input_mat <- rbind(input_mat, zero_mat)
    }

    input_mat <- input_mat[ref_genes, , drop = FALSE]
    
    # Transpose: LiblineaR expects Rows = Cells
    if(ncol(input_mat) > 0) {
        
        # Predict with proba = TRUE
        # newx must be the transposed matrix
        pred_obj <- predict(model, newx = t(input_mat), proba = TRUE)
        
        # Extract Probabilities
        probs <- pred_obj$probabilities
        
        # Get max prob and label
        max_probs <- apply(probs, 1, max)
        
        # LiblineaR returns predictions in $predictions
        labels <- as.character(pred_obj$predictions)
        
        # Rejection Threshold
        labels[max_probs < 0.7] <- "Unassigned"
        
        results <- data.frame(
            cell_id = colnames(raw_counts),
            svm_label = labels,
            svm_prob = max_probs,
            dataset = dataset_name
        )
    } else {
        results <- data.frame(cell_id=character(), svm_label=character(), svm_prob=numeric())
    }
    
    return(results)
}
