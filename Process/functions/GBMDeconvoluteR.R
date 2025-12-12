# BACKEND FUNCTIONS ------------------------------------------------------------
check_user_data <- function(input){
    
    if(ncol(input) < 3){
        
        return("Insufficient samples (< 2) in data")
        
    }else if(nrow(input) < 500){
        
        return("Insufficient genes (n < 500) in data")
        
    }else if(!is.character(rownames(input))){
        
        return("Gene symbols not detected in the first column")
        
    }else if(length(unique(colnames(input))) != ncol(input)){
        
        return("Duplicate sample names detected in the data")
        
    }else if(length(unique(rownames(input))) != nrow(input)){
        
        return("Duplicate genes symbols detected in the data")
        
    }else if(any(apply(input[,-1], 2, function(x) any(is.na(x))))){
        
        return("Missing values detected in the data")
        
    }else if(any(apply(input[,-1], 2, function(x) any(x<0)))){
        
        return("Negative values detected in the data")
        
    }else{
        
        NULL
        
    }
    
}

preprocess_data  <- function(input, filter_zero_genes = TRUE, log_matrix = TRUE){
    
    
    if(is.matrix(input) & !is.null(rownames(input))){
        
        exprs_matrix <- tinyscalop::exprs_levels(
            m = input,
            bulk = TRUE,
            log_scale = log_matrix
        )
    }else{
        
        genes <- input[[1]]
        
        exprs_matrix <- as.matrix(input[,-1])
        
        rownames(exprs_matrix) <- genes
        
        exprs_matrix <- tinyscalop::exprs_levels(
            m = exprs_matrix,
            bulk = TRUE,
            log_scale = TRUE
        ) 
    }
    
    if(filter_zero_genes){
        
        exprs_matrix[Matrix::rowSums(exprs_matrix) != 0, ,drop = FALSE]
        
    } else return(exprs_matrix)
    
}

clean_cell_labels <- function(cell_labels) {
    cell_labels <- str_replace_all(cell_labels, "s$", "")
    cell_labels <- str_replace_all(cell_labels, "_", " ")
    cell_labels <- str_replace_all(cell_labels, "TAM", "Macrophage")
    return(cell_labels)
}


check_marker_coverage <- function(input, 
                                  input_gene_col = NULL, 
                                  markers,
                                  marker_col = "marker",
                                  celllabel_col = "cell_label",
                                  celltype_col = "cell_type",
                                  conserved_min = 50
                                  ) {
    
    if(!is.data.frame(markers)){
        stop("Markers must be provided as a data.frame")
    } else if(!all(c(marker_col, celltype_col, celllabel_col) %in% colnames(markers))){
        stop("Markers data.frame must contain specified columns for marker_col, celltype_col, and celllabel_col")
    }
    
    markers <- markers %>%
        dplyr::rename(
            "marker" = !!marker_col,
            "cell_type" = !!celltype_col,
            "cell_label" = !!celllabel_col
        )
    
    if(is.matrix(input) & !is.null(rownames(input))){
        markers$missing <- !markers$marker %in% rownames(input)
    } else if(!is.null(input_gene_col)){
        markers$missing <-  !markers$marker %in% input[[input_gene_col]]
    } else{
        stop("Input data must be a matrix with rownames or a data.frame with specified input_gene_col")
    }
    
    marker_coverage <- markers %>%
        dplyr::group_by(cell_label, cell_type) %>%
        dplyr::summarise(
            total_markers = dplyr::n(),
            present_percent = length(which(!missing)),
            missing_percent = length(which(missing)), .groups = "keep"
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(min_coverage = missing_percent < conserved_min) %>%
        dplyr::arrange(dplyr::desc(missing_percent))

    if(any(!marker_coverage$min_coverage)){

        error_msg = glue::glue("More than {conserved_min}% of the genes in the uploaded data are
                           missing in the selected marker gene list.
                           <br><br> <i>Please see the Markers tab for a detailed breakdown</i>")
        
        
        return(
            list(
                error_msg = shiny::HTML(error_msg),
                coverage_table = marker_coverage
            )
        )
        
    }else return(NULL)
    
}

gather_markers <- function(marker_df,
                           cell_label_filter = NULL,
                           cell_type_filter = NULL,
                           negate = FALSE,
                           cell_label_col = "cell_label",
                           cell_type_col = "cell_type",
                           marker_col = "marker",
                           aggregate_by = c("cell_type", "cell_label"), 
                           return_df = FALSE) {
    
    ## --- Defensive checks --------------------------------------------------  
    
    if (!is.data.frame(marker_df)) {
        stop("`marker_df` must be a data.frame.")
    }
    
    aggregate_by <- match.arg(aggregate_by)
    
    if (!is.logical(negate) || length(negate) != 1) {
        stop("`negate` must be a single logical value (TRUE or FALSE).")
    }
    
    required_cols <- unique(c(cell_label_col, cell_type_col, marker_col))
    missing_cols <- required_cols[!required_cols %in% colnames(marker_df)]
    
    if (length(missing_cols) > 0) {
        stop(
            "marker_df does not contain the specified column names: ",
            paste(missing_cols, collapse = ", ")
        )
    }
    
    if (!is.null(cell_label_filter) && !is.atomic(cell_label_filter)) {
        stop("`cell_label_filter` must be an atomic vector.")
    }
    if (!is.null(cell_type_filter) && !is.atomic(cell_type_filter)) {
        stop("`cell_type_filter` must be an atomic vector.")
    }
    
    if (negate && is.null(cell_label_filter) && is.null(cell_type_filter)) {
        stop("When `negate = TRUE`, at least one filter must be supplied.")
    }
    
    if (nrow(marker_df) == 0) {
        warning("marker_df has zero rows; returning empty list.")
        return(setNames(list(), character(0)))
    }
    
    ## --- Build mask --------------------------------------------------------
    
    mask <- rep(TRUE, nrow(marker_df))
    
    if (!is.null(cell_label_filter)) {
        mask <- mask & (marker_df[[cell_label_col]] %in% cell_label_filter)
    }
    
    if (!is.null(cell_type_filter)) {
        mask <- mask & (marker_df[[cell_type_col]] %in% cell_type_filter)
    }
    
    if (negate) {
        mask <- !mask
    }
    
    filtered <- marker_df[mask, , drop = FALSE]
    
    if (nrow(filtered) == 0) {
        warning("Filtering resulted in zero rows; returning empty list.")
        return(setNames(list(), character(0)))
    }
    
    if(return_df){
        return(filtered)
    } else{
        group_col <- switch(
            aggregate_by,
            "cell_type" = cell_type_col,
            "cell_label" = cell_label_col
        )
        
        marker_list <- split(filtered[[marker_col]], filtered[[group_col]])
        
        return(marker_list)
    }

}


refine_neoplastic_markers <- function(input, 
                                      markers, 
                                      filter_threshold = 0.4, 
                                      TI_only = T, 
                                      TI_markers,
                                      neoplastic_celltype_label = "Cancer"
                                      ) {
    
    set.seed(1234)
    
    refined_markers <- tryCatch({
        
        tinyscalop::filter_signatures(
            m = input,
            sigs = gather_markers(
                marker_df = markers, 
                cell_type_filter = neoplastic_celltype_label, 
                aggregate_by = "cell_label"
                ),
            filter.threshold = filter_threshold
        )
        
    },
    error = function(e) return(NULL),
    warning = function(w) return(NULL)
    )
    
    if(!is.null(refined_markers)){
        
        refined_markers  = data.frame(
            "marker" = unlist(refined_markers, use.names = FALSE),
            "cell_label" = rep(names(refined_markers), lengths(refined_markers)),
            "cell_type" = neoplastic_celltype_label
        )
        
        if(TI_only){
            
            refined_markers = refined_markers[refined_markers$marker %in% TI_markers,]
            
        }
        
        non_neoplastic_markers <- gather_markers(
            marker_df = markers, 
            cell_type_filter = neoplastic_celltype_label, 
            aggregate_by = "cell_label",
            negate = TRUE, 
            return_df = TRUE
        )
        

        refined_markers <- refined_markers %>%
            tidyr::drop_na() %>%
            rbind(non_neoplastic_markers) %>%
            dplyr::rename(
                "HUGO symbols" = "marker",
                "Cell population" = "cell_label" 
            )
        
        return(refined_markers)
        
    } else stop("Neoplastic marker refinement failed")
  
    
}

score_data <-  function(exprs_data, markers, out_digits = 3){
    
    scores = MCPcounter::MCPcounter.estimate(exprs_data, featuresType = "HUGO_symbols", genes = markers)
    
    as.data.frame(t(scores)) %>%
        
        tibble::rownames_to_column(var = "Mixture") %>%
        
        mutate(across(-Mixture, ~round(.x, digits = out_digits)))
}

# MAIN WRAPPER FUNCTION --------------------------------------------------------
run_GBMDeconvoluteR = function(
        data,
        markers,
        filter_zero = TRUE,
        log_input_matrix = TRUE,
        min_conserved_markers = 50,
        neoplastic_filt_threshold = 0.4,
        TI_only = TRUE,
        TI_markers,
        score_out_digits = 3
) {
    
    if(is.null(check_user_data(data))){
        
        cleaned_data <- preprocess_data(data, filter_zero_genes = filter_zero, log_matrix = log_input_matrix)
        
    } else cat("User data checks failed!")
    
    initial_marker_check <- check_marker_coverage(
        input = cleaned_data,
        markers = markers,
        conserved_min = min_conserved_markers
    )
    
    
    if(!is.null(initial_marker_check)){
        
        cat(initial_marker_check$error_msg)
        
        return(initial_marker_check$coverage_table)
        
    }
    
    refined_markers = refine_neoplastic_markers(
        input = cleaned_data,
        markers = markers,
        filter_threshold = neoplastic_filt_threshold,
        TI_only = TI_only,
        TI_markers = TI_markers
    )
    
    scores = score_data(
        exprs_data = cleaned_data,
        markers = refined_markers,
        out_digits = score_out_digits
    )
    
    return(scores)
    
}