# params <- list(
#     mask_filt = c("F0001"),
#     colour_by = "cell_label",
#     labels_to_highlight = c("Oligodendrocyte", "Endothelial", "MES3"),
#     cell_fills = spe@metadata$colors$cells$fills,
#     cell_outlines = spe@metadata$colors$cells$outlines,
#     thick_outlines = TRUE,
#     cell_id_col = "cell_ID",
#     img_id_col = "fov",
#     exprs_val_name = "data",
#     out_dir = io$output$temp_out,
#     bkg_color = "grey98",
#     unlabelled_color = "grey93",
#     label_colour = "grey10",
#     mask_name_mcol_col = "patient_surgery"
# )


check_mask_dims <-  function(masks) {
    dims = sapply(masks, dim)
    # drop channel dimension if present
    dims = dims[1:2, , drop = FALSE]
    
    # are all equal?
    all_equal = all(apply(dims, 1, function(x) length(unique(x)) == 1))
    
    if (all_equal) {
        # return single integer if square, otherwise both
        w = unique(dims[2, ])
        h = unique(dims[1, ])
        if (w == h) {
            message("✅ All masks are square with dimension: ", w, " × ", w, " pixels.")
            return(w)
        } else {
            message("✅ All masks have the same rectangular size: ", w, " × ", h, " pixels.")
            return(c(width = w, height = h))
        }
    } else {
        warning("⚠️ Masks have differing dimensions. See details below.")
        print(t(dims))
        return(NULL)
    }
}


plot_cell_objects <- function(masks,
                              spe,
                              cell_id_col = "cell_ID",
                              img_id_col = "fov",
                              exprs_val_name = "data",
                              mask_filt,
                              labels_to_highlight,
                              colour_by = "cell_label",
                              cell_fills = spe@metadata$colors$cells$fills,
                              cell_outlines = spe@metadata$colors$cells$outlines,
                              thick_outlines = TRUE,
                              bkg_color = "grey98",
                              unlabelled_color = "grey93",
                              label_colour = "grey15",
                              mask_name_mcol_col = "patient_surgery",
                              image_titles = NULL,
                              scale_bar_info = NULL,
                              out_dir = io$output$temp_out
                              ){

    stopifnot(all(mask_filt %in% names(masks)))
    stopifnot(colour_by %in% colnames(spe@colData))
    stopifnot(exprs_val_name %in% assayNames(spe))
    stopifnot(all(labels_to_highlight %in% spe[[colour_by]]))

    cur_masks <- masks[names(masks) %in% mask_filt]
    cur_spe <- spe[, spe[[colour_by]] %in% labels_to_highlight]

    cell_label_fills <- cell_fills[labels_to_highlight]
    
    cur_fill <- list(cell_label_fills)
    names(cur_fill) <- colour_by
    
    cur_outline <- cell_outlines[labels_to_highlight]
    
    cur_spe$outline_col <- cur_outline[cur_spe[[colour_by]]]
    
    filename_add <- unique(cur_masks@elementMetadata[[mask_name_mcol_col]])
    img_id_add <- paste0(unique(mask_filt), collapse = "_")
    
    if(length(filename_add) >= 1){
        filename_add <-  paste0(filename_add, collapse = "_")
        out_filename <- paste0(img_id_add, "-" ,filename_add, ".png")
    }else{
        out_filename <- paste0(img_id_add, ".png")
    }
    
    if(is.null(image_titles)){
        image_title_text <- unique(cur_masks@elementMetadata[[img_id_col]])
        } else image_title_text <- image_titles 
    
    # plotting the cell masks coloured by specific cells
    plotCells(mask = cur_masks,
              object = cur_spe,
              cell_id = cell_id_col,
              img_id = img_id_col,
              exprs_values = exprs_val_name,
              outline_by = "outline_col",
              colour_by = colour_by,
              background_colour = bkg_color,
              missing_colour = unlabelled_color,
              colour = cur_fill,
              display = "single",
              image_title = list(text = image_title_text, colour = label_colour),
              scale_bar = scale_bar_info,
              thick = thick_outlines,
              save_plot = list(
                  filename = nf(
                      filename = out_filename,
                      filepath = out_dir
                  ),
                  scale = 3
              )
    )

}
