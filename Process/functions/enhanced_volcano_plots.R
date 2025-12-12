set_axis_lims <- function(df,
                          y_var = "log2FC",
                          x_var = "p_adj",
                          x_allowance = 0.1,
                          y_allowance = 0) {
    if (!y_var %in% colnames(df)) cli::cli_abort("y_var not found in df supplied")
    if (!x_var %in% colnames(df)) cli::cli_abort("x_var not found in df supplied")
    
    out <- list(
        xlims = c(
            min(df[[y_var]]) * (1 + x_allowance),
            max(df[[y_var]]) * (1 + x_allowance)
        ),
        ylims = c(
            0,
            ceiling(ceiling(-log10(min(df[[x_var]]))) * (1 + y_allowance))
        )
    )
    
    # Sensible defaults when positive or negative log2FC doesn't exist
    if (out$xlims[[1]] >= 0) out$xlims[[1]] <- -0.1
    if (out$xlims[[2]] <= 0) out$xlims[[2]] <- 0.1
    
    return(out)
}


prep_volcano_data <- function(
         stats_tbl,
         current_cells,
         signif_color = "black",
         nonsignif_color = "grey",
         # column names
         cell_type_col = "cell_types",
         foldchange_col = "log2FC",
         p_value_col = "p_adj",
         # significance
         signif_threshold = 0.05,
         # label text colour function
         label_text_fun = IMCfuncs::get_text_color,
         # axis expansion
         x_allowance = 0.1,
         y_allowance = 0,
         # title
         plot_title = NULL,
         # user-supplied method info (named list)
         method_info = list(),
         # OPTIONAL overrides for colour/shape
         colour_by = NULL,        # column name or vector
         shape_by = NULL,         # column name or vector
         # optional explicit mappings (used when colour_by/shape_by is a column)
         colour_values = NULL,    # named or unnamed vector of colours
         shape_values = NULL,     # named or unnamed vector of shapes
         # EnhancedVolcano visual parameters
         FCcutoff = 0,
         point_size = 7,
         lab_size = 6.0,
         legend_position = "top",
         legend_lab_size = 15,
         legend_icon_size = 5.0,
         boxed_labels = TRUE,
         draw_connectors = TRUE,
         arrowheads = FALSE,
         lab_col = "black",
         lab_face = "bold",
         legend_color_title = "foo",
         legend_shape_title = "bar") {
        
    # --- basic checks ----------------------------------------------------
    needed_cols <- c(cell_type_col, foldchange_col, p_value_col)
    missing_cols <- setdiff(needed_cols, colnames(stats_tbl))
    if (length(missing_cols) > 0) {
        stop(
            "prep_volcano_data: missing required columns in stats_tbl: ",
            paste(missing_cols, collapse = ", ")
        )
    }
    
    if (!all(stats_tbl[[cell_type_col]] %in% names(current_cells))) {
        missing_cells = setdiff(unique(stats_tbl[[cell_type_col]]), names(current_cells))
        warning(
            "Some cell types have no colour in current_cells: ",
            paste(missing_cells, collapse = ", ")
        )
    }
    
    n <-  nrow(stats_tbl)
    
    # --- build plot_data -------------------------------------------------
    plot_data <- stats_tbl
    plot_data$label_fill <- current_cells[plot_data[[cell_type_col]]]
    plot_data$label_text <- label_text_fun(plot_data$label_fill)
    
    # default significance: always from numeric p_value_col
    plot_data$is_significant <-  plot_data[[p_value_col]] < signif_threshold
    
    # default colours based on significance
    plot_data$point_fill <-  ifelse(
        test = plot_data$is_significant,
        yes = signif_color,
        no = nonsignif_color
    )
    
    # default shapes: leave NULL (EnhancedVolcano default)
    plot_data$point_shape <- NULL
    
    # --- optional COLOUR override ---------------------------------------
    if (!is.null(colour_by)) {
        if (is.character(colour_by) && length(colour_by) == 1L && colour_by %in% colnames(stats_tbl)) {
            
            group <- stats_tbl[[colour_by]]
            group_chr <- as.character(group)
            group_levels <- unique(group_chr)
            
            # derive a mapping from group -> colour
            if (!is.null(colour_values)) {
                # if named, match by name; if unnamed, recycle over levels
                if (!is.null(names(colour_values))) {
                    missing_map <- setdiff(group_levels, names(colour_values))
                    if (length(missing_map) > 0) {
                        stop("colour_values is missing entries for: ",
                             paste(missing_map, collapse = ", "))
                    }
                    colour_map <- colour_values[group_levels]
                } else {
                    # unnamed: recycle / slice
                    colour_map <- rep(colour_values, length.out = length(group_levels))
                    names(colour_map) <- group_levels
                }
            } else {
                # no colour_values given: generate a simple palette
                colour_map <- grDevices::rainbow(length(group_levels))
                names(colour_map) <- group_levels
            }
            
            plot_data$point_fill <- colour_map[group_chr]
            
        } else if (length(colour_by) == n) {
            # direct per-point vector of colours
            plot_data$point_fill <- as.character(colour_by)
        } else {
            stop("colour_by must be either a column name in stats_tbl or a vector of length nrow(stats_tbl).")
        }
    }
    # --- optional SHAPE override ----------------------------------------
    if (!is.null(shape_by)) {
        if (is.character(shape_by) && length(shape_by) == 1L && shape_by %in% colnames(stats_tbl)) {
            
            group <- stats_tbl[[shape_by]]
            group_chr <- as.character(group)
            group_levels <- unique(group_chr)
            
            if (!is.null(shape_values)) {
                if (!is.null(names(shape_values))) {
                    missing_shapes <- setdiff(group_levels, names(shape_values))
                    if (length(missing_shapes) > 0) {
                        stop("shape_values is missing entries for: ",
                             paste(missing_shapes, collapse = ", "))
                    }
                    shape_map <- shape_values[group_levels]
                } else {
                    # unnamed: recycle over levels
                    default_shapes <- shape_values
                    shape_map <- rep(default_shapes, length.out = length(group_levels))
                    names(shape_map) <- group_levels
                }
            } else {
                # default shape set (cycles if more groups)
                default_shapes <- c(21, 24, 22, 25, 23)
                shape_map <- rep(default_shapes, length.out = length(group_levels))
                names(shape_map) <- group_levels
            }
            
            plot_data$point_shape = shape_map[group_chr]
            
        } else if (length(shape_by) == n) {
            # direct per-point vector of shapes
            plot_data$point_shape <- shape_by
        } else {
            stop("shape_by must be either a column name in stats_tbl or a vector of length nrow(stats_tbl).")
        }
    }
    # --- output list -----------------------------------------------------
    out <- list()
    out$plot_data <- plot_data
    
    # column names to reuse later
    out$cell_type_col <- cell_type_col
    out$foldchange_col <- foldchange_col
    out$p_signif_col <- p_value_col
    
    out$title = plot_title
    
    # y-axis label expression
    is_adj <- grepl("adj", p_value_col, ignore.case = TRUE)
    if (is_adj) {
        out$ylab_exprs <- bquote(~ -Log[10] ~ italic(P[Adjusted]))
    } else {
        out$ylab_exprs <- bquote(~ -Log[10] ~ italic(P))
    }
    
    # method_info: keep user bits, then append p metadata at bottom
    if (is.null(method_info) || !is.list(method_info)) {
        method_info <- list()
    }
    
    method_info[c("p_thresh", "p_adjusted")] <- NULL
    method_info <- c(
        method_info,
        list(
            p_thresh = signif_threshold,
            p_adjusted = is_adj
        )
    )
    
    out$method_info <- method_info
    
    # colCustom: per-point colours; names used only for grouping in legend
    out$point_fill <- stats::setNames(
        object = plot_data$point_fill,
        nm = ifelse(plot_data$is_significant, "significant", "not significant")
    )
    
    # shapeCustom: per-point shapes (can be NULL)
    out$point_shape <- plot_data$point_shape
    
    # axis limits via set_axis_lims()
    lims <- set_axis_lims(
        df          = plot_data,
        y_var       = foldchange_col,
        x_var       = p_value_col,
        x_allowance = x_allowance,
        y_allowance = y_allowance
    )
    
    out$xlim <- lims$xlims
    out$ylim <- lims$ylims
    
    # EnhancedVolcano control parameters
    out$pCutoff <- signif_threshold
    out$FCcutoff <- FCcutoff
    out$pointSize <- point_size
    out$labSize <- lab_size
    out$legendPosition <- legend_position
    out$legendLabSize <- legend_lab_size
    out$legendIconSize <- legend_icon_size
    out$boxedLabels <- boxed_labels
    out$drawConnectors <- draw_connectors
    out$arrowheads <- arrowheads
    out$labCol <- lab_col
    out$labFace <- lab_face
    out$legend_color_title <- legend_color_title
    out$legend_shape_title <- legend_shape_title
    
    return(out)
    }

plot_volcano <- function(prep_list) {
    vars_col <- prep_list$cell_type_col
    foldchange_col <- prep_list$foldchange_col
    
    if (!vars_col %in% colnames(prep_list$plot_data)){
        stop("vars_col supplied is not in prep_list$plot_data")
    }
    
    subtitle_str <- ""
    if (!is.null(prep_list$method_info) && length(prep_list$method_info) > 0) {
        kv <- paste(
            paste(names(prep_list$method_info), prep_list$method_info, sep = " = "),
            collapse = "\n"
        )
        subtitle_str = kv
    }
    
    plot_labels <- list(
        title = if (is.null(prep_list$title)) "" else prep_list$title,
        subtitle = subtitle_str,
        caption = paste("Total cell types:", nrow(prep_list$plot_data))
    )
    
    ylab_exprs <- prep_list$ylab_exprs
    
    p <- EnhancedVolcano::EnhancedVolcano(
        toptable = prep_list$plot_data,
        lab = prep_list$plot_data[[vars_col]],
        x = foldchange_col,
        y = prep_list$p_signif_col,
        pCutoff = prep_list$pCutoff,
        FCcutoff = prep_list$FCcutoff,
        pointSize = prep_list$pointSize,
        labSize = prep_list$labSize,
        colCustom = prep_list$point_fill,
        shapeCustom = prep_list$point_shape,
        legendPosition = prep_list$legendPosition,
        legendLabSize = prep_list$legendLabSize,
        legendIconSize = prep_list$legendIconSize,
        boxedLabels = prep_list$boxedLabels,
        drawConnectors = prep_list$drawConnectors,
        arrowheads = prep_list$arrowheads,
        labCol = prep_list$labCol,
        labFace = prep_list$labFace,
        title = plot_labels$title,
        subtitle = plot_labels$subtitle,
        caption = plot_labels$caption,
        ylab = ylab_exprs,
        xlim = prep_list$xlim,
        ylim = prep_list$ylim
    )
    
    # styling tweaks
    p$layers[[1]]$aes_params$color <- "black"
    p$layers[[1]]$aes_params$fill <- prep_list$point_fill
    p$layers[[1]]$aes_params$alpha <- 0.6
    
    sig_idx <- prep_list$plot_data$is_significant
    p$layers[[4]]$aes_params$segment.linetype <- 2
    p$layers[[4]]$aes_params$segment.alpha <- 0.5
    p$layers[[4]]$aes_params$fill <- prep_list$plot_data$label_fill[sig_idx]
    p$layers[[4]]$aes_params$colour <- prep_list$plot_data$label_text[sig_idx]
    p$layers[[4]]$geom_params$box.padding <- grid::unit(1, "lines")
    
    p <- p + guides(
        color = guide_legend(reverse = TRUE, title = prep_list$legend_color_title, title.position = "top"),
        shape = guide_legend(reverse = TRUE, title = prep_list$legend_shape_title, title.position = "top")
    )
    
    p <- p + theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text( face = "italic", size = 14, margin = margin(b = 10))
    )
    
    return(p)
}