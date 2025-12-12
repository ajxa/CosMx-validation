cor_mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
}

adjust_p_matrix <- function(p_mat, method = "fdr", diag_value = 1) {
    
    p_mat <- as.matrix(p_mat)
    if (nrow(p_mat) != ncol(p_mat)) {
        stop("p_mat must be a square matrix.")
    }
    
    n <- nrow(p_mat)
    
    # Indices for upper triangle (excluding diagonal)
    upper_idx <- upper.tri(p_mat, diag = FALSE)
    
    pvals <- p_mat[upper_idx]
    
    padj <- stats::p.adjust(pvals, method = method)
    
    padj_mat <- matrix(NA_real_, nrow = n, ncol = n)
    dimnames(padj_mat) <- dimnames(p_mat)
    
    # Fill upper triangle
    padj_mat[upper_idx] <- padj
    
    # Mirror to lower triangle to keep symmetry
    padj_mat[lower.tri(padj_mat, diag = FALSE)] <-
        t(padj_mat)[lower.tri(padj_mat, diag = FALSE)]
    
    diag(padj_mat) <- diag_value
    
    return(padj_mat)
}

clean_corr_plot_data <- 
    function(
        corr_matrix,
        p_matrix,
        p_sig_level = 0.05,
        minimum_cor = 0.5,
        type = c("full", "lower", "upper"),
        remove_diagonal = FALSE,
        show_below_min_cor = FALSE,
        digits = 2
        ) {
    if (!is.matrix(corr_matrix)) stop("corr_matrix is not a matrix!")
    if (!is.matrix(p_matrix)) stop("p_matrix is not a matrix!")
    if (nrow(corr_matrix) != ncol(corr_matrix)) stop("corr_matrix is not a square matrix!")
    if (nrow(p_matrix) != ncol(p_matrix)) stop("p_matrix is not a square matrix!")
    if (nrow(corr_matrix) != nrow(p_matrix)) stop("corr_matrix and p_matrix have different dimensions!")
    
    type <- match.arg(type, several.ok = FALSE, choices = c("full", "lower", "upper"))
    
    corr_matrix <- base::round(x = corr_matrix, digits = digits)
    
    if (type == "lower") {
        corr_matrix[upper.tri(corr_matrix)] <- NA
        p_matrix[upper.tri(p_matrix)] <- NA
    } else if (type == "upper") {
        corr_matrix[lower.tri(corr_matrix)] <- NA
        p_matrix[lower.tri(p_matrix)] <- NA
    }
    
    if (remove_diagonal) {
        diag(corr_matrix) <- NA
        diag(p_matrix) <- NA
    }
    
    corr_matrix <- as.data.frame.table(corr_matrix, responseName = "value")
    p_matrix <- as.data.frame.table(p_matrix, responseName = "value")
    
    corr_matrix$pvalue <- p_matrix$value
    corr_matrix$p_signif <- corr_matrix$pvalue < p_sig_level
    corr_matrix$min_cor <- abs(corr_matrix$value) >= minimum_cor
    
    if (!show_below_min_cor) {
        corr_matrix$plot_value <- ifelse(corr_matrix$p_signif & corr_matrix$min_cor, corr_matrix$value, NA)
    }
    
    return(corr_matrix)
}

plot_corr <- function(corr_plot_data,
                      show_labels = TRUE,
                      label_colour = "black",
                      label_size = 4,
                      shape_lab_col = NULL,
                      shape_lab_shape = 23,
                      shape_lab_size = 5,
                      shape_lab_colour = "black",
                      shape_lab_fill = "grey",
                      tile_outline_colour = "grey25",
                      tile_linetype = 3,
                      title = "",
                      subtitle = "",
                      diverging_colours = io$inputs$colors$diverging,
                      pos_colours = io$inputs$colors$positive,
                      colour_bar_name = "Corr",
                      colour_barwidth = 2.5,
                      colour_barheight = 15,
                      colour_bar_ouline = "grey90",
                      legend_pos = c("right", "left", "top", "bottom")) {
    # Default plot theme:
    corr_plot_theme <-
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text = element_text(size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
            axis.text.y = element_text(size = 12),
            axis.ticks = element_line(color = "grey50", linewidth = 0.2),
            axis.title = element_blank(),
            legend.title = element_blank(),
            legend.key.size = ggplot2::unit(10, "mm"),
            legend.text = element_text(size = 12),
            legend.position = match.arg(legend_pos, several.ok = FALSE, choices = c("right", "left", "top", "bottom")),
            panel.grid.major = element_blank(),
            plot.title = element_text(face = "bold", size = 20),
            plot.subtitle = element_text(size = 12)
        )
    
    # Main plot code:
    
    if (!is.null(shape_lab_col)) {
        if (!shape_lab_col %in% colnames(corr_plot_data)) stop("label_column is not present in the data")
        show_labels <- FALSE
        shape_labels <- corr_plot_data[[shape_lab_col]]
    }
    
    p <- corr_plot_data %>%
        ggplot(mapping = aes(x = Var1, y = Var2)) +
        geom_tile(aes(fill = plot_value), color = tile_outline_colour, linetype = tile_linetype, na.rm = TRUE) +
        {
            if (show_labels) geom_text(aes(label = plot_value), color = label_colour, size = label_size, na.rm = TRUE)
        } +
        {
            if (!is.null(shape_lab_col)) geom_point(aes(shape = !!sym(shape_lab_col)), size = shape_lab_size, color = shape_lab_colour, fill = shape_lab_fill, na.rm = TRUE)
        } +
        labs(title = title, subtitle = subtitle)
    
    
    if (all(corr_plot_data$plot_value > 0, na.rm = TRUE)) {
        p <- p + ggplot2::scale_fill_gradientn(
            colors = pos_colours,
            guide = guide_colorbar(
                ticks = TRUE,
                ticks.colour = colour_bar_ouline,
                frame.colour = colour_bar_ouline,
                barwidth = colour_barwidth,
                barheight = colour_barheight
            ),
            name = colour_bar_name,
            na.value = "white"
        )
    } else {
        p <- p + ggplot2::scale_fill_gradientn(
            colors = diverging_colours,
            limits = c(-1, 1),
            guide = guide_colorbar(
                ticks = TRUE,
                ticks.colour = colour_bar_ouline,
                frame.colour = colour_bar_ouline,
                barwidth = colour_barwidth,
                barheight = colour_barheight
            ),
            name = colour_bar_name,
            na.value = "white"
        )
    }
    
    if (!is.null(shape_lab_col)) p <- p + scale_shape_manual(values = c(shape_lab_shape), na.translate = FALSE)
    
    p <- p + corr_plot_theme
    
    return(p)
}


add_highlight_regions <- function(baseplot,
                                  corr_matrix,
                                  groups,
                                  highlight_colors,
                                  highlight_width = 3) {
    unique_groups <- unique(groups)
    
    highlight_regions <- data.frame(
        group = factor(unique_groups, levels = names(highlight_colors)),
        xmin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        xmax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        ymin = sapply(unique_groups, function(g) min(which(groups == g))) - 0.5,
        ymax = sapply(unique_groups, function(g) max(which(groups == g))) + 0.5,
        color = highlight_colors[unique(groups)]
    )
    
    out_plot <- baseplot +
        geom_rect(
            data = highlight_regions,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = group),
            fill = NA, linewidth = highlight_width, inherit.aes = FALSE
        ) +
        scale_color_manual(name = "Groups", values = highlight_colors)
    
    return(out_plot)
}
