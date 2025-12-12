# labels: character vector of label names
# label_to_category: named vector mapping each label -> category  (names = labels)
# palette: one of "turbo","viridis","magma","plasma","inferno","cividis","mako","rocket"
# alpha, direction: passed straight through to the viridisLite function
# return_outlines: if TRUE, also returns complementary (mirrored) colours for outlines
splice_a_viridis = function(labels,
                            label_to_category,
                            palette = "turbo",
                            alpha = 1,
                            direction = 1,
                            return_outlines = TRUE) {
    # sanity
    labels = as.character(labels)
    if (!all(labels %in% names(label_to_category))) {
        miss = setdiff(labels, names(label_to_category))
        stop("Missing label_to_category entries for: ", paste(miss, collapse = ", "))
    }
    
    # choose the viridisLite function
    pal_fun = switch(
        palette,
        "turbo"   = viridisLite::turbo,
        "viridis" = viridisLite::viridis,
        "magma"   = viridisLite::magma,
        "plasma"  = viridisLite::plasma,
        "inferno" = viridisLite::inferno,
        "cividis" = viridisLite::cividis,
        "mako"    = viridisLite::mako,
        "rocket"  = viridisLite::rocket,
        stop("Unknown palette: ", palette)
    )
    
    # categories follow first appearance order (change to sort(unique(...)) if you prefer)
    cats_seen = label_to_category[labels]
    categories = unique(cats_seen)
    ncat = length(categories)
    
    # equal breakpoints across [0,1] (no margins)
    bps = seq(0, 1, length.out = ncat + 1)
    
    fills = character(length(labels))
    names(fills) = labels
    outlines = if (return_outlines) fills else NULL
    
    for (i in seq_len(ncat)) {
        cat = categories[i]
        labs_in_cat = labels[cats_seen == cat]
        k = length(labs_in_cat)
        
        # fill colours: evenly spaced inside this category’s segment
        fill_cols = pal_fun(
            n = k,
            alpha = alpha,
            begin = bps[i],
            end = bps[i + 1],
            direction = direction
        )
        
        # complementary/contrasting outlines:
        if (return_outlines) {
            # mirror the segment: t -> 1 - t
            outline_cols = pal_fun(
                n = k,
                alpha = alpha,
                begin = 1 - bps[i + 1],
                end = 1 - bps[i],
                direction = direction
            )
        }
        
        # deterministic assignment (alphabetical within category)
        ord = order(labs_in_cat)
        labs_sorted = labs_in_cat[ord]
        fills[labs_sorted] = fill_cols[seq_len(k)]
        if (return_outlines) outlines[labs_sorted] = outline_cols[seq_len(k)]
    }
    
    if (return_outlines) list(fills = fills, outlines = outlines) else fills
}




# pick black/white text for readability on a given hex fill
.best_text_for = function(hex) {
    rgb = grDevices::col2rgb(hex) / 255
    # relative luminance
    L = 0.2126*rgb[1,] + 0.7152*rgb[2,] + 0.0722*rgb[3,]
    ifelse(L > 0.55, "#000000", "#FFFFFF")
}

# pal can be:
#   1) named vector of hex (fills), or
#   2) list(fills=<named hex>, outlines=<named hex>)


plot_palette_swatches = function(pal, ncol = 6, title = "Palette preview",
                                 show_labels = TRUE, label_size = 3.2, border_width = 2) {
    
    if (is.list(pal)) {
        fills = pal$fills
        outlines = pal$outlines
    } else {
        fills = pal
        outlines = rep("#000000", length(fills))
        names(outlines) = names(fills)
    }
    
    stopifnot(!is.null(names(fills)))
    labs = names(fills)
    
    df = data.frame(
        label = labs,
        fill  = unname(fills),
        outline = unname(outlines[labs]),
        idx = seq_along(labs)
    )
    
    df$col = ((df$idx - 1) %% ncol) + 1
    df$row = -((df$idx - 1) %/% ncol)  # negative so first row is at the top
    df$text_col = .best_text_for(df$fill)
    
    p = ggplot(df, aes(x = col, y = row)) +
        geom_tile(aes(fill = fill), color = df$outline, linewidth = border_width, width = 0.85, height = 0.85) +
        scale_fill_identity() +
        coord_equal(expand = FALSE) +
        theme_void() +
        ggtitle(paste0(title, "\n")) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.margin = margin(8, 8, 8, 8)
        )
    
    if (show_labels) {
        p = p + geom_text(aes(label = label), color = df$text_col, size = label_size)
    }
    
    p
}