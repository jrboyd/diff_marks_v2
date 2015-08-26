
# plot style functions
ref_line = function(MIN, MAX) lines(c(0, 0), c(MIN, MAX), col = "dodgerblue", lty = 5, lwd = 2)

custom_plot = function(x, y, col, txt, ...) {
    x_plot = y - x
    y_plot = apply(cbind(x, y), 1, min)
    plot(x = x_plot, y = y_plot, col = col, pch = 16, ...)
    MIN = par("usr")[1]
    MAX = par("usr")[2]
    box()
    par(family = "mono")
    text(MIN + (MAX - MIN) * 0.02, MAX - (MAX - MIN) * 0.05, txt, adj = c(0, 1))
    ref_line(MIN, MAX)
    par(family = "")
}

scale_colors = function(data, scale, list_in, bg_color = "black", list_color = "red", colors = NA) {
    # input colors if modifying existing colors list
    if (is.na(colors[1])) {
        colors = rep(rgb(0, 0, 0, 0.1), nrow(data))
        names(colors) = rownames(data)
    }
    if (is.na(scale[1])) {
        scale = rep(1, nrow(data))
        names(scale) = rownames(data)
    }
    if(length(intersect(rownames(data), list_in)) == 0) return(colors)
    cr = colorRamp(colors = c(bg_color, list_color), alpha = T)
    tmp = cr(scale[list_in])/255
    colors[list_in] = rgb(tmp[, 1], tmp[, 2], tmp[, 3], tmp[, 4])
    return(colors)
}

plot_list = function(data, list_in, colors, a, b, note = "", scale = NA, ...) {
    custom_plot(data[list_in, a], data[list_in, b], colors[list_in], txt = note, cex = 0.5, ...)
}

plot_bg = function(data, lists_excluded, colors, a, b, note = "", scale = NA, ...) {
    non_ensgs = rownames(data)
    for (l in lists_excluded) {
        non_ensgs = setdiff(non_ensgs, l)
    }
    plot_list(data, non_ensgs, colors, a, b, note = note, scale = scale, ...)
}

plot_merge = function(data, list_a, list_b, colors, a = 1, b = 2, note, ...) {
    colors = colors[rownames(data)]
    toSort = rep(0, length(colors))
    names(toSort) = rownames(data)
    list_a = intersect(list_a, names(colors))
    list_b = intersect(list_b, names(colors))
    toSort[list_a] = colors[list_a]
    toSort[list_b] = colors[list_b]
    o = order(toSort, decreasing = T)
    # joined plot
    custom_plot(data[o, a], data[o, b], colors[o], txt = note, ...)
}
