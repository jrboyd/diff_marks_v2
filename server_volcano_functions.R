
# plot style functions
ref_line = function(MIN, MAX) lines(c(MIN, MAX), c(MIN, MAX), col = rgb(0, 0, 1, 1), lty = 3, lwd = 6)

custom_plot = function(x, y, col, txt, ...) {
    plot(x = x, y = y, col = col, pch = 16, ...)
    MIN = par("usr")[1]
    MAX = par("usr")[2]
    box()
    par(family = "mono")
    text(MIN + (MAX - MIN) * 0.02, MAX - (MAX - MIN) * 0.01, txt, adj = c(0, 1))
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
    toSort = rep(0, length(colors))
    names(toSort) = rownames(data)
    toSort[list_a] = colors[list_a]
    toSort[list_b] = colors[list_b]
    o = order(toSort, decreasing = T)
    # joined plot
    custom_plot(data[o, a], data[o, b], colors[o], txt = note, ...)
}
