assign_ensg = function(test_res, ref_dict) {
    # columns must be chrm, start, end, ref_dict rownames are ids sort all by start then by chrm
    o = order(ref_dict[, 2])
    ref_dict = ref_dict[o, ]
    o = order(ref_dict[, 1])
    ref_dict = ref_dict[o, ]
    
    o = order(test_res[, 2])
    test_res = test_res[o, ]
    o = order(test_res[, 1])
    test_res = test_res[o, ]
    
    keep_res = rep(F, nrow(test_res))
    ensg_res = rep("", nrow(test_res))
    # intersect ma_res with ensgs, assigning ensg or removing peaks
    j = 1  #track position in ref_dict
    last_chrm = ""
    for (i in 1:nrow(test_res)) {
        # iterate test_res
        test_chrm = test_res[i, 1]
        if (test_chrm != last_chrm) {
            last_chrm = test_chrm
            print(last_chrm)
        }
        test_start = test_res[i, 2]
        test_end = test_res[i, 3]
        checking_ref = T
        while (checking_ref) {
            if (j > nrow(ref_dict)) {
                checking_ref = F
                next
            }
            ref_ensg = rownames(ref_dict)[j]
            ref_chrm = ref_dict[j, 1]
            ref_start = ref_dict[j, 2]
            ref_end = ref_dict[j, 3]
            if (test_chrm != ref_chrm) {
                # on diff chrm
                o = order(c(ref_chrm, test_chrm))  #not clear if ref or test has gotten ahead
                if (o[1] == 2) {
                  checking_ref = F
                  
                } else {
                  j = j + 1
                }
                # print(paste('ref', ref_chrm, 'not queried', test_chrm)) print(paste('i', i, 'j', j)) print(paste(ref_chrm, 'to', ref_dict[j + 1,5]))
                next
            }
            last = j
            if (test_end < ref_start) {
                # peak still before next ref feature, check next peak print(paste('next peak:',test_end, '<', ref_start))
                checking_ref = F
                next
            }
            # test may now intersect ref or be past ref there is overlap, print(paste('MATCH:',test_start, '<', ref_end))
            if (test_start < ref_end) {
                keep_res[i] = T
                ensg_res[i] = ref_ensg
                checking_ref = F  #done with this test entry
                j = j + 1  #we don't want to use ref entries more than once
                next
            }
            # test start must be past ref_end, advance ref, keep checking peak print('next ref')
            j = j + 1
            # checking_ref = F
        }
    }
    test_res = test_res[keep_res, ]
    rownames(test_res) = ensg_res[keep_res]
    return(test_res)
}

get_ensgs = function(res, is1, dat_names) {
    id_col = 4
    keep = logical()
    if (is1) {
        keep = res[, id_col] == "cond1"
    } else {
        keep = res[, id_col] == "cond2"
    }
    ensgs_a = unique(rownames(res)[keep])
    ensgs_a = intersect(ensgs_a, dat_names)
    return(ensgs_a)
}



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

plot_merge = function(data, list_a, list_b, colors, a, b, note, ...) {
    toSort = rep(0, length(colors))
    names(toSort) = rownames(data)
    toSort[list_a] = colors[list_a]
    toSort[list_b] = colors[list_b]
    o = order(toSort, decreasing = T)
    # joined plot
    custom_plot(data[o, a], data[o, b], colors[o], txt = note, ...)
}

plot_2lists = function(data, list_a, list_b, a, b, scale = NA, plot_title = "TITLE") {
    # organize plot
    layout(matrix(c(rep(1, 3), 2:4, rep(5, 3)), ncol = 3, byrow = T), heights = c(1, 1, 3))
    par(mai = rep(0, 4), xpd = F)
    
    plot(1, type = "n", axes = F, xlim = c(0, 1), ylim = c(0, 1))
    text(0.5, 0.1, plot_title, cex = 3, adj = c(0.5, 0))
    
    # establish parameters
    MAX = max(data)
    MIN = min(data)
    
    colors = scale_colors(data, scale, list_a)
    colors = scale_colors(data, scale, list_b, list_color = "green", colors = colors)
    
    # ensures higher values in list are drawn on top
    o = order(scale[list_a])
    list_a = list_a[o]
    
    o = order(scale[list_b])
    list_b = list_b[o]
    
    cell_lines = matrix(unlist(strsplit(colnames(data), split = " ")), nrow = 2)[1, ]
    
    # list_a plot
    note = paste0("up in ", cell_lines[a], " : ", length(list_a))
    plot_list(data, list_a, colors, a, b, note, scale, xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F)
    note = paste0("up in ", cell_lines[b], " : ", length(list_b))
    plot_list(data, list_b, colors, a, b, note, scale, xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F)
    
    note = paste0("background", ":", (nrow(data) - length(list_a) - length(list_b)))
    plot_bg(data, list(list_a, list_b), colors, a, b, note, scale, xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F)
    
    note = paste0("merge", ":", nrow(data))
    plot_merge(data, list_a, list_b, colors, a, b, note, xlim = c(MIN, MAX), ylim = c(MIN, MAX), axes = F, cex = 1.5)
    # toSort = rep(0, length(colors)) names(toSort) = rownames(data) toSort[list_a] = scale[list_a] toSort[list_b] = scale[list_b] o = order(toSort, decreasing = T)
    # #joined plot custom_plot(data[o,a], data[o,b], colors[o], txt = paste0('merge', ':', nrow(data)), MIN, MAX, cex = 1.5)
    
}

passing_fc = function(a, b, data, floor = 0, thresh = 2) {
    data = ifelse(data < floor, floor, data)
    keep = (data[, a] - data[, b]) > thresh
    return(rownames(data)[keep])
} 
