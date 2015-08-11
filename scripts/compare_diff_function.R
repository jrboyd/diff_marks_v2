ma_dp_fc_compare = function(a, b, filter_detected = F, file_prefix = NA, p_thresh = 2) {
    print("constructing reference dictionary...")
    ensg_dict = parse_gtf("ref/ensg_promoters_counted.gff")
    keep = ensg_dict[, 5] != "chrM"
    ensg_dict = ensg_dict[keep, ]
    
    lines = c("MCF10A", "MCF7", "MDA231")
    lines = rep(lines, 2)
    mods = c("H3K4AC", "H3K4ME3")
    l1 = lines[a]
    l2 = lines[b]
    m = mods[1]
    if (a > 3 & b > 3) {
        m = mods[2]
    } else if (!(a <= 3 & b <= 3)) {
        stop("histone marks don't match!")
    }
    # columns of res tables will be chr,start,end,up/down,log10pval,...
    
    fname = dir(path = "ma_norm_results", full.names = T, pattern = paste0(paste(l1, m, sep = "_"), ".+", paste(l2, m, sep = "_"), ".+result.xls"))
    print(paste(fname, "from MAnorm loading..."))
    ma_res = read.table(fname, header = T, sep = "\t", comment.char = "", stringsAsFactors = F)
    ma_res = ma_res[, c(1:4, 9)]
    
    keep = !grepl(pattern = "common", ma_res[, 4])  #remove common peaks, only interesetd in diff
    ma_res = ma_res[keep, ]
    
    keep = grepl(pattern = "unique_peak1", ma_res[, 4])  #unify terminology for convenience
    ma_res[keep, 4] = "cond1"
    
    keep = grepl(pattern = "unique_peak2", ma_res[, 4])
    ma_res[keep, 4] = "cond2"
    
    
    # keep = grepl(pattern = 'unique_peak1', ma_res[, 4]) ma_res[keep, 4] = paste(l1, 'higher than', l2, 'for', m) keep = grepl(pattern = 'unique_peak2', ma_res[, 4])
    # ma_res[keep, 4] = paste(l2, 'higher than', l1, 'for', m) keep = grepl(pattern = 'common_peak1', ma_res[, 4]) ma_res[keep, 4] = paste(l1, 'shared with', l2, 'for',
    # m) keep = grepl(pattern = 'common_peak2', ma_res[, 4]) ma_res[keep, 4] = paste(l2, 'shared with', l1, 'for', m)
    
    diffpeak_res = matrix("", nrow = 0, ncol = 5)
    bed_files = dir(path = "diff_peak_results", full.names = T, pattern = paste0(paste(l1, m, sep = "_"), ".+", paste(l2, m, sep = "_"), ".+.bed"))
    
    for (fname in bed_files) {
        print(paste(fname, " from MACS2 bdgdiff loading..."))
        tmp = read.table(fname, header = F, skip = 1, stringsAsFactors = F)
        diffpeak_res = rbind(diffpeak_res, tmp)
    }
    keep = !grepl(pattern = "common", diffpeak_res[, 4])
    diffpeak_res = diffpeak_res[keep, ]
    
    keep = grepl(pattern = "cond1", diffpeak_res[, 4])
    diffpeak_res[keep, 4] = "cond1"
    
    keep = grepl(pattern = "cond2", diffpeak_res[, 4])
    diffpeak_res[keep, 4] = "cond2"
    
    
    # keep = grepl(pattern = 'common', diffpeak_res[, 4]) diffpeak_res[keep, 4] = paste(l1, 'same as', l2, 'for', m) keep = grepl(pattern = 'cond1', diffpeak_res[, 4])
    # diffpeak_res[keep, 4] = paste(l1, 'higher than', l2, 'for', m) keep = grepl(pattern = 'cond2', diffpeak_res[, 4]) diffpeak_res[keep, 4] = paste(l1, 'lower than',
    # l2, 'for', m)
    
    library("gtools")
    
    
    ref = ensg_dict[, c("chrm", "start", "end")]
    
    print("Matching MAnorm results to promoters...")
    ma_res = assign_ensg(ma_res, ref)
    print("Matching MACS2 bdgdiff results to promoters...")
    diffpeak_res = assign_ensg(diffpeak_res, ref)
    
    MIN_PVAL = p_thresh
    MAX_PVAL = 100
    is_na = is.na(ma_res[, 5])  #after consulting tracks, NA ma_norm diff peaks are exclusive to one cell line, should be high p-value
    is_na = rownames(ma_res)[is_na]
    ma_res[is_na, 5] = MAX_PVAL
    is_inf = is.infinite(ma_res[, 5])  #inf values should also be high
    is_inf = rownames(ma_res)[is_inf]
    ma_res[is_inf, 5] = MAX_PVAL
    
    # filter for pval thresh
    keep = ma_res[, 5] > p_thresh
    ma_res = ma_res[keep, ]
    
    keep = diffpeak_res[, 5] > p_thresh
    diffpeak_res = diffpeak_res[keep, ]
    
    
    # bdgdiff does not have NA or inf values
    
    load("my_fe_corrected.save")
    
    if (filter_detected) {
        keep = apply(my_fe, 1, max) > 1.5
        my_fe = my_fe[keep, ]
    }
    
    print("Plotting comparisons...")
    ensgs_ma_b_gt_a = get_ensgs(ma_res, is1 = F, rownames(my_fe))
    ensgs_ma_a_gt_b = get_ensgs(ma_res, is1 = T, rownames(my_fe))
    ensgs_dp_b_gt_a = get_ensgs(diffpeak_res, is1 = F, rownames(my_fe))
    ensgs_dp_a_gt_b = get_ensgs(diffpeak_res, is1 = T, rownames(my_fe))
    
    res2pval = function(res) {
        pval = res[, 5]
        names(pval) = rownames(res)
        pval = ifelse(pval > MAX_PVAL, MAX_PVAL, pval)
        pval = ifelse(pval < MIN_PVAL, MIN_PVAL, pval)
        pval = pval * m_scale + b_scale
        return(pval)
    }
    # convert pvalues to scale for plotting
    MIN_SCALE = 0.2
    MAX_SCALE = 1
    
    m_scale = (MAX_SCALE - MIN_SCALE)/(MAX_PVAL - MIN_PVAL)
    b_scale = MIN_SCALE - MIN_PVAL * m_scale
    ma_pval = res2pval(ma_res)
    dp_pval = res2pval(diffpeak_res)
    
    
    
    plot_2lists_wrap = function(list_a, list_b, scale = NA, plot_count = 0, plot_title = "TITLE") {
        # , scale_a = NA, scale_b = NA){ i'm assuming both scales will be set or neither
        if (length(scale) == 1 && is.na(scale)) {
            scale = rep(1, nrow(my_fe))
            names(scale) = rownames(my_fe)
        }
        
        if (!is.na(file_prefix)) 
            png(paste0(file_prefix, "_", plot_count, ".png"), width = 6, height = 10, res = 150, units = "in")
        plot_2lists(my_fe, list_a, list_b, a, b, scale, plot_title)
        if (!is.na(file_prefix)) 
            dev.off()
    }
    
    
    plot_n = 0
    plot_n = plot_n + 1
    plot_2lists_wrap(ensgs_dp_a_gt_b, ensgs_dp_b_gt_a, dp_pval, plot_count = plot_n, plot_title = "MACS2 bdgdiff")
    
    # plot_n = plot_n + 1 plot_2lists_wrap(ensgs_ma_a_gt_b, ensgs_ma_b_gt_a, ma_pval, plot_count = plot_n, plot_title = 'MAnorm') ensgs_both_a_gt_b =
    # intersect(ensgs_ma_a_gt_b, ensgs_dp_a_gt_b) ensgs_both_b_gt_a = intersect(ensgs_ma_b_gt_a, ensgs_dp_b_gt_a) comb = c(ensgs_both_a_gt_b, ensgs_both_b_gt_a)
    # comb_pval = (ma_pval[comb] + dp_pval[comb]) / 2 plot_n = plot_n + 1 plot_2lists_wrap(ensgs_both_a_gt_b, ensgs_both_b_gt_a, comb_pval, plot_count = plot_n,
    # plot_title = 'MACS2 bdgdiff\nMAnorm') ensgs_fc_a_gt_b = passing_fc(a, b, my_fe) ensgs_fc_b_gt_a = passing_fc(b, a, my_fe) plot_n = plot_n + 1
    # plot_2lists_wrap(ensgs_fc_a_gt_b, ensgs_fc_b_gt_a, plot_count = plot_n, plot_title = '4 FC') ensgs_all_a_gt_b = intersect(ensgs_both_a_gt_b, ensgs_fc_a_gt_b)
    # ensgs_all_b_gt_a = intersect(ensgs_both_b_gt_a, ensgs_fc_b_gt_a) plot_n = plot_n + 1 plot_2lists_wrap(ensgs_all_a_gt_b, ensgs_all_b_gt_a, comb_pval, plot_count =
    # plot_n, plot_title = 'MACS2 bdgdiff\nMAnorm\n4 FC') print('all done!')
} 
