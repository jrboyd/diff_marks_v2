source("scripts/functions_movingAverage.R")

load_ngsprofiles = function(mark_data) {
  #if(exists('ngs_profiles')) return(ngs_profiles)
  ngs_file = "data/ngs_profiles.save"
  if (!file.exists(ngs_file)) {
    print("pre-calculated ngs data not found, calculating profiles...")
    ngs_profiles = list()
    for (n in sub(" ", "_", colnames(mark_data))) {
      
      in_dir = dir("data/ngsplot_data", full.names = T, pattern = n)
      
      fname = paste(in_dir, "/heatmap.RData", sep = "")
      print(paste("loading ngs data from", fname, "..."))
      load(fname)
      tmp = enrichList[[1]]
      ensgs = rownames(tmp)
      ensgs = unlist(lapply(strsplit(ensgs, ":"), function(x) return(x[1])))
      strand = tmp[2:nrow(tmp), 4]
      dat = tmp
      rownames(dat) = ensgs
      dat = dat[intersect(rownames(dat), rownames(mark_data)), ]
      # rownames(dat) = ensg2cut[rownames(dat)]
      ngs_profiles[[n]] = dat
    }
    save(ngs_profiles, file = ngs_file)
    print(paste("done! file is saved as", ngs_file))
  } else {
    print("loading pre-calculated ngs profiles")
    load(ngs_file)
  }
  return(ngs_profiles)
}


applyWindow = function(dat, win = 10) {
  if (win < 2) 
    return(dat)
  out = matrix(0, nrow = nrow(dat), ncol = ncol(dat)/win)
  for (i in 1:ncol(out)) {
    start = (i - 1) * win + 1
    end = i * win
    # out[,i] = apply(dat[,start:end],1,median)
    out[, i] = rowMeans(dat[, start:end])
  }
  return(out)
}

# plotNGS_geneList = function(geneList, ymax = 4){ }
plotNGS_heatmap = function(sel, samples_toPlot){
  toPlot = matrix(0, ncol = 0, nrow = length(sel))
  for(samp in samples_toPlot){
    dat = ngs_profiles[[samp]][sel,, drop = F]
    toPlot = cbind(toPlot, dat)
  }
  MAX = max(toPlot)
  MIN = min(toPlot)
  absMAX = max(abs(c(MAX, MIN)))
  resltn = 100
  start = round((.5 - (abs(MIN) / absMAX) * .5) * resltn)
  end = round((.5 + (abs(MAX) / absMAX) * .5) * resltn)
  
  hmap_res = NULL
  cr = colorRamp(c('blue', 'black', 'yellow'))
  colors = rgb(cr(start:end/resltn)/255)
  
  xlabs = rep('', ncol(toPlot))
  xlabs[1:length(samples_toPlot) * 100 - 50] = samples_toPlot
  if(nrow(toPlot) > 2){
    hmap_res = heatmap.3(toPlot, nsplits = length(samples_toPlot), classCount = 3, col = colors, labCol = xlabs, cexCol = 3, srtCol = 0, adjCol = c(.5,0), key.title = '', key.xlab = 'logFE')
  } else {
    plot0()
    text(.5,.5, 'select at least 3 genes')
  }
  
  return(hmap_res)
}


plotNGS_wBG = function(fg_ENSGcut_list, bg_ENSGcut_list = NA, list_name, sel_name = "selected", invert = F, ymax = NA, linesToPlot = c("MCF10A", "MCF7", "MDA231"), marksToPlot = c('H3K4AC', 'H3K4ME3'), smoothing = 1, plotColors = NA) {
  # plot ngs profile style plots ENSGcut_list : is a character vector of cut ensg ids, cut means version number removed list_name : name or description of input list,
  # used in title invert : if T, everything not in list will be plotted ymax : the upper ylim for plots
  seg = length(marksToPlot)
  if(is.na(plotColors)){
    plotColors = RColorBrewer::brewer.pal(3, 'Dark2')
    names(plotColors) = cell_lines
  }
  layout(matrix(c(1:(seg+2)), ncol = 1, byrow = T), heights = c(1,rep(3/seg, seg),1))
  par(mai = c(0, 1, 0, 0.2))
  fg_lwd = 3
  bg_lwd = 3
  bg_lty = 2
  bg_pch = 21
  
  xs = 0:100
  xs = (200 * xs) - 10000
  
  bg_keep = character()
  if (!is.na(bg_ENSGcut_list)) {
    bg_keep = bg_ENSGcut_list
    if (length(bg_keep) < 2 && is.na(bg_keep)) {
      bg_keep = ID_SET
    }
    if (invert) {
      tmp = rep(T, length(ID_SET))
      names(tmp) = ID_SET
      tmp[bg_keep] = F
      bg_keep = ID_SET[tmp]
    }
    bg_keep = intersect(bg_keep, ID_SET)
  }
  # if(length(fg_keep) < 2 && is.na(fg_keep)){ fg_keep = ID_SET }
  fg_keep = fg_ENSGcut_list
  if (invert) {
    tmp = rep(T, length(ID_SET))
    names(tmp) = ID_SET
    tmp[fg_keep] = F
    fg_keep = ID_SET[tmp]
  }
  
  fg_keep = intersect(fg_keep, ID_SET)
  
  plot(c(0, 1), c(0, 1), type = "n", axes = F, xlab = "", ylab = "")
  # if(!is.na(bg_ENSGcut_list)){
  text(0.7, 0.5, paste(list_name, "\n", length(bg_keep) + length(fg_keep), " genes", sep = ""), cex = 1.5)
  # }
  legend(x = 0, y = 0.7, legend = linesToPlot, fill = plotColors[linesToPlot], bty = "n")
  if (!is.na(bg_ENSGcut_list)) {
    legend(x = 0, y = 0.4, legend = c(paste(length(fg_keep), "genes in", sel_name), paste(length(bg_keep), "other genes in list")), lty = c(3, 1), bty = "n", lwd = 2)
  }
  
  
  if (is.na(ymax)) {
    ymax = 0
    ymin = Inf
    for (l in linesToPlot) {
      for(m in marksToPlot){
        dat = ngs_profiles[[paste(l, m, sep = "_")]]
        ymin = min(ymin, min(colMeans(dat[fg_keep, , drop = F])))
        ymax = max(ymax, max(colMeans(dat[fg_keep, , drop = F])))
      }
      
      
    }
  }
  
#   plot(c(0, 1), type = "n", xlim = c(-1000, 1000), ylim = c(ymin, ymax), ylab = "H3K4ac", lwd = 2, xaxt = "n")
  
  plot.line_style = function(keep, dat, l) {
    toPlot = colMeans(dat[keep, , drop = F])
    if (smoothing > 1) {
      toPlot = movingAverage(toPlot, n = smoothing)
    }
    lines(xs, toPlot, col = plotColors[l], lwd = fg_lwd, lty = 1)
    
  }
  plot.point_style = function(keep, dat, l) {
    toPlot = colMeans(dat[keep, , drop = F])
    if (smoothing > 1) {
      toPlot = movingAverage(toPlot, n = smoothing)
    }
    points(xs, toPlot, col = plotColors[l], pch = 19, cex = 1)
  }
  
  for(m in marksToPlot){
    plot(c(0, 1), type = "n", xlim = c(-10000, 10000), ylim = c(ymin, ymax), ylab = m, lwd = 2)
    for (l in linesToPlot) {
      dat = ngs_profiles[[paste(l, m, sep = "_")]]
      plot.line_style(fg_keep, dat, l)
      # points(xs, colMeans(dat[fg_keep,]), col = l2col[l], pch = 19, cex = 1)
    }
    if (!is.na(bg_ENSGcut_list)) {
      for (l in linesToPlot) {
        dat = ngs_profiles[[paste(l, m, sep = "_")]]
        plot.point_style(bg_keep, dat, l)
        # lines(xs, colMeans(dat[bg_keep,]), col = l2col[l], lwd = fg_lwd, lty = 1)
      }
    }
    lines(c(-10000,10000), c(0,0), lty = 2, col = 'gray', lwd = 2)
    
  }
#   for (l in linesToPlot) {
#     me_d = ngs_profiles[[paste(l, "H3K4ME3", sep = "_")]]
#     plot.line_style(fg_keep, me_d, l)
#     # points(xs, colMeans(me_d[fg_keep,]), col = l2col[l], pch = 19, cex = 1)
#   }
#   if (!is.na(bg_ENSGcut_list)) {
#     for (l in linesToPlot) {
#       me_d = ngs_profiles[[paste(l, "H3K4ME3", sep = "_")]]
#       plot.point_style(bg_keep, me_d, l)
#       # lines(xs, colMeans(me_d[bg_keep,]), col = l2col[l], lwd = fg_lwd, lty = 1)
#     }
#   }
#   plot(c(0, 1), c(0, 1), type = "n", axes = F, xlab = "", ylab = "")
}


plotNGS = function(ENSGcut_list, list_name, invert = F, ymax = 4, linesToPlot = c("MCF10A", "MCF7", "MDA231")) {
  # plot ngs profile style plots ENSGcut_list : is a character vector of cut ensg ids, cut means version number removed list_name : name or description of input list,
  # used in title invert : if T, everything not in list will be plotted ymax : the upper ylim for plots
  seg = length(marksToPlot)
  layout(matrix(c(1:(seg+2)), ncol = 1, byrow = T), heights = c(1,rep(3/seg, seg),1))
  par(mai = c(0, 1, 0, 0.2))
  xs = 0:100
  xs = (20 * xs) - 1000
  keep = ENSGcut_list
  if (length(keep) < 2 && is.na(keep)) {
    keep = ID_SET
  }
  if (invert) {
    tmp = rep(T, length(ID_SET))
    names(tmp) = ID_SET
    tmp[keep] = F
    keep = ID_SET[tmp]
  }
  plot(c(0, 1), c(0, 1), type = "n", axes = F, xlab = "", ylab = "")
  text(0.5, 0.5, paste(list_name, "\n", length(keep), " genes", sep = ""), cex = 1.5)
  legend(x = "left", legend = lines, fill = l2col[lines], bty = "n")
  plot(c(0, 1), type = "n", xlim = c(-1000, 1000), ylim = c(0, ymax), ylab = "H3K4ac", lwd = 2, xaxt = "n")
  for (l in linesToPlot) {
    ac_d = ac_dat[[l]]
    keep = intersect(keep, rownames(ac_d))
    lines(xs, colMeans(ac_d[keep, , drop = F]), col = l2col[l], lwd = 2)
  }
  plot(c(0, 1), type = "n", xlim = c(-1000, 1000), ylim = c(0, ymax), ylab = "H3K4me3", lwd = 2)
  for (l in linesToPlot) {
    me_d = me_dat[[l]]
    keep = intersect(keep, rownames(me_d))
    lines(xs, colMeans(me_d[keep, , drop = F]), col = l2col[l], lwd = 2)
  }
  plot(c(0, 1), c(0, 1), type = "n", axes = F, xlab = "", ylab = "")
}

sym2cut = function(sym) {
  keep = sapply(ensg2sym, function(x) return(any(x == sym)))
  cut = ensg2cut[names(ensg2sym)[keep]]
  return(cut)
} 
