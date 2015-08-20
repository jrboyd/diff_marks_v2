if(!exists('loaded')){
  library("shiny")
  library("xtable")
  library("RColorBrewer")
  library("xlsx")
  library("gage")
  source("scripts//heatmap.3-split.R")
  source("scripts/heatmap.ngsplots_kmeans_with_sideplot.R")
  source("scripts/heatmap.res_lists.R")
  source("scripts/functions_movingAverage.R")
  
  pre_calc_dir = 'data/precalc_results/'
  prostate = T
  if(prostate) pre_calc_dir = "data/precalc_prostate/"
  if(!dir.exists(pre_calc_dir)){
    pre_calc_dir = '/slipstream/home/joeboyd/data/precalc_results'  
  }
  if(!dir.exists(pre_calc_dir)){
   stop('no valid pre calc results found! please update pre_calc_dir') 
  }
  to_load = dir(pre_calc_dir, full.names = T)
  for(tl in to_load){
    print(basename(tl))
    load(tl)
  }
  if(prostate){
    promoter_FE_matched = matched_promoter_FE
    rna_norm_matched = matched_rna
    enst_ref = matched_enst_ref
    ensg_ref = enst_ref
    rownames(ensg_ref) = ensg_ref$gene_id
    load("ref/ensg_dicts.save")
    promoter_wide_matched_prof = matched_ngs_promoters
    
  }
  my_fe = log2(promoter_FE_matched)
  if(!prostate){
    colnames(my_fe) = gsub('_prom', '', colnames(my_fe))
    tmp = c(3,4,1,2,5,6)
    my_fe = my_fe[,c(tmp, tmp + 6, tmp + 12)]
  }else{
    tmp = c(6,2,4,5,1,3)
    my_fe = my_fe[,c(tmp)]
  }
  
  my_rna = log2(rna_norm_matched + 4)
  
  ngs_profiles = lapply(promoter_wide_matched_prof,  function(x){
    MIN = quantile(x, .1)
    MAX = quantile(x, .98)
    x = ifelse(x < MIN, MIN, x)
    x = ifelse(x > MAX, MAX, x)
    u = mean(x)
    sdev = sd(x)
    x = (x - u) / sdev
    return(x)
  })
  
  debug = F
  xy_type_choices = c("ChIPseq", 'RNAseq')
  display_filter_choices = c("Background", "Up", "Down")
  selection_filter_choices = c("Up or Down", "Up", "Down", "Unchanged", 'No Filter')
  selection_method_choices = c("Fold Change", "MAnorm", "MACS2 bdgdiff")
  detail_plot_types = c("None", "ngsplots - profiles", "ngsplots - heatmap", "FE heatmap")
  deseq_groups = c('none', 'MCF10A vs MCF7', 'MCF10A vs MDA231', 'MCF7 vs MDA231')
  exDat_choices = c("Lines", "Heatmaps", "Barplots")
  as_FC = F
  exDat_linePlot = function(dat, ylim){
    xs = 1:ncol(dat)
    if(as_FC){
      dat = dat - dat[,1]
      ylim = c(-5,5)
      dat = ifelse(dat < -5, -5, dat)
      dat = ifelse(dat > 5, 5, dat)
    }
    plot(0, type = 'n', axes = F, xlim = c(.8,ncol(dat) + .2), ylim = ylim)
    hidden = apply(dat, 1, function(x){
      lines(xs, x, col = 'gray')
    })
    
    lines(xs, apply(dat, 2, median), col = "#f62626", lwd = 3)
    
  }
  exDat_hmap = function(dat, ylim){
    par(mai = rep(.1,4))
    cr = colorRamp(c('green', 'red'))
    colors = rgb(cr(0:100/100)/255)
    nclust = min(4, nrow(dat)-1)
    
    if(as_FC){
      dat = dat - dat[,1]
      ylim = c(-5,5)
    }
    if(nclust > 3){
      kclust = kmeans(dat - dat[,1], centers = nclust)
      o = order(kclust$cluster)
      dat = dat[o,, drop = F]
      kclust$cluster = kclust$cluster[o]
      for(i in 1:nclust){
        keep = kclust$cluster == i
        o = order(rowSums(dat[keep,, drop = F]), decreasing = F)
        dat[keep,] = dat[keep,, drop = F][o,, drop = F]
      }
    }
#     sub_hclust = hclust(dist(dat - dat[,1]))
#     dat = dat[sub_hclust$order,]
    breaks = (-1:100 + .5)/100 * (ylim[2] - ylim[1]) + ylim[1]
    breaks[1] = -1000000
    breaks[length(breaks)] = 10000000
    image(0:3, 0:nrow(dat), t(dat), xlim = c(0, 3), ylim = c(0, nrow(dat)), axes = FALSE, xlab = "", ylab = "", col = colors, breaks = breaks)
  }
  exDat_default = NA
  exDat_functions = c(exDat_linePlot, exDat_hmap, exDat_default)
  names(exDat_functions) = exDat_choices
  
  
  column_choices = colnames(my_fe)
  lines = sapply(strsplit(column_choices, '_'), function(x)return(x[1]))
  mods = sapply(strsplit(column_choices, '_'), function(x)return(x[2]))
  
  cell_lines = unique(lines)
  histone_mods = unique(mods)
  
  name2index = 1:ncol(my_fe)
  names(name2index) = column_choices
  
  rna_name2index = 1:ncol(my_rna)
  names(rna_name2index) = colnames(my_rna)
  
  source("scripts//functions_ngsplot.R")
  if(!exists('ngs_profiles')) ngs_profiles = NA
  
  ID_SET = rownames(my_fe)
  l2col = RColorBrewer::brewer.pal(n = 3, "Dark2")
  names(l2col) = lines[1:3] 
  
  plot0 = function(width = 1, height = 1){
    fudge = 0.037037037037
    plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F)
  }
  marks_mismatch_message = "FC (marks don't match)"
  
  loaded = T
}
