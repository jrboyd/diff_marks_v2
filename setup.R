if(!exists('loaded')){
  library("shiny")
  library("xtable")
  library("RColorBrewer")
  source("scripts/diffpeaks_vs_manorm_vs_fe.R")
  source("scripts//heatmap.3-split.R")
  source("scripts/heatmap.ngsplots_kmeans_with_sideplot.R")
  source("scripts/heatmap.res_lists.R")
  
  pre_calc_dir = 'data/precalc_results/'
  if(!dir.exists(pre_calc_dir)){
    pre_calc_dir = '/slipstream/home/joeboyd/data/precalc_4shiny'  
  }
  if(!dir.exists(pre_calc_dir)){
   stop('no valid pre calc results found! please update pre_calc_dir') 
  }
  to_load = dir(pre_calc_dir, full.names = T)
  for(tl in to_load){
    print(basename(tl))
    load(tl)
  }
  my_fe = log2(promoter_FE_matched)
  colnames(my_fe) = gsub('_prom', '', colnames(my_fe))
  tmp = c(3,4,1,2,5,6)
  my_fe = my_fe[,c(tmp, tmp + 6, tmp + 12)]
  #colnames(my_fe) = gsub('_', ' ', colnames(my_fe))
  
  my_rna = rna_norm_matched
  colnames(my_rna) = sub('_rna', '', colnames(my_rna))
  
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
  
  debug = T
  xy_type_choices = c("ChIPseq", 'RNAseq')
  display_filter_choices = c("Background", "Up", "Down")
  selection_filter_choices = c("Up or Down", "Up", "Down", "Unchanged", 'No Filter')
  selection_method_choices = c("Fold Change", "MAnorm", "MACS2 bdgdiff")
  detail_plot_types = c("None", "ngsplots - profiles", "ngsplots - heatmap", "FE heatmap")
  
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
    plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F, )
  }
  marks_mismatch_message = "FC (marks don't match)"
  
  loaded = T
}
