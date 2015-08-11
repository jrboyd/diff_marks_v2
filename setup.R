if(!exists('loaded')){
  library("shiny")
  library("xtable")
  library("RColorBrewer")
  source("scripts/diffpeaks_vs_manorm_vs_fe.R")
  source("scripts//heatmap.3-split.R")
  
  
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
  
  ngs_profiles = promoter_wide_matched_prof
  
  debug = F
  
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
  
  source("scripts//functions_ngsplot.R")
  if(!exists('ngs_profiles')) ngs_profiles = NA
  
  ID_SET = rownames(my_fe)
  l2col = RColorBrewer::brewer.pal(n = 3, "Dark2")
  names(l2col) = lines[1:3] 
  
  plot0 = function(drawBorder = T){
    plot(c(0,1), c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
    if(drawBorder)box()
  }
  marks_mismatch_message = "FC (marks don't match)"
  
  loaded = T
}
