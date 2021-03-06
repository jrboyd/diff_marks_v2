#contains meatier server code
plot_volcano = function(disp_data, list_up, list_dn, sel, name_a, name_b, bg_opacity, fg_opacity, detect_thresh){

  
  scale = rep(1, nrow(disp_data)) #max(disp_data)
  names(scale) = rownames(disp_data)
  colors = rep(rgb(0,0,0,bg_opacity), nrow(disp_data))
  names(colors) = rownames(disp_data)
  if(length(list_up) > 0){
    colors = scale_colors(data = disp_data, scale = scale, list_in = list_up, bg_color = rgb(0,0,0,bg_opacity), list_color = rgb(1,0,0,fg_opacity), colors = colors)
  }
  if(length(list_dn) > 0){
    colors = scale_colors(data = disp_data, scale = scale, list_in = list_dn, bg_color = rgb(0,0,0,bg_opacity), list_color = rgb(0,1,0,fg_opacity), colors = colors)
  }
  if(length(sel) > 0){
    colors = scale_colors(data = disp_data, scale = scale, list_in = sel, bg_color = rgb(0,0,0,bg_opacity), list_color = rgb(0,0,1,fg_opacity), colors = colors)
  }
  max_str = max(nchar(name_a), nchar(name_b))
  note = paste0(format(name_a, width = max_str), ' - ', length(list_dn), '\n', format(name_b, width = max_str), ' - ', length(list_up))
  MIN = min(disp_data)
  MAX = max(disp_data)
  YMAX = max(apply(disp_data, 1, min))
  plot_merge(data = disp_data, list_a = list_up, list_b = list_dn, colors = colors, note = '',
             xlab = paste(name_b, 'log2 -', name_a, 'log2'), ylab = paste('minimum of', name_a, 'and', name_b), xlim = c(-MAX, MAX), ylim = c(MIN, YMAX), cex = .8)
  text(-MAX + 2*MAX*.02,MIN + (YMAX-MIN)*.98, note, adj = c(0,1))
  
  if(detect_thresh > 0){
    lines(c(-MAX, MAX), c(detect_thresh, detect_thresh), col = 'yellow', lwd = 2)
  }
}

plot_details = function(disp_data, list_up, list_dn, sel, lines2plot, marks2plot, plot_type, smoothing_window, cluster_plot_type, nclust){
  hmap_res = NULL
  if(is.null(lines2plot) || is.null(marks2plot)){
    plot0()
    text(.5,.5, 'please\nselect at least 1 cell line and 1 histone mark')
    return()
  }
  if(length(sel) < 1){
    plot0()
    text(.5,.5, 'no points selected')
    return()
  }
  detail_desc = paste(paste(lines2plot, collapse = ', '), ":", paste(marks2plot, collapse = ', '))
  
  to_plot = unlist(lapply(marks2plot, function(x)paste(lines2plot, x, sep = '_')))
  if(debug) print(paste("plotting in detail", paste(to_plot, collapse = ', ')))
  if(plot_type == detail_plot_types[2]){#ngs profiles
    plotNGS_wBG(sel, bg_ENSGcut_list = NA, list_name = detail_desc, sel_name = 'Selected', linesToPlot = lines2plot, marksToPlot = marks2plot, smoothing = smoothing_window)
  }else if(plot_type == detail_plot_types[3]){#ngs heatmaps
    if(length(sel) < 3){
      plot0()
      text(.5,.5, 'selection too small for ngsheatmap!')
      return()
    }
    sel_prof = lapply(ngs_profiles, function(x){
      return(x[sel,])
    })
    #only do side plot if it won't be confusing
    doSidePlot = min(c(length(marks2plot), length(lines2plot))) == 1
    nclust = min(nclust, length(sel)-1)
    nr = 4 + nclust
    nc = 8
    if(doSidePlot) nc = nc + 1
    lmat_custom = matrix(0, ncol = nc, nrow = nr)
    lmat_custom[2,nc-2] = 1
    lmat_custom[nr-1,nc-1] = 2
    lmat_custom[nr,-2:-1+nc] = 3
    lmat_custom[2,nc] = 4
    lmat_custom[2,nc-1] = 5
    rna_toplot = paste(lines2plot, '_rna', sep = '')
    line2i = 1:3
    names(line2i) = cell_lines
    ed_colors = RColorBrewer::brewer.pal(8, 'Set1')[line2i]
    # if(debug) print(cluster_plot_type)
    hmap_res = heatmap.ngsplots(sel_prof, 
                                nclust = nclust, 
                                cex.col = 1.5, 
                                doSidePlot = doSidePlot, 
                                labelWithCounts = T, 
                                extraData = my_rna,#[,rna_toplot, drop = F], 
                                extraData_colors = ed_colors,
                                extraData_plotFunction = exDat_functions[[cluster_plot_type]],
                                lmat_custom = lmat_custom,cex.row = 2, 
                                labels_right = character(),
                                detail_desc, 
                                profiles_to_plot = to_plot, 
                                forPDF = F, globalScale = 1.2, 
                                labels_below = rep(lines2plot, length(marks2plot)), 
                                labels_above = marks2plot,
                                sidePlot_smoothing = smoothing_window)
    
    plot0(); if(doSidePlot) text(.5,.1, 'average profile', adj = c(.5,0))
    plot0();
    if(cluster_plot_type == exDat_choices[3]){
      text(.5,.7, 'each dash is 2-fold')
    }else{
      for(i in 1:length(cell_lines)){
        cl = cell_lines[i]
        xpos = (i-.5)/length(cell_lines)
        
        text(xpos, 1, cl, adj = c(1,.5), srt = 90)
      }
    }
    
    plot0()
    if(cluster_plot_type == exDat_choices[3]){
      legend('bottom', legend = colnames(my_rna)[line2i], fill = ed_colors, horiz = T, bty = 'n')
    }else if(cluster_plot_type == exDat_choices[2]){
      text(.75,.8, "red is high\ngreen is low")
    }else if(cluster_plot_type == exDat_choices[1]){
      text(.75,.8, "red line is median")
    }
    if(doSidePlot) legend('center', legend = to_plot, fill = RColorBrewer::brewer.pal(8, 'Dark2')[1:length(to_plot)], horiz = T, bty = 'n')
    plot0();text(.5,.1, 'cluster size', adj = c(.5,0))
    plot0();text(.5,.3, 'RNA expression', adj = c(.5,0))
    
  }else if(plot_type == detail_plot_types[4]){#heatmap of all cell lines and mods
    if(length(sel) == 1){
      par(mai = c(2,1,1,1))
      plot(1:ncol(disp_data), disp_data[sel,], axes = F, xlab = '', ylab = 'log2 FE')
      box()
      axis(side = 2)
      axis(side = 1, at = 1:ncol(my_fe), labels = colnames(my_fe), las = 2)
      title(ensg_dict[sel,]$gene_name)
    }else{
      if(length(to_plot) < 2){
        plot0()
        text(.5,.5,'please\nselect more lines and marks to plot')
      }else{
        hmap_res = heatmap.3(disp_data[sel,to_plot,drop = F], nsplits = length(marks2plot), classCount = min(6, length(sel)), main = paste(length(sel), 'selected genes'), key.xlab = 'log2 FE', key.title = '')
      }
      
    }
    
  }else{
    plot0()
    text(.5,.5, 'no detail plot type selected')
  }
  return(hmap_res)
  
}

#depending on selection filter type selected, selection is intersected/excluded using up and down lists
filter_selections = function(filter, sel, list_up, list_dn){
  if(length(sel) > 0){
    if(filter == selection_filter_choices[1]){#up or down
      sel = intersect(sel, union(list_up, list_dn))
    }else if(filter == selection_filter_choices[2]){#up
      sel = intersect(sel, list_up)
    }else if(filter == selection_filter_choices[3]){#down
      sel = intersect(sel, list_dn)
    }else if(filter == selection_filter_choices[4]){#unchanged
      sel = setdiff(sel, list_up)
      sel = setdiff(sel, list_dn)
    }
  }
  return(sel)
}

#table of selected genes is fetched and then formatted for output to xlsx file
content_table = function(file, sel, hmap_res){
  out = get_sel_table(sel, hmap_res)
  html_key = 'Promoter in UCSC'
  html = out[,html_key]
  links = sapply(strsplit(html, '"'), function(x)return(x[4]))
  labels = sapply(strsplit(html, '"'), function(x)return(x[5]))
  labels = sub(">", '', labels)
  labels = sub("</a>", '', labels)
  out[,html_key] = paste('=HYPERLINK(', links,',', labels,')', sep = '"')
  cluster_key = 'Cluster'
  if(any(cluster_key == colnames(out))){
    cluster_nums = sub('</td', '', sapply(strsplit(out[,cluster_key], '>'), function(x)return(x[2])))
    cluster_fill = sapply(strsplit(out[,cluster_key], '"'), function(x)return(x[2]))
    out[,cluster_key] = cluster_nums
  }
  my_wb <- createWorkbook()
  my_wb1 <- createSheet(wb=my_wb, sheetName="selected list")
  my_wb2 <- createSheet(wb=my_wb, sheetName="parameters")
  addDataFrame(x=out, sheet=my_wb1, row.names = F, col.names = T)
  nr = length(getRows(my_wb1))
  nc = length(getCells(getRows(my_wb1)[1]))
  rows = getRows(my_wb1)[2:nr]
  html_i = (1:ncol(out))[colnames(out) == html_key]
  html_col = getCells(rows, html_i)
  for(i in 1:length(links)){
    cell = html_col[[i]]
    setCellValue(cell, labels[i])
    addHyperlink(cell, links[i])
  }
  if(any(cluster_key == colnames(out))){
    cluster_CB = CellBlock(my_wb1, 2, html_i+1, nr-1, 1, create = F)
    colors = unique(cluster_fill)
    for(i in 1:length(colors)){
      my_fill = Fill(foregroundColor = colors[i])
      fill_rows = (1:(nr-1))[cluster_fill == colors[i]]
      CB.setFill(cluster_CB, my_fill, fill_rows, 1)
    }
  }
  
  autoSizeColumn(my_wb1, colIndex = 1:nc)
  addDataFrame(x='parameters go here', sheet=my_wb2)
  saveWorkbook(my_wb, file)
  #write.xlsx2(out, file = file, sheetName = '1', row.names = F, col.names = T)#, quote = F, sep =',')
}

sel_msig2gene_set = function(sel_msig){
  set_regex = NULL
  invert_regex = F
  gene_set = sel_msig#most msig_choices match msig keys exactly
  if(debug) print(sel_msig)
  if(sel_msig == 'all') gene_set = NULL
  else if(sel_msig == 'c2-cgp'){#most complicated since biocarta, kegg and reactome must be removed
    gene_set = 'c2'
    set_regex = c("KEGG", "BIOCARTA", "REACTOME")
    invert_regex = T
  }else if(sel_msig == 'c2-biocarta'){
    gene_set = 'c2'
    set_regex = "BIOCARTA"
  }else if(sel_msig == 'c2-kegg'){
    gene_set = 'c2'
    set_regex = "KEGG"
  }else if(sel_msig == 'c2-reactome'){
    gene_set = 'c2'
    set_regex = "REACTOME"
  }else if(sel_msig == 'GO:BP'){
    gene_set = sel_msig
  }
  out = list(set_regex = set_regex, invert_regex = invert_regex, gene_set = gene_set)
  return(out)
}

get_gene_set = function(set_name){
  keep = sapply(msigdb_categories, function(x){
    return(any(x$names == set_name))
  })
  gene_set = character()
  if(any(keep)){#set_name is one of msig
    cat = names(keep)[keep]
    cat = msigdb_categories[[cat]]
    gene_set = cat$sets[cat$names == set_name][[1]]
  }else{#must be from go
    gene_set = gos2sym[[set_name]]
  }
  return(gene_set)
}

get_goTable = function(sel_msig, clust_sel, hmap_res){
  cluster_members = hmap_res$cluster_members
  go_input = unlist(hmap_res$cluster_members[clust_sel])
  if(is.null(go_input)){
    return((as.data.frame("no data has been selected for enrichment testing!")))
  }
  params = sel_msig2gene_set(sel_msig)
  set_regex = params$set_regex
  invert_regex = params$invert_regex
  gene_set = params$gene_set

  if(sel_msig == 'GO:BP'){
    binom_res = binom_go_enrich(test_ensg = go_input)
  }else{
    go_input = ensg_dict[go_input,]$gene_name
    binom_res = binom_msig_enrich(gene_list = go_input, gene_set = gene_set, set_regex = set_regex, invert_regex = invert_regex)
  }
  return(binom_res)
}

content_goTable = function(file, out){
  #out = get_goTable(sel_msig, clust_sel, hmap_res)
  my_wb <- createWorkbook()
  my_wb1 <- createSheet(wb=my_wb, sheetName="enrichment results")
  my_wb2 <- createSheet(wb=my_wb, sheetName="parameters")
  addDataFrame(x=out, sheet=my_wb1, row.names = F, col.names = T)
  nc = ncol(out)
  autoSizeColumn(my_wb1, colIndex = 1:nc)
  addDataFrame(x='parameters go here', sheet=my_wb2)
  saveWorkbook(my_wb, file)
}

#selecte genes from scatterplot assembled in table with related information
get_sel_table = function(sel, hmap_res ){
  if(length(sel) < 1){
    return(matrix('no data selected'))
  }
  if(is.null(hmap_res)){
    sel_as_symbols = ensg_dict[sel,]$gene_name
    sel_as_position = ensg_dict[sel,]$ucsc
    base_url = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jrboyd&hgS_otherUserSessionName=TM_K4_with_peaks'
    if(prostate) base_url = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jrboyd&hgS_otherUserSessionName=prostate_old_FE'
    sel_as_urls = paste0(base_url, '&position=', sel_as_position)
    sel_as_urls = paste0('<a target="_blank" href="', sel_as_urls, '">On UCSC</a>') 
    out_table = cbind(sel, sel_as_symbols, sel_as_position, sel_as_urls)
    colnames(out_table) = c('ENSG ID', 'Gene Symbol', 'Position', 'Promoter in UCSC')
    return(out_table)
  }else{
    res = hmap_res
    colors = res[['colors']]
    asPlotted = rownames(res[['as_plotted']])
    classSizes = res[['class_sizes']]
    ensg2colors = rep(colors[1], length(asPlotted))
    names(ensg2colors) = asPlotted
    for(i in 2:length(colors)){
      start = sum(classSizes[1:(i-1)]) + 1
      end = sum(classSizes[1:i])
      ensg2colors[start:end] = colors[i]
    }
    num = 1:length(unique(ensg2colors))
    names(num) = unique(ensg2colors)
    ensg2num = num[ensg2colors]
    sel = asPlotted
    sel_as_symbols = ensg_dict[sel,]$gene_name
    sel_as_position = ensg_dict[sel,]$ucsc
    base_url = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jrboyd&hgS_otherUserSessionName=TM_K4_with_peaks'
    if(prostate) base_url = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jrboyd&hgS_otherUserSessionName=prostate_old_FE'
    sel_as_urls = paste0(base_url, '&position=', sel_as_position)
    sel_as_urls = paste0('<a target="_blank" href="', sel_as_urls, '">On UCSC</a>') 
    out_table = cbind(sel, sel_as_symbols, sel_as_position, sel_as_urls, paste0('<td bgcolor="', ensg2colors, '">',ensg2num,'</td>'))
    colnames(out_table) = c('ENSG ID', 'Gene Symbol', 'Position', 'Promoter in UCSC', 'Cluster')
    rownames(out_table) = NULL
    return(out_table)
  }
}

#based on brush, points in scatterplot are selected and returned as list of rownames from data table
get_selected = function(x, y, list_up, list_dn, filter, v){
  new_selection = character()
  x_plot = y - x
  y_plot = apply(cbind(x,y), 1, min)
  x = x_plot
  y = y_plot
  if(!is.null(v$brush)){
    keep = (x > v$brush$xmin & x < v$brush$xmax) &
      (y > v$brush$ymin & y < v$brush$ymax)
    new_selection = names(x)[keep]
    new_selection = filter_selections(filter, new_selection, list_up, list_dn)
  }else if(!is.null(v$click1)){
    data_dist = abs(x - v$click1$x) + abs(y - v$click1$y)
    closest = names(sort(data_dist))
    closest = filter_selections(filter, closest, list_up, list_dn)
    new_selection = closest[1]
    
  }else{
    new_selection = filter_selections(filter, new_selection, list_up, list_dn)
  }
}