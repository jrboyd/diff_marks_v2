source('setup.R')

shinyServer(function(input, output, session) {
  
  
  output$detail_lines = renderUI({
    default = union(lines[react_index_x()], lines[react_index_y()])
    keep = c(input$x_type, input$y_type) == xy_type_choices[1]
    default = default[keep]
    checkboxGroupInput(inputId = 'detail_lines', label = 'Detail Cell Lines', choices = cell_lines, selected = default)
  })
  
  output$detail_marks = renderUI({
    default = union(mods[react_index_x()], mods[react_index_y()])
    keep = c(input$x_type, input$y_type) == xy_type_choices[1]
    default = default[keep]
    checkboxGroupInput(inputId = 'detail_marks', label = 'Detail Histone Marks', choices = histone_mods, selected = default)
  })
  
  output$detail_plot = renderPlot({
    i_x = react_index_x()
    i_y = react_index_y()
    
    disp_data = my_fe
    
    list_up = react_list_up()
    list_up = intersect(rownames(disp_data), list_up)
    list_dn = react_list_dn()
    list_dn = intersect(rownames(disp_data), list_dn)
    sel = react_selected()
    sel = intersect(rownames(disp_data), sel)
    
    clear_hmap_res = T
    if(is.null(input$detail_lines) || is.null(input$detail_marks)){
      plot0()
      text(.5,.5, 'please\nselect at least 1 cell line and 1 histone mark')
      return()
    }
    if(length(sel) < 1){
      plot0()
      text(.5,.5, 'no points selected')
      return()
    }
    detail_desc = paste(paste(input$detail_lines, collapse = ', '), ":", paste(input$detail_marks, collapse = ', '))
    
    to_plot = unlist(lapply(input$detail_marks, function(x)paste(input$detail_lines, x, sep = '_')))
    if(debug) print('to_plot')
    if(debug) print(to_plot)
    if(input$detail_type == detail_plot_types[2]){#ngs profiles
      plotNGS_wBG(sel, bg_ENSGcut_list = NA, list_name = detail_desc, sel_name = 'Selected', linesToPlot = input$detail_lines, marksToPlot = input$detail_marks, smoothing = input$smoothing_window)
    }else if(input$detail_type == detail_plot_types[3]){#ngs heatmaps
      if(length(sel) < 3){
        plot0()
        text(.5,.5, 'selection too small for ngsheatmap!')
        return()
      }
      sel_prof = lapply(ngs_profiles, function(x){
        return(x[sel,])
      })
      #only do side plot if it won't be confusing
      doSidePlot = min(c(length(input$detail_marks), length(input$detail_lines))) == 1
      nclust = min(6, length(sel)-1)
      nr = 4 + nclust
      nc = 6
      if(doSidePlot) nc = nc + 2
      lmat_custom = matrix(0, ncol = nc, nrow = nr)
      lmat_custom[nr-1,nc-2] = 1
      lmat_custom[nr-1,nc-1] = 2
      lmat_custom[nr,-2:-1+nc] = 3
      lmat_custom[1,nc] = 4
      res = heatmap.ngsplots(sel_prof, nclust = nclust, cex.col = 3.3, doSidePlot = doSidePlot, labelWithCounts = T, extraData = my_rna, lmat_custom = lmat_custom,cex.row = 2.5, labels_right = character(),
                             detail_desc, profiles_to_plot = to_plot, 
                             forPDF = F, globalScale = .6, 
                             labels_below = rep(input$detail_lines, length(input$detail_marks)), 
                             labels_above = input$detail_marks)
      
      plot0();text(.5,.5, 'average profile')
      plot0();text(.5,.5, 'log gene expression')
      plot0();legend('center', legend = c('MCF10A', 'MCF7', 'MDA231'), fill = RColorBrewer::brewer.pal(3, 'Set1'), horiz = T, bty = 'n')
      plot0();text(.5,.5, 'cluster size')
      v$hmap_res = res
      clear_hmap_res = F
    }else if(input$detail_type == detail_plot_types[4]){#heatmap of all cell lines and mods
      if(length(sel) == 1){
        par(mai = c(2,1,1,1))
        plot(1:ncol(disp_data), disp_data[sel,], axes = F, xlab = '', ylab = 'log2 FE')
        box()
        axis(side = 2)
        axis(side = 1, at = 1:ncol(my_fe), labels = colnames(my_fe), las = 2)
        title(ensg_dict[sel,]$gene_name)
      }else{
        print(colnames(disp_data))
        print(to_plot)
        res = heatmap.3(disp_data[sel,to_plot,drop = F], nsplits = length(input$detail_marks), classCount = min(6, length(sel)), main = paste(length(sel), 'selected genes'), key.xlab = 'log2 FE', key.title = '')
        v$hmap_res = res
        clear_hmap_res = F
      }
      
    }else{
      plot0()
      text(.5,.5, 'no detail plot type selected')
    }
    if(clear_hmap_res) v$hmap_res = NULL
  })
  
  output$volcano =  renderPlot({
    i_x = react_index_x()
    i_y = react_index_y()
    
    disp_data = react_displayed()
    
    list_up = react_list_up()
    list_up = intersect(rownames(disp_data), list_up)
    list_dn = react_list_dn()
    list_dn = intersect(rownames(disp_data), list_dn)
    
    sel = react_selected()
    sel = intersect(rownames(disp_data), sel)
    
    #print(lines)
    try({
      #print(i_x)
      name_a = column_choices[i_x]
      #print(name_a)
      name_b = column_choices[i_y]
      
      scale = rep(1, nrow(disp_data)) #max(disp_data)
      names(scale) = rownames(disp_data)
      colors = rep(rgb(0,0,0,input$bg_opacity), nrow(disp_data))
      names(colors) = rownames(disp_data)
      if(length(list_up) > 0){
        colors = scale_colors(data = disp_data, scale = scale, list_in = list_up, bg_color = rgb(0,0,0,input$bg_opacity), list_color = rgb(1,0,0,input$fg_opacity), colors = colors)
      }
      if(length(list_dn) > 0){
        colors = scale_colors(data = disp_data, scale = scale, list_in = list_dn, bg_color = rgb(0,0,0,input$bg_opacity), list_color = rgb(0,1,0,input$fg_opacity), colors = colors)
      }
      
      #print(sel)
      if(length(sel) > 0){
        colors = scale_colors(data = disp_data, scale = scale, list_in = sel, bg_color = rgb(0,0,0,input$bg_opacity), list_color = rgb(0,0,1,input$fg_opacity), colors = colors)
      }
      max_str = max(nchar(name_a), nchar(name_b))
      
      note = paste0(format(name_a, width = max_str), ' - ', length(list_up), '\n', format(name_b, width = max_str), ' - ', length(list_dn))
      MIN = min(my_fe[,c(i_x, i_y)])
      MAX = max(my_fe[,c(i_x, i_y)])
      plot_merge(data = disp_data, list_a = list_up, list_b = list_dn, colors = colors, note = note,
                 xlab = paste(name_a, 'log2 FE'), ylab = paste(name_b, 'log2 FE'), xlim = c(MIN, MAX), ylim = c(MIN, MAX), cex = .8)
      detect_thresh = input$detect_threshold
      if(detect_thresh > 0){
        lines(c(MIN, detect_thresh), c(detect_thresh, detect_thresh), col = 'yellow')
        lines(c(detect_thresh, detect_thresh), c(MIN, detect_thresh),  col = 'yellow')
      }
    }, silent = F)
  })
  
  output$select_gene_list = renderUI({
    all_choices = c('a', 'b')
    
    return(selectInput(width = '50%',
                       inputId = 'selected_gene_list', 
                       label = 'Select group for export', 
                       choices = all_choices,
                       selected = all_choices[1]
    ))    
  })
  
  react_displayed = reactive({
    if(debug) print('react_displayed')
    
    
    
    keep = apply(my_fe,1,max) > input$detect_threshold
    displayed_data = react_xy_dat()[keep,]
    displayed_groups = input$display_filter
    #print(displayed_groups)
    list_up = react_list_up()
    list_dn = react_list_dn()
    if(!is.null(react_deseq())){
      deseq_updown = react_deseq()
      list_up = deseq_updown$up
      list_dn = deseq_updown$down
      print(list_up)
    }
    list_up = intersect(list_up, rownames(displayed_data))
    list_dn = intersect(list_dn, rownames(displayed_data))
    
    for(disp_grp in display_filter_choices){
      if(!any(disp_grp == displayed_groups)){#group not selected for display, remove from displayed_data.
        if(disp_grp == display_filter_choices[1]){#bg
          kept = union(list_up, list_dn)
          displayed_data = displayed_data[kept,]
        }else if(disp_grp == display_filter_choices[2]){#up
          kept = setdiff(rownames(displayed_data), list_up)
          #print(kept)
          displayed_data = displayed_data[kept,]
        }else if(disp_grp == display_filter_choices[3]){#dn
          kept = setdiff(rownames(displayed_data), list_dn)
          displayed_data = displayed_data[kept,]
        }else{
          stop('react_displayed : unrecognized display grp')
        }
      }
    }
    return(displayed_data)
  })
  
  react_FC = reactive({
    if(debug) print('react_FC')
    fc_thresh = input$fc_threshold
    out = list()
    i_x = react_index_x()
    i_y = react_index_y()
    dat = react_xy_dat()
    x = dat[,1]
    y = dat[,2]
    keep = y > (x + fc_thresh)
    out$up = rownames(dat)[keep]
    keep = x > (y + fc_thresh)
    out$down = rownames(dat)[keep]
    return(out)
  })
  
  
  react_loadMAnorm = reactive({
    if(debug) print('react_loadMAnorm')
    if(debug) print(paste(react_index_x(), react_index_y()))
    a = colnames(my_fe)[react_index_x()]
    b = colnames(my_fe)[react_index_y()]
    key = paste(a, b, sep = '_')
    out = MAnorm_res[[key]]
    return(out)
  })
  
  react_MAnorm = reactive({
    if(debug) print('react_MAnorm')
    out = react_loadMAnorm()
    pval_thresh = input$pval_threshold
    keep = out$res[out$up,5] > pval_thresh
    out$up = out$up[keep]
    keep = out$res[out$down,5] > pval_thresh
    out$down = out$down[keep]
    return(out)
  })
  
  react_loadMACS2 = reactive({
    if(debug) print('react_loadMACS2')
    if(debug) print(paste(react_index_x(), react_index_y()))
    a = colnames(my_fe)[react_index_x()]
    b = colnames(my_fe)[react_index_y()]
    out = MACS2_bdgdiff_res[[paste(a, b, sep = '_')]]
    return(out)
  })
  
  react_MACS2 = reactive({
    if(debug) print('react_MACS2')
    out = react_loadMACS2()
    pval_thresh = input$pval_threshold
    keep = out$res[out$up,5] > pval_thresh
    out$up = out$up[keep]
    keep = out$res[out$down,5] > pval_thresh
    out$down = out$down[keep]
    return(out)
  })
  
  process_lists = function(direction, sel_methods){
    new_list = character()
    list_fun = list(
      react_FC(),
      react_MAnorm(),
      react_MACS2())
    names(list_fun) = selection_method_choices#hardcoded from ui.R, room for improvement
    if(is.null(sel_methods)){
      print('no lists selected')
      
    }else{
      for(i in 1:length(sel_methods)){
        list_i = list_fun[[sel_methods[i]]][[direction]]
        if(i == 1){
          new_list = list_i
        }else{
          new_list = intersect(new_list, list_i )
        }
      }
    }
    return(new_list)
  }
  
  react_list_up = reactive({
    if(debug) print('react_list_up')
    direction = 'up'
    sel_methods = input$available_methods
    if(length(sel_methods) > 0 && sel_methods == marks_mismatch_message){
      sel_methods = selection_method_choices[1]
    }
    return(process_lists(direction, sel_methods))
  })
  
  
  
  react_list_dn = reactive({
    if(debug) print('react_list_dn')
    direction = 'down'
    sel_methods = input$available_methods
    if(length(sel_methods) > 0 && sel_methods == marks_mismatch_message){
      sel_methods = selection_method_choices[1]
    }
    return(process_lists(direction, sel_methods))
  })
  
  #   react_x_name = reactive({
  #     if(debug) print('react_x_name')
  #     
  #     if(input$x_type == xy_type_choices[1]){
  #       sel = name2index[input$x_values]  
  #     }else{
  #       print(input$x_values)
  #       sel = rna_name2index[input$x_values]  
  #     }
  #     
  #     if(is.null(sel)) sel = 1
  #     if(length(sel) < 1) sel = 1
  #     print(paste('lkasjdflajsd', sel))
  #     return(input$x_type)
  #   })
  
  react_x_values = reactive({
    if(debug) print('react_x_vals')
    
    if(input$x_type == xy_type_choices[1]){
      dat = my_fe[,input$x_values]
    }else{
      dat = my_rna[,input$x_values]
    }
    #print(dat)
    return(dat)
  })
  
  react_y_values = reactive({
    if(debug) print('react_y_vals')
    
    if(input$y_type == xy_type_choices[1]){
      dat = my_fe[,input$y_values]
    }else{
      dat = my_rna[,input$y_values]
    }
    return(dat)
  })
  
  react_deseq = reactive({
    dpair = input$deseq_pair
    if(dpair == deseq_groups[1]){#none option selected
      return(NULL)
    }
    dpair = gsub(' ', '_', dpair)
    line_a = strsplit(dpair, '_')[[1]][1]
    line_b = strsplit(dpair, '_')[[1]][3]
    dpair = paste(dpair, c('UP', 'DOWN'), sep = '_')
    up_res = deseq_results[[dpair[2]]]
    down_res = deseq_results[[dpair[1]]]

    keep = -log10(up_res) > input$pval_threshold
    
    up_res = names(up_res)[keep]
    keep = (my_rna[up_res,line_b] - my_rna[up_res,line_a]) > input$fc_threshold
    print(range(my_rna[up_res,line_b] - my_rna[up_res,line_a]))
    up_res = up_res[keep]
    
    keep = -log10(down_res) > input$pval_threshold
    down_res = names(down_res)[keep]
    keep = my_rna[down_res,line_b] - my_rna[down_res,line_a] < -input$fc_threshold
    down_res = down_res[keep]
    print(up_res)
    out = list(up = up_res, down = down_res)
    
    return(out)
  })
  
  react_xy_dat = reactive({
    if(debug) print('react_xy_dat')
    
    dat = cbind(react_x_values(), react_y_values())
    colnames(dat) = c(input$x_values, input$y_values)
    #print(dat)
    return(dat)
  })
  
  react_index_x = reactive({
    if(debug) print('react_index_x')
    
    if(input$x_type == xy_type_choices[1]){
      sel = name2index[input$x_values]  
    }else{
      sel = rna_name2index[input$x_values]  
    }
    
    if(is.null(sel)) sel = 1
    if(length(sel) < 1) sel = 1
    return(sel)
  })
  
  react_index_y = reactive({
    if(debug) print('react_index_y')
    sel = name2index[input$y_values]
    if(is.null(sel)) sel = 2
    if(length(sel) < 1) sel = 2
    return(sel)
  })
  
  filter_selections = function(filter, sel, list_up, list_dn){
    filter = input$selection_filter
    if(length(sel) > 0){
      if(filter == selection_filter_choices[1]){#up or down
        sel = intersect(sel, union(list_up, list_dn))
        #print(sel)
      }else if(filter == selection_filter_choices[2]){#up
        sel = intersect(sel, list_up)
        #print(sel)
      }else if(filter == selection_filter_choices[3]){#down
        sel = intersect(sel, list_dn)
      }else if(filter == selection_filter_choices[4]){#unchanged
        sel = setdiff(sel, list_up)
        sel = setdiff(sel, list_dn)
      }
    }
    return(sel)
  }
  
  react_selected = reactive({
    if(debug) print('react_selected')
    new_selection = character()
    list_up = react_list_up()
    list_dn = react_list_dn()
    i_x = react_index_x()
    i_y = react_index_y()
    filter = input$selection_filter
    disp_data = react_displayed()
    x = disp_data[,1]
    y = disp_data[,2]
    if(!is.null(v$brush)){
      keep = (x > v$brush$xmin & x < v$brush$xmax) &
        (y > v$brush$ymin & y < v$brush$ymax)
      new_selection = rownames(disp_data)[keep]
      new_selection = filter_selections(filter, new_selection, list_up, list_dn)
    }else if(!is.null(v$click1)){
      data_dist = abs(x - v$click1$x) + abs(y - v$click1$y)
      closest = names(sort(data_dist))
      closest = filter_selections(filter, closest, list_up, list_dn)
      new_selection = closest[1]
      
    }else{
      new_selection = filter_selections(filter, new_selection, list_up, list_dn)
    }
    
    
    return(new_selection)
  })
  
  v <- reactiveValues(
    click1 = NULL,  # Represents the first mouse click, if any
    n = 0,
    brush = NULL,
    hmap_res = NULL
    #selected = character()
  )
  
  
  
  # Handle clicks on the plot
  observeEvent(input$volcano_dblclick, {
    if(debug) print('dblclick')
    v$click1 <- input$volcano_dblclick
    v$brush = NULL
  })
  
  observeEvent(input$volcano_click, {
    if(debug) print('click')
    #v$selected = character()
    #v$brush = NULL
  })
  
  
  
  # Handle bush on the plot
  observeEvent(input$volcano_brush, {
    if(debug) print('brush')
    v$brush = input$volcano_brush
    v$n <- v$n + 1
    vb = v$brush
    v$click1 = NULL
  })
  
  # Handle bush on the plot
  observeEvent(input$volcano_hover, {
    v$hover = input$volcano_hover
  })
  
  observeEvent(input$reset, {
    # Reset both the range and the first click, if any.
    v$range <- NULL
    v$click1 <- NULL
  })
  
  output$available_methods = renderUI({
    n_hist = length(unique(histone_mods))
    if(all(c(input$x_type, input$y_type) == xy_type_choices[1]) && react_index_x() %% n_hist == react_index_y() %% n_hist){
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = selection_method_choices, selected = selection_method_choices[1]))  
    }else{
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = marks_mismatch_message, selected = marks_mismatch_message))  
    }
    
    
  })
  
  
  
  output$x_select = renderUI({
    if(input$x_type == xy_type_choices[[2]]){
      choices = colnames(my_rna)
    }else{
      choices = colnames(my_fe)
    }
    return(selectInput(inputId = 'x_values', label = 'Select X value ', choices = choices, selected = choices[1]))
    
  })
  
  output$y_select = renderUI({
    choices = colnames(my_fe)
    return(selectInput(inputId = 'y_values', label = 'Select Y value ', choices = choices, selected = choices[7]))
  })
  
  get_sel_table = function(sel, hmap_res ){
    if(length(sel) < 1){
      return(matrix('no data selected'))
    }
    if(is.null(hmap_res)){
      sel_as_symbols = ensg_dict[sel,]$gene_name
      sel_as_position = ensg_dict[sel,]$ucsc
      base_url = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jrboyd&hgS_otherUserSessionName=TM_K4_with_peaks'
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
      sel_as_urls = paste0(base_url, '&position=', sel_as_position)
      sel_as_urls = paste0('<a target="_blank" href="', sel_as_urls, '">On UCSC</a>') 
      out_table = cbind(sel, sel_as_symbols, sel_as_position, sel_as_urls, paste0('<td bgcolor="', ensg2colors, '">',ensg2num,'</td>'))
      colnames(out_table) = c('ENSG ID', 'Gene Symbol', 'Position', 'Promoter in UCSC', 'Cluster')
      rownames(out_table) = NULL
      return(out_table)
    }
  }
  
  output$selTable = renderTable({
    disp_data = react_displayed()
    sel = react_selected()
    sel = intersect(rownames(disp_data), sel)
    
    if(length(sel) < 1){
      return(xtable(as.data.frame('no data selected')))
    }
    out_table = xtable(as.data.frame(get_sel_table(sel, v$hmap_res)))
    return(out_table)
  }, sanitize.text.function = force)    
  
  #download handlers
  dl_name = function(){
   
    return('test')
  }
  
  
  dl_tablename = reactive({
    fname = dl_name()
    fname = paste('table_',fname, '.xlsx', sep = '')
  })
  
  content_table = function(file){
    sel = react_selected()
#     ac = active_columns()
#     out = cbind(ensg, ensg2sym[ensg])
#     colnames(out) = c('ensg', 'gene_symbol')
#     for(i in 1:length(ac)){
#       out = cbind(out, padj[ensg,ac[i], drop = F], fc[ensg,ac[i], drop = F], maxes[ensg,ac[i], drop = F])
#       colnames(out)[(3+3*(i - 1)):(5+3*(i - 1))] = paste(ac[i], c('log10 padj', 'log2 fc', 'log2 max'))
#     }
    out = get_sel_table(sel, v$hmap_res)
    html_key = 'Promoter in UCSC'
    html = out[,html_key]
    links = sapply(strsplit(html, '"'), function(x)return(x[4]))
    labels = sapply(strsplit(html, '"'), function(x)return(x[5]))
    labels = sub(">", '', labels)
    labels = sub("</a>", '', labels)
    out[,html_key] = paste('=HYPERLINK(', links,',', labels,')', sep = '"')
     cluster_key = 'Cluster'
     cluster_nums = sub('</td', '', sapply(strsplit(out[,cluster_key], '>'), function(x)return(x[2])))
     cluster_fill = sapply(strsplit(out[,cluster_key], '"'), function(x)return(x[2]))
     out[,cluster_key] = cluster_nums
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
    cluster_CB = CellBlock(my_wb1, 2, html_i+1, nr-1, 1, create = F)
    colors = unique(cluster_fill)
    for(i in 1:length(colors)){
      my_fill = Fill(foregroundColor = colors[i])
      fill_rows = (1:(nr-1))[cluster_fill == colors[i]]
      CB.setFill(cluster_CB, my_fill, fill_rows, 1)
    }
    
    autoSizeColumn(my_wb1, colIndex = 1:nc)
    addDataFrame(x='parameters go here', sheet=my_wb2)
    saveWorkbook(my_wb, file)
    #write.xlsx2(out, file = file, sheetName = '1', row.names = F, col.names = T)#, quote = F, sep =',')
  }
  
  output$dl_table = downloadHandler(
    filename = dl_tablename,
    content = content_table
  )
  
  ###volcano plot download
  dl_volcname = reactive({
    fname = dl_name()
    fname = paste('volcano_',fname, '.pdf', sep = '')
  })
  content_volc = function(file){
    pdf(file)
    plot0()
    text(.5,.5,'empty plot')
    dev.off()
  }
  output$dl_volcano = downloadHandler(
    filename = dl_volcname,
    content = content_volc
  )
  
})