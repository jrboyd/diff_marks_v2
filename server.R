source('setup.R')
source('server_functions.R')
source('server_volcano_functions.R')

shinyServer(function(input, output, session) {
  v <- reactiveValues(
    click1 = NULL,  # Represents the first mouse click, if any
    n = 0,
    brush = NULL,
    hmap_res = NULL,
    data_ready = F,
    volcano_ready = F,
    detail_ready = F
    #selected = character()
  )
  
  react_mark_x = reactive({
    if(debug) print("react_mark_x")
    if(v$data_ready){
      dat_names = colnames(react_xy_dat())
      str = strsplit(dat_names, '_')[[1]][2]
      return(str)
    }
  })
  
  react_mark_y = reactive({
    if(debug) print("react_mark_y")
    if(v$data_ready){
      dat_names = colnames(react_xy_dat())
      str = strsplit(dat_names, '_')[[2]][2]
      return(str)
    }
  })
  
  react_line_x = reactive({
    if(debug) print("react_line_x")
    if(v$data_ready){
      dat_names = colnames(react_xy_dat())
      str = strsplit(dat_names, '_')[[1]][1]
      return(str)
    }
  })
  
  react_line_y = reactive({
    if(debug) print("react_line_y")
    if(v$data_ready){
      dat_names = colnames(react_xy_dat())
      str = strsplit(dat_names, '_')[[2]][1]
      return(str)
    }
  })
  
  output$detail_lines = renderUI({
    if(debug) print("detail_lines")
    react_xy_dat()
    default = character()
    if(v$detail_ready){
      default = union(react_line_x(), react_line_y())
      keep = c(input$x_type, input$y_type) == xy_type_choices[1]
      default = default[keep]
    }
    checkboxGroupInput(inputId = 'detail_lines', label = 'Detail Cell Lines', choices = cell_lines, selected = default)
  })
  
  output$detail_marks = renderUI({
    if(debug) print("detail_marks")
    default = character()
    if(v$detail_ready){
      default = union(react_mark_x(), react_mark_y())
      keep = c(input$x_type, input$y_type) == xy_type_choices[1]
      default = default[keep]
    }
    checkboxGroupInput(inputId = 'detail_marks', label = 'Detail Histone Marks', choices = histone_mods, selected = default)
  })
  
  
  output$detail_plot_ui = renderUI({
    plotOutput('detail_plot', width = input$detail_width * 100, height = 650)
  })
  
  output$detail_plot = renderPlot({
    if(debug) print("detail_plot")
    if(!v$volcano_ready){
      plot0()
      text(.5,.5, 'waiting for volcano...')
      return()
    }
    v$detail_ready <<- T
    disp_data = my_fe
    list_up = react_list_up()
    list_up = intersect(rownames(disp_data), list_up)
    list_dn = react_list_dn()
    list_dn = intersect(rownames(disp_data), list_dn)
    sel = react_get_selected()
    sel = intersect(rownames(disp_data), sel)
    
    lines2plot = input$detail_lines
    marks2plot = input$detail_marks
    plot_type = input$detail_type
    smoothing_window = input$smoothing_window
    
    
    hmap_res = plot_details(disp_data, list_up, list_dn, sel, lines2plot, marks2plot, plot_type, smoothing_window, input$cluster_plot_type)
    if(!is.null(hmap_res)){ v$hmap_res = hmap_res}else{v$hmap_res = NULL}
  })
  
  output$volcano_ui = renderUI({
    plotOutput('volcano', 
               dblclick = "volcano_dblclick",
               click = 'volcano_click',
               brush = brushOpts(id = 'volcano_brush', delay = 600, delayType = 'debounce', resetOnNew = T),
               hover = 'volcano_hover', width = 100*input$volcano_size, height = 100*input$volcano_size)
  })
  
  output$volcano =  renderPlot({
    if(!v$data_ready){
      plot0()
      text(.5,.5, 'waiting for data...')
      return()
    }
    v$volcano_ready <<- T
    disp_data = react_get_displayed_data()
    list_up = react_list_up()
    list_dn = react_list_dn()
    
    sel = react_get_selected()
    sel = intersect(rownames(disp_data), sel)
    
    try({
      name_a = colnames(disp_data)[1]
      name_b = colnames(disp_data)[2]
#       print('name')
#       print(name_a)
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
      if(length(sel) > 0){
        colors = scale_colors(data = disp_data, scale = scale, list_in = sel, bg_color = rgb(0,0,0,input$bg_opacity), list_color = rgb(0,0,1,input$fg_opacity), colors = colors)
      }
      max_str = max(nchar(name_a), nchar(name_b))
      note = paste0(format(name_a, width = max_str), ' - ', length(list_dn), '\n', format(name_b, width = max_str), ' - ', length(list_up))
      MIN = min(disp_data)
      MAX = max(disp_data)
      YMAX = max(apply(disp_data, 1, min))
      plot_merge(data = disp_data, list_a = list_up, list_b = list_dn, colors = colors, note = note,
                 xlab = paste(name_b, 'log2 -', name_a, 'log2'), ylab = paste('minimum of', name_a, 'and', name_b), xlim = c(-MAX, MAX), ylim = c(MIN, YMAX), cex = .8)
      detect_thresh = input$detect_threshold
      if(detect_thresh > 0){
        lines(c(-MAX, MAX), c(detect_thresh, detect_thresh), col = 'yellow')
        #lines(c(detect_thresh, detect_thresh), c(MIN, detect_thresh),  col = 'yellow')
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
  
  #return rownames of detectable datapoints
  react_get_detectable = reactive({
    keep = apply(my_fe,1,max) > input$detect_threshold
    return(names(keep)[keep])
  })
  
  #displayed data has been selected as xy, passes detectable test, 
  #and is in a group (up down, bg) selected for display
  react_get_displayed_data = reactive({
    if(debug) print('react_get_displayed_data')
    keep = react_get_detectable()
    displayed_data = react_xy_dat()[keep,]
    displayed_groups = input$display_filter
    list_up = react_list_up()
    list_dn = react_list_dn()
    
    for(disp_grp in display_filter_choices){
      if(!any(disp_grp == displayed_groups)){#group not selected for display, remove from displayed_data.
        if(disp_grp == display_filter_choices[1]){#bg
          kept = union(list_up, list_dn)
          displayed_data = displayed_data[kept,]
        }else if(disp_grp == display_filter_choices[2]){#up
          kept = setdiff(rownames(displayed_data), list_up)
          displayed_data = displayed_data[kept,]
        }else if(disp_grp == display_filter_choices[3]){#dn
          kept = setdiff(rownames(displayed_data), list_dn)
          displayed_data = displayed_data[kept,]
        }else{
          stop('react_get_displayed_data : unrecognized display grp')
        }
      }
    }
    return(displayed_data)
  })
  
  react_FC = reactive({
    if(debug) print('react_FC')
    fc_thresh = input$fc_threshold
    out = list()
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
  
  #process_lists intersects selected diff methods for up or down
  process_lists = function(direction, sel_methods){
    new_list = character()
    list_fun = list(
      react_FC(),
      react_MAnorm(),
      react_MACS2())
    names(list_fun) = selection_method_choices#hardcoded from ui.R, room for improvement
    if(is.null(sel_methods)){
      if(debug) print('no lists selected')
      
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
    if(is.null(sel_methods)){
      sel_methods = selection_method_choices[1]
    }
    out_list = process_lists(direction, sel_methods)
    out_list = intersect(out_list, react_get_detectable())
    return(out_list)
  })
  
  react_list_dn = reactive({
    if(debug) print('react_list_dn')
    direction = 'down'
    sel_methods = input$available_methods
    if(length(sel_methods) > 0 && sel_methods == marks_mismatch_message){
      sel_methods = selection_method_choices[1]
    }
    if(is.null(sel_methods)){
      sel_methods = selection_method_choices[1]
    }
    out_list = process_lists(direction, sel_methods)
    out_list = intersect(out_list, react_get_detectable())
    return(out_list)
  })
  
  react_x_values = reactive({
    if(debug) print('react_x_vals')
    
    if(input$x_type == xy_type_choices[1]){
      dat = my_fe[,input$x_col_name]
    }else{
      dat = my_rna[,input$x_col_name]
    }
    return(dat)
  })
  
  react_y_values = reactive({
    if(debug) print('react_y_vals')
    
    if(input$y_type == xy_type_choices[1]){
      dat = my_fe[,input$y_col_name]
    }else{
      dat = my_rna[,input$y_col_name]
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
    up_res = up_res[keep]
    keep = -log10(down_res) > input$pval_threshold
    down_res = names(down_res)[keep]
    keep = my_rna[down_res,line_b] - my_rna[down_res,line_a] < -input$fc_threshold
    down_res = down_res[keep]
    out = list(up = up_res, down = down_res)
    
    return(out)
  })
  
  react_xy_dat = reactive({
    if(debug) print('react_xy_dat')
    if(!v$data_ready){
      if(is.null(input$x_col_name) | is.null(input$y_col_name)){
        return(NULL)
      }else{
        v$data_ready <<- T
      }
    }
    dat = cbind(react_x_values(), react_y_values())
    colnames(dat) = c(input$x_col_name, input$y_col_name)
    return(dat)
  })
  
  react_index_x = reactive({
    if(debug) print('react_index_x')
    if(input$x_type == xy_type_choices[1]){
      sel = name2index[input$x_col_name]  
    }else{
      sel = rna_name2index[input$x_col_name]  
    }
    if(is.null(sel)) sel = 1
    if(length(sel) < 1) sel = 1
    return(sel)
  })
  
  react_index_y = reactive({
    if(debug) print('react_index_y')
    sel = name2index[input$y_col_name]
    if(is.null(sel)) sel = 2
    if(length(sel) < 1) sel = 2
    return(sel)
  })
  
  
  
  
  
  react_get_selected = reactive({
    if(debug) print('react_get_selected')
    list_up = react_list_up()
    list_dn = react_list_dn()
    filter = input$selection_filter
    disp_data = react_get_displayed_data()
    x = disp_data[,1]
    y = disp_data[,2]
    new_selection = character()
    if(!is.null(v)){
      new_selection = get_selected(x, y, list_up, list_dn, filter, v) 
    }
    return(new_selection)
  })
  
  # Handle clicks on the plot
  observeEvent(input$volcano_dblclick, {
    if(debug) print('dblclick')
    v$click1 <- input$volcano_dblclick
    v$brush = NULL
  })
  
  observeEvent(input$volcano_click, {
    if(debug) print('click')
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
    if(is.null(react_xy_dat())){
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = 'waiting on data...'))
    }
    mark1 = react_mark_x()
    mark2 = react_mark_y()
    if(all(c(input$x_type, input$y_type) == xy_type_choices[1]) && mark1 == mark2){
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
    return(selectInput(inputId = 'x_col_name', label = 'Select "From" Value ', choices = choices, selected = choices[1]))
    
  })
  
  output$y_select = renderUI({
    choices = colnames(my_fe)
    return(selectInput(inputId = 'y_col_name', label = 'Select "To" value ', choices = choices, selected = choices[7]))
  })
  
  output$goTable = renderTable({
    if(!v$volcano_ready){
      return(xtable(data.frame('waiting on volcano plot...')))
    }
    if(debug) print("goTable")
    sel = react_get_selected()
    writeClipboard((sel))
    if(length(sel) < 1){
      return(xtable(as.data.frame('no data selected')))
    }
    out_table = xtable(as.data.frame(get_sel_table(sel, v$hmap_res)))
    return(out_table)
  }, sanitize.text.function = force)   
  
  
  output$selTable = renderTable({
    if(!v$volcano_ready){
      return(xtable(data.frame('waiting on volcano plot...')))
    }
    if(debug) print("selTable")
    sel = react_get_selected()
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
  
  output$dl_table = downloadHandler(
    filename = dl_tablename,
    content = function(file){
      sel = react_get_selected()
      content_table(file, sel, v$hmap_res)
    }
  )
  
  ###volcano plot download
  dl_detail_name = reactive({
    fname = dl_name()
    fname = paste('detail_',fname, '.pdf', sep = '')
  })
  content_detail = function(file){
    pdf(file, width = input$detail_width*100/50, height = 650/50)
    disp_data = my_fe
    list_up = react_list_up()
    list_up = intersect(rownames(disp_data), list_up)
    list_dn = react_list_dn()
    list_dn = intersect(rownames(disp_data), list_dn)
    sel = react_get_selected()
    sel = intersect(rownames(disp_data), sel)
    
    lines2plot = input$detail_lines
    marks2plot = input$detail_marks
    plot_type = input$detail_type
    smoothing_window = input$smoothing_window
    
    
    hmap_res = plot_details(disp_data, list_up, list_dn, sel, lines2plot, marks2plot, plot_type, smoothing_window, input$cluster_plot_type)
    dev.off()
  }
  output$dl_detail = downloadHandler(
    filename = dl_detail_name,
    content = content_detail
  )
  
})