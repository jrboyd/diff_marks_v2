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
    selection_plot_ready = F,
    detail_ready = F,
    current_mode = modes$chip,
    fresh_start_up = T,
    filter = rownames(my_fe)
    
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
      #       keep = c(x_type(), y_type()) == xy_type_choices[1]
      #       default = default[keep]
    }
    checkboxGroupInput(inputId = 'detail_lines', label = 'Detail Cell Lines', choices = cell_lines, selected = default)
  })
  
  output$detail_marks = renderUI({
    if(debug) print("detail_marks")
    default = character()
    if(v$detail_ready){
      default = union(react_mark_x(), react_mark_y())
      keep = c(x_type(), y_type()) == xy_type_choices[1]
      default = default[keep]
    }
    checkboxGroupInput(inputId = 'detail_marks', label = 'Detail Histone Marks', choices = histone_mods, selected = default)
  })
  output$go_clust = renderUI({
    return(selectInput("go_clust", "GO Cluster input", choices = c(1:input$nclust), selected =  c(1:input$nclust), multiple = T))
  })
  output$filter_clust = renderUI({
    return(selectInput("filter_clust", "Filter to cluster", choices = c(1:input$nclust), selected =  c(1:input$nclust), multiple = T))
  })
  
  #respond to actionButton apply_filter_clust by apply cluster filter
  observeEvent(input$apply_filter_clust, {
    print("apply_filter")
    apply_filter_fun()
  })
  #wraps up actions of apply_filter, protecting it from dependency on too many reactives
  apply_filter_fun = function(){
    if(is.null(v$hmap_res)) return()
    clust_sel = as.numeric(input$filter_clust)
    new_filter = unlist(v$hmap_res$cluster_members[clust_sel])
    v$filter <<- new_filter
  }
  #respond to actionButton release_filter_clust by reseting filter
  observeEvent(input$release_filter_clust, {
    print("release_filter")
    v$filter <<- rownames(my_fe)
  })
  
  
  output$detail_plot_ui = renderUI({
    plotOutput('detail_plot', width = input$detail_width * 100, height = 650)
  })
  
  output$detail_plot = renderPlot({
    if(debug) print("detail_plot")
    if(!v$selection_plot_ready){
      plot0()
      text(.5,.5, 'waiting for selection_plot...')
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
    
    
    hmap_res = plot_details(disp_data, list_up, list_dn, sel, lines2plot, marks2plot, plot_type, smoothing_window, input$cluster_plot_type, as.numeric(input$nclust))
    if(!is.null(hmap_res)){ v$hmap_res = hmap_res}else{v$hmap_res = NULL}
  })
  
  output$selection_plot_ui = renderUI({
    plotOutput('selection_plot', 
               dblclick = "selection_plot_dblclick",
               click = 'selection_plot_click',
               brush = brushOpts(id = 'selection_plot_brush', delay = 600, delayType = 'debounce', resetOnNew = T),
               hover = 'selection_plot_hover', width = 100*input$selection_plot_size, height = 100*input$selection_plot_size)
  })
  
  output$selection_plot =  renderPlot({
    if(!v$data_ready){
      plot0()
      text(.5,.5, 'waiting for data...')
      return()
    }
#     apply_filter()
#     release_filter()
    v$selection_plot_ready <<- T
    disp_data = react_get_displayed_data()
    list_up = react_list_up()
    list_dn = react_list_dn()
    sel = react_get_selected()
    sel = intersect(rownames(disp_data), sel)
    name_a = colnames(disp_data)[1]
    name_b = colnames(disp_data)[2]
    bg_opacity = input$bg_opacity
    fg_opacity = input$fg_opacity
    detect_thresh = input$detect_threshold
    
    plot_volcano(disp_data, list_up, list_dn, sel, name_a, name_b, bg_opacity, fg_opacity, detect_thresh)
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
    kept = rownames(my_fe)
    if(v$current_mode == modes$chip){
      keep = apply(my_fe,1,max) > input$detect_threshold
      kept = names(keep)[keep]
    }else if(v$current_mode == modes$rna){
      keep = apply(my_rna,1,max) > input$detect_threshold
      kept = names(keep)[keep]
    }else{#mode must be mixed
      keep = apply(cbind(my_rna, my_fe),1,max) > input$detect_threshold
      kept = names(keep)[keep]
    }
    kept = intersect(kept, v$filter)
    return(kept)
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
    a = input$x_col_name
    b = input$y_col_name
    key = paste(a, b, sep = '_')
    out = MAnorm_res[[key]]
    if(is.null(out)){
      key = paste(b, a, sep = '_')
      out = MAnorm_res[[key]]
      tmp = out$up
      out$up = out$down
      out$down = tmp
    }
    return(out)
  })
  
  react_DESeq2 = reactive({
    if(debug) print('react_MAnorm')
    out = react_loadDESeq2()
    if(length(out) > 0){
      pval_thresh = input$pval_threshold
      keep = -log10(out$up) > pval_thresh
      out$up = names(out$up)[keep]
      keep = -log10(out$down) > pval_thresh
      out$down = names(out$down)[keep]
    }
    return(out)
    
  })
  
  react_loadDESeq2 = reactive({
    if(debug) print('react_DESeq2')
    a = input$x_col_name
    b = input$y_col_name
    
    key = paste(a, 'vs', b, 'UP', sep = '_')
    out = list()
    out$down = deseq_results[[key]]
    out$up = deseq_results[[sub("UP", "DOWN", key)]]
    if(is.null(out$up)){
      key = paste(b, a, sep = '_')
      out$down = deseq_results[[key]]
      out$up = deseq_results[[sub("UP", "DOWN", key)]]
    }
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
    a = input$x_col_name
    b = input$y_col_name
    key = paste(a, b, sep = '_')
    out = MACS2_bdgdiff_res[[key]]
    if(is.null(out)){
      key = paste(b, a, sep = '_')
      out = MACS2_bdgdiff_res[[key]]
      tmp = out$up
      out$up = out$down
      out$down = tmp
    }
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
    list_fun = list(react_FC())#if mode is mixed, will just use FC
    names(list_fun) = marks_mismatch_message
    if(v$current_mode == modes$chip){
      list_fun = list(
        react_FC(),
        react_MAnorm(),
        react_MACS2())
      names(list_fun) = selection_method_choices#hardcoded from ui.R, room for improvement
    }else if(v$current_mode == modes$rna){
      list_fun = list(
        react_FC(),
        react_DESeq2())
      names(list_fun) = selection_method_choices_rna#hardcoded from ui.R, room for improvement
    }
    if(is.null(sel_methods)){
      if(debug) 
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
    out_list = character()
    if(!is.null(sel_methods)){
      out_list = process_lists(direction, sel_methods)
      out_list = intersect(out_list, react_get_detectable())
    }else if(v$fresh_start_up){
      v$fresh_start_up = F
      out_list = process_lists(direction, selection_method_choices[1])
      out_list = intersect(out_list, react_get_detectable())
    }
    return(out_list)
  })
  
  react_list_dn = reactive({
    if(debug) print('react_list_dn')
    direction = 'down'
    sel_methods = input$available_methods
    if(length(sel_methods) > 0 && sel_methods == marks_mismatch_message){
      sel_methods = selection_method_choices[1]
    }
    out_list = character()
    if(!is.null(sel_methods)){
      out_list = process_lists(direction, sel_methods)
      out_list = intersect(out_list, react_get_detectable())
    }
    return(out_list)
  })
  
  react_x_values = function(){
    if(debug) print('react_x_vals')
    dat = NULL
    if(any(input$x_col_name == colnames(my_fe))){
      dat = my_fe[,input$x_col_name]
    }else if(any(input$x_col_name == colnames(my_rna))){
      dat = my_rna[,input$x_col_name]
    }else{
      stop(paste("input$x_col_name", input$x_col_name, "not in any loaded data source"))
    }
    return(dat)
  }
  
  react_y_values = function(){
    if(debug) print('react_y_vals')
    dat = NULL
    if(any(input$y_col_name == colnames(my_fe))){
      dat = my_fe[,input$y_col_name]
    }else if(any(input$y_col_name == colnames(my_rna))){
      dat = my_rna[,input$y_col_name]
    }else{
      stop(paste("input$y_col_name", input$y_col_name, "not in any loaded data source"))
    }
    return(dat)
  }
  
  #   react_deseq = reactive({
  #     dpair = input$deseq_pair
  #     if(dpair == deseq_groups[1]){#none option selected
  #       return(NULL)
  #     }
  #     dpair = gsub(' ', '_', dpair)
  #     line_a = strsplit(dpair, '_')[[1]][1]
  #     line_b = strsplit(dpair, '_')[[1]][3]
  #     dpair = paste(dpair, c('UP', 'DOWN'), sep = '_')
  #     up_res = deseq_results[[dpair[2]]]
  #     down_res = deseq_results[[dpair[1]]]
  #     
  #     keep = -log10(up_res) > input$pval_threshold
  #     
  #     up_res = names(up_res)[keep]
  #     keep = (my_rna[up_res,line_b] - my_rna[up_res,line_a]) > input$fc_threshold
  #     up_res = up_res[keep]
  #     keep = -log10(down_res) > input$pval_threshold
  #     down_res = names(down_res)[keep]
  #     keep = my_rna[down_res,line_b] - my_rna[down_res,line_a] < -input$fc_threshold
  #     down_res = down_res[keep]
  #     out = list(up = up_res, down = down_res)
  #     
  #     return(out)
  #   })
  
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
  observeEvent(input$selection_plot_dblclick, {
    if(debug) print('dblclick')
    v$click1 <- input$selection_plot_dblclick
    v$brush = NULL
  })
  
  observeEvent(input$selection_plot_click, {
    if(debug) print('click')
  })
  
  
  
  # Handle bush on the plot
  observeEvent(input$selection_plot_brush, {
    if(debug) print('brush')
    v$brush = input$selection_plot_brush
    v$n <- v$n + 1
    vb = v$brush
    v$click1 = NULL
  })
  
  # Handle bush on the plot
  observeEvent(input$selection_plot_hover, {
    v$hover = input$selection_plot_hover
  })
  
  observeEvent(input$reset, {
    # Reset both the range and the first click, if any.
    v$range <- NULL
    v$click1 <- NULL
  })
  
  output$available_methods = renderUI({
    if(debug) print(v$current_mode)
    n_hist = length(unique(histone_mods))
    if(is.null(react_xy_dat())){
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = 'waiting on data...'))
    }
    mark1 = react_mark_x()
    mark2 = react_mark_y()
    if(v$current_mode == modes$chip){
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = selection_method_choices, selected = selection_method_choices[1]))  
    }else if(v$current_mode == modes$rna){
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = selection_method_choices_rna, selected = selection_method_choices_rna[1]))  
    }else{
      return(checkboxGroupInput(inputId = 'available_methods', label = 'Differential Methods', choices = marks_mismatch_message, selected = marks_mismatch_message))  
    }
  })
  outputOptions(output, "available_methods", suspendWhenHidden = FALSE)
  
  x_type = reactive({
    if(input$x_type != input$y_type){
      v$current_mode = modes$mixed
    }else if(input$x_type == xy_type_choices[1]){
      v$current_mode = modes$chip
    }else if(input$x_type == xy_type_choices[2]){
      v$current_mode = modes$rna
    }else{
      stop("unrecognized x_type encountered")
    }
    
    
    # input$detail_type
    if(is.null(input$x_type)){
      return(xy_type_choices[1])
    }else{
      return(input$x_type)
    }
  })
  
  y_type = reactive({
    if(input$x_type != input$y_type){
      v$current_mode = modes$mixed
    }else if(input$x_type == xy_type_choices[1]){
      v$current_mode = modes$chip
    }else if(input$x_type == xy_type_choices[2]){
      v$current_mode = modes$rna
    }else{
      stop("unrecognized x_type encountered")
    }
    # input$detail_type
    if(is.null(input$y_type)){
      return(xy_type_choices[1])
    }else{
      return(input$y_type)
    }
  })
  
  output$x_select = renderUI({
    
    input$detail_type
    choices = tryCatch(expr = {
      if(x_type() == xy_type_choices[[2]]){
        cn = colnames(my_rna)
      }else{
        cn = colnames(my_fe)
      }
      cn
    }, error = function(e){
      return(colnames(my_fe))
    })
    return(selectInput(inputId = 'x_col_name', label = 'Select "From" Value ', choices = choices, selected = choices[1]))
    
  })
  
  output$y_select = renderUI({
    input$detail_type
    choices = tryCatch(expr = {
      if(y_type() == xy_type_choices[[2]]){
        cn = colnames(my_rna)
      }else{
        cn = colnames(my_fe)
      }
      cn
      
    }, error = function(e){
      return(colnames(my_fe))
    })
    return(selectInput(inputId = 'y_col_name', label = 'Select "To" value ', choices = choices, selected = choices[2]))
  })
  
  
  
  output$goTable = renderTable({
    if(!v$selection_plot_ready){
      return(xtable(data.frame('waiting on selection_plot plot...')))
    }
    if(debug) print("goTable")
    
    clust_sel = as.numeric(input$go_clust)
    sel_msig = input$msig_choices
    binom_res = get_goTable(sel_msig, clust_sel, v$hmap_res)
    if(nrow(binom_res) == 0){
      return(xtable(as.data.frame("no significant enrichment found!")))
    }
    out_table = xtable(as.data.frame(binom_res))
    return(out_table)
  }, sanitize.text.function = force)   
  
  
  output$selTable = renderTable({
    if(!v$selection_plot_ready){
      return(xtable(data.frame('waiting on selection_plot plot...')))
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
    fname = paste('selectedTable_',fname, '.xlsx', sep = '')
  })
  
  output$dl_table = downloadHandler(
    filename = dl_tablename,
    content = function(file){
      sel = react_get_selected()
      content_table(file, sel, v$hmap_res)
    }
  )
  
  ###selection_plot plot download
  dl_selection_plot_name = reactive({
    fname = dl_name()
    fname = paste('selection_plot_',fname, '.pdf', sep = '')
  })
  content_selection_plot = function(file){
    pdf(file, width = input$selection_plot_size*100/50, height = input$selection_plot_size*100/50)
    
    disp_data = react_get_displayed_data()
    
    list_up = react_list_up()
    list_dn = react_list_dn()
    sel = react_get_selected()
    sel = intersect(rownames(disp_data), sel)
    name_a = colnames(disp_data)[1]
    name_b = colnames(disp_data)[2]
    bg_opacity = input$bg_opacity
    fg_opacity = input$fg_opacity
    detect_thresh = input$detect_threshold
    
    plot_volcano(disp_data, list_up, list_dn, sel, name_a, name_b, bg_opacity, fg_opacity, detect_thresh)
    dev.off()
  }
  output$dl_selection_plot = downloadHandler(
    filename = dl_selection_plot_name,
    content = content_selection_plot
  )
  
  ##detail plot download
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
    
    
    hmap_res = plot_details(disp_data, list_up, list_dn, sel, lines2plot, marks2plot, plot_type, smoothing_window, input$cluster_plot_type, as.numeric(input$nclust))
    dev.off()
  }
  output$dl_detail = downloadHandler(
    filename = dl_detail_name,
    content = content_detail
  )
  
  ##goTable download
  dl_goTable_name = reactive({
    fname = dl_name()
    fname = paste('goTable_',fname, '.xlsx', sep = '')
  })
  
  output$dl_goTable = downloadHandler(
    filename = dl_goTable_name,
    content = function(file){
      clust_sel = as.numeric(input$go_clust)
      sel_msig = input$msig_choices
      content_goTable(file, sel_msig, clust_sel, v$hmap_res)
    }
  )
  
})