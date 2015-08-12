source('setup.R')

shinyUI(fluidPage(
  headerPanel('Differential ChIPSeq methods comparison'),
  fluidRow(
    column(width = 6, checkboxGroupInput(inputId = 'display_filter', label = 'Point Display Filtering', choices = display_filter_choices, selected = display_filter_choices)),
    column(width = 6, radioButtons(inputId = 'selection_filter', label = 'Point Selection Filtering', choices = selection_filter_choices, selected = selection_filter_choices[1]))
    
  ),
  fluidRow(
    column(width = 6, uiOutput(outputId = 'available_methods')),
    column(width = 6, radioButtons(inputId = 'deseq_pair', label = 'Deseq2 Results Filtering', choices = deseq_groups, selected = deseq_groups[1]))
  ),
  fluidRow(
    column(width = 6, sliderInput(inputId = 'bg_opacity', label = 'Background Opacity', min = 0.05, max = 1, step = .05, value = .15)),
    column(width = 6, sliderInput(inputId = 'fg_opacity', label = 'Up/down Opacity', min = 0.05, max = 1, step = .05, value = .8))
  ),
  fluidRow(
    column(width = 6, radioButtons(inputId = 'x_type', label = 'x-axis datatype', choices = xy_type_choices, selected = xy_type_choices[1])),
    column(width = 6, radioButtons(inputId = 'y_type', label = 'y-axis datatype', choices = xy_type_choices, selected = xy_type_choices[1]))
  ),
  fluidRow(
    column(width = 6, uiOutput(outputId = 'x_select')),
    column(width = 6, uiOutput(outputId = 'y_select'))
  ),
  fluidRow(
    column(width = 4, radioButtons(inputId = 'detail_type', label = 'Detail Plot Type', choices = detail_plot_types, selected = detail_plot_types[4])),
    column(width = 4, uiOutput(outputId = 'detail_lines')), #checkboxGroupInput(inputId = 'detail_lines', label = 'Detail Cell Lines', choices = cell_lines, selected = detail_plot_types[4])),
    column(width = 4, uiOutput(outputId = 'detail_marks'))#checkboxGroupInput(inputId = 'detail_marks', label = 'Detail Histone Marks', choices = histone_mods, selected = detail_plot_types[4]))
    
  ),
  #   fluidRow(
  #     column(width = 3,radioButtons(inputId = 'updown', label = 'Fold change direction', choices = c('up', 'down', 'either', 'no change'), selected = 'either'))
  #   ),
  fluidRow(
    column(width = 3, sliderInput('detect_threshold', label = 'log2 detection threshold', min = 0, max = 6, value = 2, step = .25)),
    column(width = 3, sliderInput('pval_threshold', label = '-log10 p-value threshold', min = 0, max = 150, value = 9)),
    column(width = 3, sliderInput('fc_threshold', label = 'log2 fold-change threshold', min = 0, max = 6, value = 2, step = .5)),
    column(width = 3, sliderInput('smoothing_window', label = 'Smoothing Window', min = 1, max = 20, value = 5, step = 1))
    #column(width = 4, sliderInput('maxes_threshold', label = 'log2 max threshold', min = 0, max = 16, value = 2))
  ),
  fluidRow(
    column(plotOutput('volcano',
                      dblclick = "volcano_dblclick",
                      click = 'volcano_click',
                      brush = brushOpts(id = 'volcano_brush', delay = 600, delayType = 'debounce', resetOnNew = T),
                      hover = 'volcano_hover'),
           width = 6),
    column(plotOutput('detail_plot'), width = 6)
    
  ),
  fluidRow(
    downloadButton('dl_table', 'Download Selected Table'),
    downloadButton('dl_volcano', 'Download Volcano Plot')
  ),
  tableOutput('selTable')
  #,
  #uiOutput('select_gene_list')
)

)