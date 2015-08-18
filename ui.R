source('setup.R')
plot_size = "600px"

shinyUI(
  fluidPage(
    titlePanel("Diff ChIP/RNAseq"),
    sidebarLayout(
      sidebarPanel( width = 2,
                    tabsetPanel(
                      tabPanel("Main",
                               radioButtons(inputId = 'x_type', label = 'x-axis datatype', choices = xy_type_choices, selected = xy_type_choices[1]),
                               radioButtons(inputId = 'y_type', label = 'y-axis datatype', choices = xy_type_choices, selected = xy_type_choices[1]),
                               uiOutput(outputId = 'x_select'),
                               uiOutput(outputId = 'y_select'),
                               
                               radioButtons(inputId = 'detail_type', label = 'Detail Plot Type', choices = detail_plot_types, selected = detail_plot_types[4]),
                               uiOutput(outputId = 'detail_lines'),
                               uiOutput(outputId = 'detail_marks')
                               
                      ),
                      tabPanel("Diff Methods",
                               uiOutput(outputId = 'available_methods'),
                               sliderInput('pval_threshold', label = '-log10 p-value threshold', min = 0, max = 150, value = 9),
                               sliderInput('fc_threshold', label = 'log2 fold-change threshold', min = 0, max = 6, value = 2, step = .5),
                               radioButtons(inputId = 'deseq_pair', label = 'Deseq2 Results Filtering', choices = deseq_groups, selected = deseq_groups[1])  
                      ),
                      tabPanel("Filters",
                               checkboxGroupInput(inputId = 'display_filter', label = 'Point Display Filtering', choices = display_filter_choices, selected = display_filter_choices),
                               radioButtons(inputId = 'selection_filter', label = 'Point Selection Filtering', choices = selection_filter_choices, selected = selection_filter_choices[1]),
                               sliderInput('detect_threshold', label = 'log2 detection threshold', min = 0, max = 6, value = 2, step = .25)
                               
                      ),
                      tabPanel("Plot Appearance",
                               sliderInput('volcano_size', label = 'Volcano Plot Size', min = 3, max = 8, value = 4.5, step = .5),
                               sliderInput(inputId = 'bg_opacity', label = 'Background Opacity', min = 0.05, max = 1, step = .05, value = .15),
                               sliderInput(inputId = 'fg_opacity', label = 'Up/down Opacity', min = 0.05, max = 1, step = .05, value = .8),
                               sliderInput('detail_width', label = 'Detail Plot Width', min = 6.5, max = 16, value = 8, step = .5),
                               radioButtons('cluster_plot_type', label = 'Cluster Plot Type', choices = exDat_choices, exDat_choices[1]),
                               sliderInput('smoothing_window', label = 'Smoothing Window', min = 1, max = 20, value = 5, step = 1)
                      )
                    )
      ),
      
      mainPanel(
        fluidRow(
          column(width = 6,
                 uiOutput("volcano_ui")
          ),
          column(width = 6,
            uiOutput("detail_plot_ui"),
            downloadButton('dl_detail', 'Download Detail Plot'),
            downloadButton('dl_table', 'Download Selected Table')
          )
        )
      )
      #     
      #     
      #     fluidRow(
      #       downloadButton('dl_table', 'Download Selected Table'),
      #       downloadButton('dl_volcano', 'Download Volcano Plot')
      #     ),
      #     tableOutput('selTable'),
      #     tableOutput('goTable')
    ),
    tableOutput('selTable')
  )
)


