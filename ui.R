source('setup.R')
plot_size = "600px"

shinyUI(
  fluidPage(
    tags$div(title="\nHover over other UI elements for additional help.\n
            Purpose: to allow users to easily access and visualize the results of several statistical approaches
             to determine genes that are differentially marked and/or expressed.\n
             To Start: click and drag within the Selection Plot to select some points (genes), and a Detail Plot for those genes will appear.\n
            The most fundamental controls for plotting are under the Main (starting) tab.\n 
             Multiple methods and parameters for determining differential binding/expression are available in the Diff Methods tab.\n
             The Filters tabs controls which points are displayed and which can be selected in the Selection Plot\n
             Additional control of each plot is available under the Plot Appearance tab.",
             
             titlePanel("Diff ChIP/RNAseq (hover here for help)")),
    sidebarLayout(
      sidebarPanel( width = 2,
                    tabsetPanel(type = 'pills',
                                tabPanel("Main",
                                         tags$div(title="Select the data type (ChIPseq or RNAseq) that will be highest in the DOWN condition, appearing on the LEFT half of the plot.",
                                                  radioButtons(inputId = 'x_type', label = '"From" datatype', choices = xy_type_choices, selected = xy_type_choices[1])),
                                         tags$div(title="Select the data type (ChIPseq or RNAseq) that will be highest in the UP condition, appearing on the RIGHT half of the plot.",
                                                  radioButtons(inputId = 'y_type', label = '"To" datatype', choices = xy_type_choices, selected = xy_type_choices[1])),
                                         tags$div(title="Click the dropdown menu to select the data that will be highest in the DOWN condition, appearing on the LEFT half of the Selection Plot.",
                                                  uiOutput(outputId = 'x_select')),
                                         tags$div(title="Click the dropdown menu to select the data that will be highest in the UP condition, appearing on the RIGHT half of the Selection Plot.",
                                                  uiOutput(outputId = 'y_select')),
                                         tags$div(title="Select the type of plot to make afte you've made a selection\n(click and drag on plot to select points)",
                                                  radioButtons(inputId = 'detail_type', label = 'Detail Plot Type', choices = detail_plot_types, selected = detail_plot_types[3])),
                                         tags$div(title="These CELL LINES are available to add to the detail plot\n(updates automatically as \"From\" and \"To\" change)",
                                                  uiOutput(outputId = 'detail_lines')),
                                         tags$div(title="These HISTONE MARKS are available to add to the detail plot\n(updates automatically as \"From\" and \"To\" change)",
                                                  uiOutput(outputId = 'detail_marks'))
                                         
                                ),
                                tabPanel("Diff Methods",
                                         tags$div(title="These DIFFERENTIAL METHODS are available for the pair of data sets selected.",
                                                  uiOutput(outputId = 'available_methods')),
                                         tags$div(title="The MINIMUM -log10 PVALUE considered statistically significant.",
                                                  sliderInput('pval_threshold', label = '-log10 p-value threshold', min = 0, max = 150, value = 9)),
                                         tags$div(title="The MINIMUM log2 FOLD-CHANGE considered biologically significant",
                                                  sliderInput('fc_threshold', label = 'log2 fold-change threshold', min = 0, max = 6, value = 2, step = .5))
                                ),
                                tabPanel("Filters",
                                         tags$div(title="Which groups of genes should be plotted?  By default, all of them.",
                                                  checkboxGroupInput(inputId = 'display_filter', label = 'Point Display Filtering', choices = display_filter_choices, selected = display_filter_choices)),
                                         tags$div(title="Which groups of genes should be selectable?  By default, only those that change significantly.",
                                                  radioButtons(inputId = 'selection_filter', label = 'Point Selection Filtering', choices = selection_filter_choices, selected = selection_filter_choices[1])),
                                         tags$div(title="A minimum detection threshold prevents too many genes with no real signal being plotted.",
                                                  sliderInput('detect_threshold', label = 'log2 detection threshold', min = 0, max = 6, value = 2, step = .25))
                                         
                                         
                                ),
                                tabPanel("Plot Appearance",
                                         tags$div(title="How much of the interface should the Selection Plot take up.",
                                                  sliderInput('selection_plot_size', label = 'Selection Plot Size', min = 3, max = 8, value = 4.5, step = .5)),
                                         tags$div(title="How opaque (opposite of transparent) should points for genes that DON'T CHANGE be?\nBy default, mostly transparent.",
                                                  sliderInput(inputId = 'bg_opacity', label = 'Background Opacity', min = 0.05, max = 1, step = .05, value = .15)),
                                         tags$div(title="How opaque (opposite of transparent) should points for genes that DO CHANGE be?\nBy default, mostly transparent.",
                                                  sliderInput(inputId = 'fg_opacity', label = 'Up/down Opacity', min = 0.05, max = 1, step = .05, value = .8)),
                                         tags$div(title="How wide should the detail plot be.\nIf multiple cell lines/histone marks are displayed, wider is better.",
                                                  sliderInput('detail_width', label = 'Detail Plot Width', min = 6.5, max = 16, value = 8, step = .5)),
                                         tags$div(title="The type of plot made for each cluster to show RNA-expression.",
                                                  radioButtons('cluster_plot_type', label = 'Cluster Plot Type', choices = exDat_choices, exDat_choices[1])),
                                         tags$div(title="How much smoothing to use if aggregated profiles are being plotted\n(if no aggregate plots, this does nothing).",
                                                  sliderInput('smoothing_window', label = 'Smoothing Window', min = 1, max = 20, value = 5, step = 1))
                                )
                    )
      ),
      
      mainPanel(
        fluidRow(
          column(width = 6,
                 tags$div(title="Select points (genes) here to view in the table and detail plot.  \nClick and drag to select a region of points or double-click to select the single nearest",
                          titlePanel("Selection Plot")),
                 uiOutput("selection_plot_ui"),
                 downloadButton('dl_selection_plot', 'Download Selection Plot')
          ),
          column(width = 6,
                 tags$div(title="This plot is meant to reveal additional relationship within genes selected via the Selection Plot.
                          \nAdditional plot type are under the Main tab.\nFor the ngsplot - heatmap, different types of cluster plots are under Plot Appearance.",
                          titlePanel("Detail Plot")),
                 uiOutput("detail_plot_ui"),
                 selectInput("nclust", "Cluster Count", choices = 3:8, selected = 6),
                 uiOutput("filter_clust"),
                 actionButton("apply_filter_clust", "Apply Cluster Filter"),
                 actionButton("release_filter_clust", "Release Cluster Filter"),
                 downloadButton('dl_detail', 'Download Detail Plot')
                 
          )
        )
      )
    ),
    column(width = 6, 
           tags$div(title="Gene currently selected in the Selection Plot appear here.\nThis table may be downloaded as an excel file.",
                    titlePanel("Selected Genes")),
           downloadButton('dl_table', 'Download Selected Genes'),
           tableOutput('selTable')
    ),
    column(width = 6,
           tags$div(title="Select points (genes) here to view in the table and detail plot.  \nClick and drag to select a region of points or double-click to select the single nearest",
                    titlePanel("Gene Set Enrichment")),
           fluidRow(
             selectInput("msig_choices", label = "MSigDB collections", choices = msig_choices, msig_choices[5]),
             uiOutput("go_clust"),
             downloadButton('dl_goTable', "Download Enrichment Results")),
           tableOutput('goTable')
    )
  )
)


