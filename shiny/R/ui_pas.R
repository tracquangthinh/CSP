pas_ui = tabPanel("CSP analysis", 
  div(
    id = "loading-content",
    h4(style = "position: relative; top: 50%; transform: translateY(-50%);","Due to the limitation of the bandwidth, sometimes we need to reload large datasets."),
    h4(style = "position: relative; top: 50%; transform: translateY(-50%);","Please do not refresh and wait for few seconds. Sorry for this inconvinience!")
  ),
  
  hidden(
    div(id = "app-content",
  sidebarPanel(width = 3,
    fluidRow(
      # Select cancer disease to plot
      selectizeInput(inputId = "cancer", label = "Select a disease:",
                  choices = c(cancers, "All diseases"),
                  multiple=TRUE, options = list(maxItems = 1),
                  selected = "LAML")
    ),

    fluidRow(     
      # Select drug to plot
      selectizeInput(inputId = "drug", label = "Select a drug:",
                  choices = c("All drugs", drugs), 
                  selected = "Quizartinib",
                  multiple=TRUE, options = list(maxItems = 1)
                  )
      
    ),
    
    fluidRow(
      h5(tags$strong("FDR condition:")), 
      column(width = 12, tags$div(id = "inline", class="t-test-input", numericInput(inputId = "ttestFDRtext", label = withMathJax(
        sprintf('\\(FDR\\ t\\) \\(\\leq\\)')),  value = 0.01, min = 0, max = 0.05, step = 0.01))),
      column(width = 12, tags$div(id = "inline", class="chi-square-test-input", numericInput(inputId = "X2testFDRtext", label = withMathJax(
        sprintf('\\(Top\\ \\chi^2\\):')), value = 0.25, min = 0, max = 1, step = 0.01))),
      column(width = 12, tags$div(id = "inline", class="fold-change-input", numericInput(inputId = "fcThres", label = withMathJax(
        sprintf('\\(Fold-change\\) \\(\\geq\\)')),  value = 1, min = 0, max = 500, step = 0.1)))
    ),
    
    fluidRow(
      h5(tags$strong("Options:")),
      # checkboxInput("cbPlotOutliers", label = "Plot outliers", value = FALSE),
      checkboxInput("cbPlotlogscale", label = "Use log scale", value = FALSE)
    ),
    hr(),
    fluidRow(
      h5(tags$strong("Drug Name:"), textOutput("drug_name")),
      h5(tags$strong("Description:")),
      textOutput("drug_description"),
      h5(tags$strong("Source:"), textOutput("drug_source")),
      h5(tags$strong("Link:"), textOutput("drug_link"))
    )
  ),
  
    mainPanel(width = 9,
      fluidRow(
        h3(textOutput("cancer_specific_pathway_title_1"), align="center")
      ),
      fluidRow(
        DT::dataTableOutput("PathwayCancerTest_UpStream")
      ),
      tags$hr(),
      plotlyOutput("PathwayBoxplot_UpStream"),
      plotlyOutput("TCGA_Invididual_Up_Boxplot"),
      DT::dataTableOutput("median_table")
    ) 
    ))
)
