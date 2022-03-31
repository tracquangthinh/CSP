drug_response_ui= tabPanel("Analysis in AML", 
  div(
      id = "loading-content-2",
      h4(style = "position: relative; top: 50%; transform: translateY(-50%);","Due to the limitation of the bandwidth, sometimes we need to reload large datasets."),
      h4(style = "position: relative; top: 50%; transform: translateY(-50%);","Please do not refresh and wait for few seconds. Sorry for this inconvinience!")
  ),
  fluidRow(
    fluidRow(column(width = 10, offset = 1,
                h4(textOutput("RDR_fig_title"), style = "word-wrap: break-word;"))
             ),
    fluidRow(column(width= 6, plotlyOutput("left_RDR_plt"), offset = 3))
  ),
  hr(),
  fluidRow(
    # fluidRow(verbatimTextOutput("click_info")),
    fluidRow(h4(textOutput("pas_auc_title"), align="center")),
    fluidRow(column(width = 10, fluidRow(DT::dataTableOutput("pas_auc_df")), offset = 1)),
    fluidRow(
      column(width=5, plotlyOutput("pas_auc_GDSC_fig"), offset = 1),
      column(width=5, plotlyOutput("pas_auc_BeatAML_fig"))
    )
      
  )
)
                  
