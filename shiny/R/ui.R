# load("./data/disease_with_30groups.RData")
# load("./data/drug_names_all.RData")
# load("./data/pathway_names.RData")
#
# source("./ui/individual_drug/cancer_specific_pathway_ui.R", local=environment())
# source("./ui/individual_drug/aml_drug_response_ui.R", local=environment())
# source("ui/across_drug/across_drug_boxplot.R", local=environment())
source("R/ui_pas.R")
source("R/ui_drug_response.R")
appCSS <- "
  #loading-content, #loading-content-2{
    position: absolute;
    background: #FFFFFF;
    opacity: 0.9;
    z-index: 100;
    left: 0;
    right: 0;
    height: 100%;
    text-align: center;
    color: #000000;
    font-size: 16px;
  }
"

intro_ui = tabPanel("About",
      column(6, offset = 3,
      includeHTML("introduction.html"),
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "ui.css")
      )),
      useShinyjs(),
      inlineCSS(appCSS)
)

ui = navbarPage("CSP",
                intro_ui,
                pas_ui,
                drug_response_ui,
                # across_drug_ui,
                theme = shinytheme("cerulean")
  )
