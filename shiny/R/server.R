source("R/packages.R")
source("R/data.R")
source("R/filter.R")
source("R/boxplot.R")

truncate_options = list(
                     # order = list(list(2, 'desc')),
                     columnDefs = list(list(
                     targets = "_all",
                     render = JS(
                       "function(data, type, row, meta) {",
                       "return type === 'display' && data != null && data.length > 50 ?",
                       "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                       "}")
                   )))
truncate_options_sort = list(
                     order = list(list(2, 'desc')),
                     columnDefs = list(list(
                     targets = "_all",
                     render = JS(
                       "function(data, type, row, meta) {",
                       "return type === 'display' && data != null && data.length > 50 ?",
                       "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                       "}")
                   )))

server = shinyServer(function(input, output, session) {
  hide(id = "loading-content", anim = TRUE, animType = "fade")
  hide(id = "loading-content-2", anim = TRUE, animType = "fade")
  shinyjs::show("app-content")
  shinyjs::show("app-content-2")

  stat = filter(gdsc_stat, input)
  
  output$cancer_specific_pathway_title_1 = renderText({ 
    if(length(input$cancer) > 0){
      abbre = cancer_abbreviations[input$cancer, ][1]
      return(paste(abbre, "|", input$drug, sep=" "))
    }
  })

  output$PathwayCancerTest_UpStream = DT::renderDataTable({
      if(length(input$drug) > 0 & length(input$cancer) > 0){
        stat()
      }
    },
    options = truncate_options,
    class = "display",
    selection = list(mode = "single", target = "row", selected = c(1))
  )


  output$PathwayBoxplot_UpStream <- renderPlotly({
    if(length(input$drug) > 0 & length(input$cancer) > 0){
      table = stat()
      pick = input$PathwayCancerTest_UpStream_rows_selected
      if(length(pick) > 0 & nrow(table) > 1){
        row = table[pick, ]
        pathway = as.character(row[1])
        disease = input$cancer
        id = 2
        if(input$cancer == "All diseases"){
          id = 3
          disease = as.character(row[2])
        }

        drug = input$drug 
        if(drug == "All drugs"){
          if(input$cancer != "All diseases"){
            id = 3 
          } else{
            id = 4
            disease = as.character(row[3])
          }
          drug = as.character(row[2])
        }

        gdsc_pas_boxplot(drug=drug, pathway=pathway,
                         disease, tolog = input$cbPlotlogscale)
      }
    }
  })


  output$TCGA_Invididual_Up_Boxplot <- renderPlotly({
    if(length(input$drug) > 0 & length(input$cancer) > 0){
      table = stat()
      pick = input$PathwayCancerTest_UpStream_rows_selected
      if(length(pick) > 0 & nrow(table) > 1){
        row = table[pick, ]
        pathway = as.character(row[1])
        disease = input$cancer
        id = 2
        if(input$cancer == "All diseases"){
          id = 3
          disease = as.character(row[2])
        }

        drug = input$drug 
        if(drug == "All drugs"){
          if(input$cancer != "All diseases"){
            id = 3 
          } else{
            id = 4
            disease = as.character(row[3])
          }
          drug = as.character(row[2])
        }

        tcga_pas_boxplot(drug=drug, pathway=pathway,
                         disease, tolog = input$cbPlotlogscale)
      }
    }
  })


  output$pas_auc_df = DT::renderDataTable({
    top_csps <- top_csps[, c(1, 2, 4, 6, 11, 12)]
    top_csps <- top_csps[order(top_csps$t_stat, decreasing = T), ]
    rownames(top_csps) <- NULL
    for(i in 3:6) {
      top_csps[, i] = round(top_csps[, i], 2)
    }
    top_csps

  }, 
    options = truncate_options,
    class = "display",
    selection = list(mode = "single", target = "row", selected = c(1))
  )

  output$pas_auc_title = renderText({
      "For top 20% CSPs from GDSC, 40 CSPs are significant in BeatAML"
  })


  output$RDR_fig_title = renderText({
    "Discovery and validation of the significant association between pathway activation score (PAS)
      and drug response(AUC) in GDSC (discovery) and BeatAML (validation). The rediscovery rate (RDR)
      is calculated by comparing top 5%, 10%, 20%, 30%, 40% ,50%, 100% CSPs in the discovery set to 
      significant CSPs in the validation set"
  })

  
  output$drug_description = renderText({
    drug = input$drug
    if(length(drug) > 0){
      if(drug == "All drugs"){
        table = stat()
        pick = input$PathwayCancerTest_UpStream_rows_selected
        if(length(pick) > 0 & nrow(table) > 1){
          row = table[pick, ]
          drug = as.character(row[2])
        }
      }

      pick = which(drug == beatAML_drug_description$Drug)
      if(length(pick) == 0){
        return("Updating...")
      } else{
        return(beatAML_drug_description[pick, 2])
      }
    }
  })

  output$drug_source = renderText({
    drug = input$drug
    if(length(drug) > 0){
      if(drug == "All drugs"){
        table = stat()
        pick = input$PathwayCancerTest_UpStream_rows_selected
        if(length(pick) > 0 & nrow(table) > 1){
          row = table[pick, ]
          drug = as.character(row[2])
        }
      }

      pick = which(drug == beatAML_drug_description$Drug)
      if(length(pick) == 0){
        return("Updating...")
      } else{
        return(beatAML_drug_description[pick, 3])
      }
    }
  })

  output$drug_link = renderText({
    drug = input$drug
    if(length(drug) > 0){
      if(drug == "All drugs"){
        table = stat()
        pick = input$PathwayCancerTest_UpStream_rows_selected
        if(length(pick) > 0 & nrow(table) > 1){
          row = table[pick, ]
          drug = as.character(row[2])
        }
      }

      pick = which(drug == beatAML_drug_description$Drug)
      if(length(pick) == 0){
        return("Updating...")
      } else{
        return(beatAML_drug_description[pick, 4])
      }
    }
  })

  output$drug_name = renderText({
    drug = input$drug
    if(length(drug) > 0){
      if(drug == "All drugs"){
        table = stat()
        pick = input$PathwayCancerTest_UpStream_rows_selected
        if(length(pick) > 0 & nrow(table) > 1){
          row = table[pick, ]
          drug = as.character(row[2])
        }
      }
    }
    return(drug)
  })


  output$median_table <- DT::renderDataTable({
    if(length(input$cancer) > 0 & length(input$drug) > 0){
      table = stat()
      pick = input$PathwayCancerTest_UpStream_rows_selected
      if(length(pick) > 0 & nrow(table) > 1){
        pathway = as.character(table[pick, 1])
        drug = input$drug
        if(drug == "All drugs"){
          drug = as.character(table[pick, 2])
        }
        individual_boxplot_table(drug = drug, pathway = pathway,
                                 disease = input$cancer,
                                 to_log=input$cbPlotlogscale)
      }
    }
  },
    options = truncate_options_sort,
    class = "display",
    selection = list(mode = "single", target = "row", selected = c(1))
  )

  RDR_res = plot_RDR(left_tail = TRUE)

  output$left_RDR_plt = renderPlotly({
    plt = ggplotly(RDR_res, source="RDR_plt")
    plt
  })


  output$pas_auc_GDSC_fig = renderPlotly({
    click_event = event_data(event = "plotly_click", source = "RDR_plt")
    pick = input$pas_auc_df_rows_selected
    if(length(pick) > 0){
      return(plot_pas_auc_fig(percent, p, pick, is_gdsc=TRUE))
    }
  })

  output$pas_auc_BeatAML_fig = renderPlotly({
    pick = input$pas_auc_df_rows_selected
    if(length(pick) > 0){
      return(plot_pas_auc_fig(percent, p, pick, is_gdsc=FALSE))
    }
  })
})

