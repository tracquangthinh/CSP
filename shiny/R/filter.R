z_transform <- function(cors) {
  sapply(cors, function(c) {
    0.5 * (log(1 + c) - log(1 - c))
  })
}

nscore <- function(x) {
  n <- length(x)
  rx <- rank(x)
  ns <- qnorm(rx / (n + 1))
  return(ns)
}

filter <- function(stat, input) {
  res <- reactive({  
    cancer_input <- as.character(input$cancer)
    drug_input <- as.character(input$drug)
    t_threshold <- as.numeric(input$ttestFDRtext)
    chi_threshold <- as.numeric(input$X2testFDRtext)
    fc_threshold <- as.numeric(input$fcThres)

    cancer <- cancer_input
    drug <- drug_input

    if(cancer == "All diseases")
      cancer <- cancers

    if(drug == "All drugs")
      drug <- drugs

    stat <- stat[stat$pt <= t_threshold, ]
    n_sample <- round(chi_threshold * nrow(stat))
    stat <- stat[order(stat$chisq_stat), ]
    stat <- stat[1:n_sample, ]
    stat <- stat[order(stat$t_stat, decreasing = TRUE), ]
    stat$t_stat <- round(stat$t_stat, 2)
    stat$chisq_stat <- round(stat$chisq_stat, 2)
    stat$pt <- round(stat$pt, 2)
    stat$fold_change <- round(stat$fold_change, 2)
    colnames(stat) <- c("Drug", "Pathway", "Disease",
                        "t-statistics", "p-value of t", 
                        "chisq-statistics","Fold change")

    table <- tryCatch({
      stat <- stat[stat$Disease %in% cancer, ]
      stat <- stat[stat$Drug %in% drug, ]

      if(cancer_input  == "All diseases") {
        if(drug_input == "All drugs") {
          stat <- stat[, c("Pathway", "Drug", "Disease",
                           "t-statistics", "p-value of t", 
                           "chisq-statistics","Fold change")]
        } else {
          stat <- stat[, c("Pathway", "Disease",
                           "t-statistics", "p-value of t", 
                           "chisq-statistics","Fold change")]
        }
      } else {
        if(drug_input == "All drugs") {
          stat <- stat[, c("Pathway", "Drug",
                           "t-statistics", "p-value of t", 
                           "chisq-statistics","Fold change")]
        } else {
          stat <- stat[, c("Pathway",
                           "t-statistics", "p-value of t", 
                           "chisq-statistics","Fold change")]
        }
      }
      rownames(stat) <- NULL
      stat
    }, error = function(cond){
         out_table <- matrix("No specific pathway found!", 1, 1)
         rownames(out_table) <- NULL  
         colnames(out_table) <- "Result"
         return(out_table)
    })

    table

  })
}
