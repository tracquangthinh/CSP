gdsc_pas_boxplot <- function(drug, pathway, disease, tolog = FALSE){
  command = paste("gdsc_pas_boxplot_", drug, sep="")
  mat = NULL
  command = paste("mat=", command,sep="")
  command = gsub("-", "_", command)
  command = gsub(" ", "_", command)
  command = gsub("\\(", "", command)
  command = gsub("\\)", "", command)
  command = gsub("/", "", command)
  eval(parse(text=command))
  mat = mat[[pathway]]
  
  if(is.null(mat)){
    # par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("This pathway does not express in this drug."),
       cex = 1.6, col = "black")
  } else{
    mat <- unlist(mat)
    groups <- sapply(names(mat), function(s) unlist(strsplit(s, "\\."))[1])
    names(groups) <- NULL

    df = data.frame("PAS" = mat, "Diseases" = groups)
    col = rep("red", nrow(df))
    pick = which(df$Diseases == disease)
    col[pick] = "blue"
    df$col = col
    
    if(tolog) {
      df$PAS = log2(df$PAS + 1)
    }

    plt_GDSC = ggplot() +
      geom_boxplot(data = df, aes(x=Diseases, y=PAS, fill=col)) +
      # geom_point(data = out_df, aes(x=Diseases, y=PAS), color="black", size=0.5) +
      scale_fill_manual(values=c("#f97750", "#57b894"))+
      # geom_rangeframe() +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, size = 18, vjust = 0.7)) +
      theme(axis.text.y = element_text(size = 16))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.title = element_text(size=16))+
      theme(legend.position="none")+
      labs(title="GDSC")+
      theme(plot.title = element_text(hjust = 0.5, size=16))
    return(ggplotly(plt_GDSC))

  }
}

tcga_pas_boxplot <- function(drug, pathway, disease, tolog = FALSE){
  command = paste("tcga_pas_boxplot_", drug, sep="")
  mat = NULL
  command = paste("mat=", command,sep="")
  command = gsub("-", "_", command)
  command = gsub(" ", "_", command)
  command = gsub("\\(", "", command)
  command = gsub("\\)", "", command)
  command = gsub("/", "", command)
  eval(parse(text=command))
  mat = mat[[pathway]]
  
  if(is.null(mat)){
    # par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste("This pathway does not express in this drug."),
       cex = 1.6, col = "black")
  } else{
    mat <- unlist(mat)
    groups <- sapply(names(mat), function(s) unlist(strsplit(s, "\\."))[1])
    names(groups) <- NULL

    df = data.frame("PAS" = mat, "Diseases" = groups)
    col = rep("red", nrow(df))
    pick = which(df$Diseases == disease)
    col[pick] = "blue"
    df$col = col
    
    if(tolog) {
      df$PAS = log2(df$PAS + 1)
    }

    plt_tcga = ggplot() +
      geom_boxplot(data = df, aes(x=Diseases, y=PAS, fill=col)) +
      scale_fill_manual(values=c("#f97750", "#57b894"))+
      # geom_point(data = out_df, aes(x=Diseases, y=PAS), color="black", size=0.5) +
      # geom_rangeframe() +
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, size = 18, vjust = 0.7)) +
      theme(axis.text.y = element_text(size = 16))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(axis.title = element_text(size=16))+
      theme(legend.position="none")+
      labs(title="TCGA")+
      theme(plot.title = element_text(hjust = 0.5, size=16))
    return(ggplotly(plt_tcga))

  }
}

individual_boxplot_table = function(drug, pathway, disease, to_log = FALSE){
  command = paste("gdsc_pas_boxplot_", drug, sep="")
  gdsc_mat = NULL
  command = paste("gdsc_mat=", command,sep="")
  command = gsub("-", "_", command)
  command = gsub(" ", "_", command)
  command = gsub("\\(", "", command)
  command = gsub("\\)", "", command)
  command = gsub("/", "", command)
  eval(parse(text=command))

  command = paste("tcga_pas_boxplot_", drug, sep="")
  tcga_mat = NULL
  command = paste("tcga_mat=", command,sep="")
  command = gsub("-", "_", command)
  command = gsub(" ", "_", command)
  command = gsub("\\(", "", command)
  command = gsub("\\)", "", command)
  command = gsub("/", "", command)
  eval(parse(text=command))

  gdsc_mat = gdsc_mat[[pathway]]
  tcga_mat = tcga_mat[[pathway]]

  med = sapply(gdsc_mat, function(s) s[3])
  TCGA_med = sapply(tcga_mat, function(s) s[3])
  
  if(to_log){
    med = log2(med + 1)
    TCGA_med = log2(TCGA_med + 1)
  }
  
  med <- round(med, 2)
  TCGA_med <- round(TCGA_med, 2)
  cancers <- sapply(names(med), function(s) unlist(strsplit(s, "\\."))[1])
  names(cancers) <- NULL
  pick <- which(cancers == "LAML")
  cancers[pick] = "AML"
  
  df = data.frame(cancers, med, TCGA_med)
  rownames(df) = NULL
  colnames(df) = c("Disease", "GDSC - Median", "TCGA - Median")
  return(df)
}

plot_RDR = function(left_tail = TRUE){
  if ("cor_gdsc" %in% colnames(cor_pas_auc)) {
    cor_pas_auc$cor_ccle <- cor_pas_auc$cor_gdsc
  }
  pick <- c(
    which(is.na(cor_pas_auc$cor_ccle)),
    which(is.na(cor_pas_auc$cor_beataml))
  )
  cor_pas_auc <- cor_pas_auc[-pick, ]
  t_gdsc <- cor_pas_auc$cor_ccle
  m <- mean(t_gdsc)
  std <- sd(t_gdsc)
  p_left_gdsc <- pnorm(t_gdsc, mean = m, sd = std, lower.tail = TRUE)
  p_right_gdsc <- pnorm(t_gdsc, mean = m, sd = std, lower.tail = FALSE)

  t_beataml <- cor_pas_auc$cor_beataml
  t_beataml <- z_transform(t_beataml)
  m <- mean(t_beataml)
  std <- sd(t_beataml)
  p_left_beataml <- pnorm(t_beataml, mean = m, sd = std, lower.tail = TRUE)
  p_right_beataml <- pnorm(t_beataml, mean = m, sd = std, lower.tail = FALSE)

  temp <- cbind(cor_pas_auc,
    p_left_gdsc = p_left_gdsc,
    p_right_gdsc = p_right_gdsc,
    p_left_beataml = p_left_beataml,
    p_right_beataml = p_right_beataml
  )

  temp1 <- temp
  percents <- c(5, 10, 20, 30, 40, 50, 100) / 100
  pick <- round(percents * nrow(temp1))

  if (left_tail) {
    temp1 <- temp1[order(temp1$cor_ccle), ]
  } else {
    temp1 <- temp1[order(temp1$cor_ccle, decreasing = TRUE), ]
  }
  if (left_tail) {
    x <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_left_beataml <= 0.05) / k
    })
  } else {
    x <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_right_beataml <= 0.05) / k
    })
  }
  x <- round(x, 2)

  if (left_tail) {
    x_0.01 <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_left_beataml <= 0.01) / k
    })
  } else {
    x_0.01 <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_right_beataml <= 0.01) / k
    })
  }
  x_0.01 <- round(x_0.01, 2)

  n_samples <- round(percents * nrow(temp1))

  plt_df <- data.frame(
    percents = paste(percents * 100, "%\n (", n_samples, ")", sep = ""),
    RDR = x, group = "0.05"
  )
  extra_df <- data.frame(
    percents = paste(percents * 100, "%\n (", n_samples, ")", sep = ""),
    RDR = x_0.01, group = "0.01"
  )
  plt_df <- rbind(plt_df, extra_df)

  plt_df$percents <- factor(plt_df$percents, levels = unique(plt_df$percents))

  breaks <- c(seq(0.1, 0.5, by = 0.1), 0.01, 0.05)
  labels <- as.character(breaks)

  plt1 <- ggplot(data = plt_df, aes(x = percents, y = RDR, group = group)) +
    geom_line(aes(color = group), size = 3) +
    geom_point(shape = 24, size = 4, fill = "#153354") +
    labs(colour = "p-value", x = "Number of top CSPs") +
    theme_classic() +
    # theme(legend.position = "none") +
    theme(legend.position = c(0.8, 0.8)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 14)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
    ylim(c(0, 0.3)) +
    # scale_color_manual(values=c("#f97750", "#57b894"))+
    scale_color_manual(values = c("#1f77b4", "#ed665d")) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "#1f77b4", size = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "#ed665d", size = 1) +
    scale_y_continuous(
      limits = c(0, 0.3), breaks = breaks, labels = labels,
      name = "RDR"
    ) +
    theme(legend.key.width = unit(3, "line"))

  plt1

}

plot_pas_auc_fig = function(percent, p, pick, is_gdsc=TRUE){
  df = top_csps 
  # drug = df$drug[pick]
  # pathway = df$pathway[pick]
  # pick = paste(drug, pathway, sep="@")
  if(is_gdsc){
    pas = top_pas_auc[[pick]]$pas_gdsc
    auc = top_pas_auc[[pick]]$auc_gdsc
  }else{
    pas = top_pas_auc[[pick]]$pas_beataml
    auc = top_pas_auc[[pick]]$auc_beataml
  }
  df = data.frame(PAS=nscore(pas), AUC=nscore(auc))
  # m = lm(df$PAS ~ df$AUC)
  plt = ggplot(df, aes(x=PAS, y=AUC))+
    geom_point(size=2, color = "#1f77b4", alpha = 0.8)+
    theme_classic()+
    theme(axis.text.x = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 16))+
    # theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.title = element_text(size=16))+
    # geom_abline(intercept = m$coefficients[1], slope = m$coefficients[2], color="#F8766D", size=2)+
    stat_smooth(method="lm", formula = "y~x", se=FALSE, color="#F8766D", size=2)
  
  if(is_gdsc){
    plt = plt + labs(title = "GDSC")
  } else{
    plt = plt + labs(title = "BeatAML")
  }
  plt
}
