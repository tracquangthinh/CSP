z_transform <- function(cors) {
  sapply(cors, function(c) {
    0.5 * (log(1 + c) - log(1 - c))
  })
}

# figure 2B-C
draw_rdr_figure <- function(cor_pas_auc, left_tail = TRUE) {
  pick <- c(which(is.na(cor_pas_auc$cor_gdsc)),
            which(is.na(cor_pas_auc$cor_beataml)))
  if (length(pick) > 0) {
    cor_pas_auc <- cor_pas_auc[-pick, ]
  }
  t_gdsc <- cor_pas_auc$cor_gdsc
  m <- mean(t_gdsc)
  std <- sd(t_gdsc)
  p_left_gdsc <-
    pnorm(t_gdsc,
          mean = m,
          sd = std,
          lower.tail = TRUE)
  p_right_gdsc <-
    pnorm(t_gdsc,
          mean = m,
          sd = std,
          lower.tail = FALSE)
  
  t_beataml <- cor_pas_auc$cor_beataml
  t_beataml <- z_transform(t_beataml)
  m <- mean(t_beataml)
  std <- sd(t_beataml)
  p_left_beataml <-
    pnorm(t_beataml,
          mean = m,
          sd = std,
          lower.tail = TRUE)
  p_right_beataml <-
    pnorm(t_beataml,
          mean = m,
          sd = std,
          lower.tail = FALSE)
  
  temp <- cbind(
    cor_pas_auc,
    p_left_gdsc = p_left_gdsc,
    p_right_gdsc = p_right_gdsc,
    p_left_beataml = p_left_beataml,
    p_right_beataml = p_right_beataml
  )
  
  temp1 <- temp
  percents <- c(5, 10, 20, 30, 40, 50, 100) / 100
  pick <- round(percents * nrow(temp1))
  
  if (left_tail) {
    temp1 <- temp1[order(temp1$cor_gdsc), ]
  } else {
    temp1 <- temp1[order(temp1$cor_gdsc, decreasing = TRUE), ]
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
  x <- round(x, 3)
  
  if (left_tail) {
    x_0.01 <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_left_beataml <= 0.01) / k
    })
  } else {
    x_0.01 <- sapply(pick, function(k) {
      sum(temp1[1:k, ]$p_right_beataml <= 0.01) / k
    })
  }
  x_0.01 <- round(x_0.01, 3)
  
  n_samples <- round(percents * nrow(temp1))
  
  plt_df <- data.frame(
    percents = paste(percents * 100, "%\n (", n_samples, ")", sep = ""),
    RDR = x,
    group = "0.05"
  )
  extra_df <- data.frame(
    percents = paste(percents * 100, "%\n (", n_samples, ")", sep = ""),
    RDR = x_0.01,
    group = "0.01"
  )
  plt_df <- rbind(plt_df, extra_df)
  
  plt_df$percents <-
    factor(plt_df$percents, levels = unique(plt_df$percents))
  
  breaks <- c(seq(0.1, 0.5, by = 0.1), 0.01, 0.05)
  labels <- as.character(breaks)
  
  plt1 <-
    ggplot(data = plt_df, aes(x = percents, y = RDR, group = group)) +
    geom_line(aes(color = group), size = 3) +
    geom_point(aes(shape = group), size = 4, fill = "#153354") +
    labs(colour = "p-value", x = "Number of top CSPs") +
    theme_classic() +
    # theme(legend.position = "none") +
    theme(legend.position = c(0.8, 0.8)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 26)) +
    theme(legend.text = element_text(size = 22),
          legend.title = element_text(size = 22)) +
    ylim(c(0, 0.25)) +
    # scale_color_manual(values=c("#f97750", "#57b894"))+
    # scale_color_manual(values = c("#1f77b4", "#ed665d")) +
    scale_colour_manual(
      "p-value",
      labels = c(0.01, 0.05),
      values = c("#1f77b4", "#ed665d")
    ) +
    scale_shape_manual("p-value",
                       labels =  c(0.01, 0.05),
                       values = c(15, 17)) +
    
    geom_hline(
      yintercept = 0.01,
      linetype = "dashed",
      color = "#1f77b4",
      size = 1
    ) +
    geom_hline(
      yintercept = 0.05,
      linetype = "dashed",
      color = "#ed665d",
      size = 1
    ) +
    scale_y_continuous(
      limits = c(0, 0.25),
      breaks = breaks,
      labels = labels,
      name = "RDR"
    ) +
    theme(legend.key.width = unit(3, "line"))
  
  plt1
}

# figure 2A
draw_n_csps_figure <- function(gdsc_stat, tcga_stat) {
  cancers <- c(
    "BLCA",
    "BRCA",
    "CESC",
    "COAD/READ",
    "DLBC",
    "ESCA",
    "GBM",
    "HNSC",
    "KIRC",
    "LAML",
    "LGG",
    "LIHC",
    "LUAD",
    "LUSC",
    "MESO",
    "NB",
    
    "OV",
    "PAAD",
    "PRAD",
    "SKCM",
    "STAD",
    "THCA",
    "UCEC"
  )
  
  cancers <- cancers[order(cancers)]
  res <- lapply(cancers, function(cancer) {
    stat1 <- gdsc_stat[gdsc_stat$disease == cancer, ]
    stat2 <- tcga_stat[tcga_stat$disease == cancer, ]
    pick1 <- paste0(stat1$drug, "@", stat1$pathway)
    pick2 <- paste0(stat2$drug, "@", stat2$pathway)
    r <-
      c(nrow(stat1), nrow(stat2), length(intersect(pick1, pick2)))
  })
  res <- do.call(rbind, res)
  colnames(res) <-
    c("CSPs in GDSC", "CSPs in TCGA", "Validated CSPs")
  res <- data.frame(res)
  res$Cancer <- cancers
  res$ID <- seq(length(cancers))
  dat <- res
  dat$CSPs.in.GDSC.log <- log2(dat$CSPs.in.GDSC + 1)
  dat$CSPs.in.TCGA.log <- log2(dat$CSPs.in.TCGA + 1)
  dat$Validated.CSPs.log <- log2(dat$Validated.CSPs + 1)
  
  dat$Cancer[dat$Cancer == "LAML"] <- "AML"
  dat <- dat[order(dat$CSPs.in.GDSC, decreasing = TRUE), ]
  dat <- dat[order(dat$CSPs.in.TCGA, decreasing = TRUE), ]
  
  myorder <- order(dat$Validated.CSPs.log, decreasing = TRUE)
  
  data <- data.frame(
    name = c(dat$Cancer, dat$Cancer, dat$Cancer),
    value = c(
      dat$CSPs.in.GDSC.log,
      dat$Validated.CSPs.log,
      dat$CSPs.in.TCGA.log
    ),
    Group = c(
      rep("GDSC", nrow(dat)),
      rep("Validation", nrow(dat)),
      rep("TCGA", nrow(dat))
    )
  )
  data$name <- factor(data$name, levels = dat$Cancer[myorder])
  data$Group <-
    factor(data$Group, levels = c("GDSC", "TCGA", "Validation"))
  
  tickVal <- c(0, 5, 50, 100, 500, 1000, 5000)
  breakVal <- log2(tickVal)
  
  palettes <-
    ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
  palname <- names(palettes)[1]
  pal <- tableau_color_pal(palname)
  max_n <- attr(pal, "max_n")
  myPalette <- pal(max_n)
  
  
  colorVal <- c(myPalette[1], myPalette[2], myPalette[3])
  plt0 <- ggplot(data, aes(x = name, y = value, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab("Cancer") +
    ylab("Number of CSPs") +
    scale_fill_manual(values = colorVal) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title = element_text(size = 26),
      axis.text = element_text(size = 22),
      legend.title = element_blank(),
      legend.text = element_text(size = 22),
      legend.position = "top"
    ) +
    theme(legend.spacing.x = unit(2, "mm")) +
    scale_y_continuous(breaks = breakVal, labels = tickVal) +
    theme(axis.title.x = element_text(margin = margin(t = -45))) +
    theme(plot.margin = unit(c(0, 0, 3, 0), "lines"))
  plt0
}

# figure 2D
get_permutation_pas <- function(pas, drug, pathway) {
  drugs <- names(pas)
  pathways <- rownames(pas[[1]])
  pas_1 <- pas[[drug]][pathway, ]
  
  groups <- sapply(names(pas_1),
                   function(s) unlist(strsplit(s, "_"))[1])
  names(groups) <- NULL
  pick_aml <- which(groups == "LAML")
  
  pas_1 <- pas_1[pick_aml]
  
  pick_d <- which(drug == drugs)
  pas_2 <- sapply(drugs[-pick_d], function(d) {
    pas[[d]][pathway, pick_aml]
  })
  pas_2 <- as.vector(pas_2)
  
  pick_p <- which(pathway == pathways)
  pas_3 <- as.vector(pas[[drug]][-pick_p, pick_aml])
  
  pas_4 <- sapply(drugs[-pick_d], function(d) {
    pas[[d]][-pick_p, pick_aml]
  })
  pas_4 <- as.vector(pas_4)
  
  res <- list(pas_1, pas_2, pas_3, pas_4)
  res
}

draw_permutation_figure <- function(pas) {
  drug <- "Quizartinib"
  pathway <- "MARTENS_BOUND_BY_PML_RARA_FUSION"
  
  pas_list <- get_permutation_pas(pas, drug, pathway)
  df <- lapply(1:4, function(i) {
    b <- as.vector(pas_list[[i]])
    data.frame(PAS = as.vector(b),
               Drugs = i,
               col = "#57b894")
  })
  df <- do.call(rbind, df)
  pick <- which(df$Drugs == 1)
  df[pick,]$col <- "#f97750"
  
  df$Drugs = as.factor(df$Drugs)
  
  x_axis <- list()
  x_axis[[1]] <- bquote(atop(D[i] * "," * P[j], "(27)"))
  x_axis[[2]] <- bquote(atop(bar(D[i]) * "," * P[j], "(6,750)"))
  x_axis[[3]] <- bquote(atop(D[i] * "," * bar(P[j]), "(128,547)"))
  x_axis[[4]] <-
    bquote(atop(bar(D[i]) * "," * bar(P[j]), "(32,136,750)"))
  
  # this p-value is only for MARTEN_PML_RARA drugable by Quizartinib
  t_test_df <- data.frame(
    group1 = 1,
    group2 = 2:4,
    label = c("p-value = 1e-4", "1e-4", "1e-4"),
    y.position = c(125, 140, 155)
  )
  
  plt <- ggplot(data = df,
                lwd = 0.8,
                fatten = 1) +
    geom_boxplot(aes(Drugs, y = PAS, fill = col),
                 outlier.shape = NA) +
    coord_cartesian(ylim = c(-140, 160)) +
    scale_fill_manual(values = c("#d9e6eb", "#f97750")) +
    add_pvalue(t_test_df, tip.length = 0, label.size = 7) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") +
    labs(y = "PAS", x = "Drug, Pathway") +
    scale_x_discrete(labels = x_axis) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.title = element_text(size = 26))
  plt
}



draw_tstat_chistat_figure <- function(stat) {
  cancers <- c(
    "BLCA",
    "BRCA",
    "CESC",
    "COAD/READ",
    "DLBC",
    "ESCA",
    "GBM",
    "HNSC",
    "KIRC",
    "LAML",
    "LGG",
    "LIHC",
    "LUAD",
    "LUSC",
    "MESO",
    "NB",
    "OV",
    "PAAD",
    "PRAD",
    "SKCM",
    "STAD",
    "THCA",
    "UCEC"
  )
  stat <- stat[stat$disease %in% cancers, ]
  
  stat$color <- 0
  stat$alpha <- 0.9
  stat$label <- ""
  stat$shape <- 1
  stat$size <- 2
  stat$size <- as.integer(stat$size)
  
  plt_1 <- ggplot(stat, aes(
    x = as.numeric(chisq_stat),
    y = as.numeric(t_stat),
    label = label
  )) +
    geom_point(aes(
      color = as.factor(color),
      shape = as.factor(shape),
      size = size
    ), alpha = 0.3) +
    scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
    scale_size(range = c(1, 2)) +
    # scale_size_identity()+
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "chi2-statistics", y = "t-statistics") +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text.x = element_text(size = 18)) +
    theme(axis.text.y = element_text(size = 18))
  
  stat_aml <- stat[stat$disease == "LAML", ]
  plt_2 <- ggplot(stat_aml, aes(
    x = as.numeric(chisq_stat),
    y = as.numeric(t_stat),
    label = label
  )) +
    geom_point(aes(
      color = as.factor(color),
      shape = as.factor(shape),
      size = size
    ), alpha = 0.3) +
    scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
    scale_size(range = c(1, 2)) +
    # scale_size_identity()+
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "chi2-statistics", y = "t-statistics") +
    theme(axis.title = element_text(size = 18)) +
    theme(axis.text.x = element_text(size = 18)) +
    theme(axis.text.y = element_text(size = 18))
  
  plt <- ggarrange(
    plt_1,
    plt_2,
    nrow = 1,
    ncol = 2,
    common.legend = F,
    labels = c("A", "B"),
    font.label = list(size = 20)
  )
  plt
  # ggsave(filename = "aml_statistics.png", plt, width = 15, height = 7)
}

draw_pie_chart_figure <- function(tbl) {
  diseases <- tbl$cancer
  gdsc <- tbl$CSPs.in.GDSC
  tcga <- tbl$CSPs.in.TCGA
  gdsc <- gdsc / sum(gdsc)
  tcga <- tcga / sum(tcga)
  
  n <- 11
  qual_col_pals <-
    brewer.pal.info[brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(
    brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
  titles <- c("GDSC", "TCGA")
  
  input <- list(gdsc, tcga)
  
  paths <- lapply(1:2, function(i) {
    value <- input[[i]]
    title <- titles[i]
    data <- data.frame(group = diseases, value = value)
    data$value <- data$value * 100
    data <- data[order(data$value, decreasing = T), ]
    rownames(data) <- 1:nrow(data)
    other <- sum(data$value[11:23])
    data <- rbind(data[1:10, ], c("Other", other))
    data$value <- as.numeric(data$value)
    data$ypos <- cumsum(data$value) - 0.5 * data$value
    
    path <- paste0("output/piechart_", title, ".pdf")
    pdf(file = path,
        width = 20,
        height = 20)
    pie(
      data$value,
      labels = paste(data$group, " - ",
                     as.character(round(data$value, 2)),
                     "%",
                     sep = ""),
      col = col_vector[1:n],
      radius = 0.7,
      main = title,
      cex.main = 5,
      cex = 3
    )
    dev.off()
  })
  
}

draw_pas_figure <- function(gdsc_pas, tcga_pas) {
  groups <-
    sapply(names(gdsc_pas), function(s)
      unlist(strsplit(s, "_"))[1])
  names(groups) <- NULL
  
  df <- data.frame("pas" = gdsc_pas, "disease" = groups)
  df$pas <- log2(df$pas + 1)
  col <- rep("red", nrow(df))
  pick <- which(groups == "LAML")
  col[pick] <- "blue"
  df$col <- col
  
  pick <- which(groups == "LAML")
  df$disease[pick] <- "AML"
  
  cancer_order <- c(
    "BLCA",
    "BRCA",
    "CESC",
    "COAD/READ",
    "DLBC",
    "ESCA",
    "GBM",
    "HNSC",
    "KIRC",
    "AML",
    "LGG",
    "LIHC",
    "LUAD",
    "LUSC",
    "MESO",
    "NB",
    "OV",
    "PAAD",
    "PRAD",
    "SKCM",
    "STAD",
    "THCA",
    "UCEC"
  )
  
  df$disease <- factor(df$disease, cancer_order)
  
  plt_gdsc <- ggplot() +
    geom_boxplot(
      data = df,
      aes(x = disease,
          y = pas, fill = col),
      color = "black",
      lwd = 0.6,
      fatten = 1
    ) +
    scale_fill_manual(values = c("#f97750", "#57b894")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 18,
      vjust = 0.7
    )) +
    theme(axis.text.y = element_text(size = 18)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 18)) +
    theme(legend.position = "none") +
    ylab("PAS") +
    xlab("Cancer")
  
  groups <-
    sapply(names(tcga_pas), function(s)
      unlist(strsplit(s, "_"))[1])
  names(groups) <- NULL
  pick <- which(groups %in% c("COAD", "READ"))
  groups[pick] <- "COAD/READ"
  
  df <- data.frame("pas" = tcga_pas, "disease" = groups)
  df$pas <- log2(df$pas + 1)
  col <- rep("red", nrow(df))
  pick <- which(groups == "LAML")
  col[pick] <- "blue"
  df$col <- col
  
  pick <- which(groups == "LAML")
  df$disease[pick] <- "AML"
  
  df$disease <- factor(df$disease, level = cancer_order)
  
  plt_tcga <- ggplot() +
    geom_boxplot(
      data = df,
      aes(x = disease,
          y = pas, fill = col),
      color = "black",
      lwd = 0.6,
      fatten = 1
    ) +
    scale_fill_manual(values = c("#f97750", "#57b894")) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      size = 18,
      vjust = 0.7
    )) +
    theme(axis.text.y = element_text(size = 18)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 18)) +
    theme(legend.position = "none") +
    ylab("PAS") +
    xlab("Cancer")
  
  plt <- ggarrange(
    plt_gdsc,
    plt_tcga,
    nrow = 2,
    ncol = 1,
    common.legend = F,
    labels = c("A", "B"),
    font.label = list(size = 20)
  )
  plt
}

nscore <- function(x) {
  n <- length(x)
  rx <- rank(x)
  ns <- qnorm(rx / (n + 1))
  return(ns)
}

draw_pas_auc_fig <- function(pas_auc) {
  df <- data.frame(PAS = nscore(pas_auc$gdsc$pas),
                   AUC = nscore(pas_auc$gdsc$auc))
  plt_gdsc <- ggplot(df, aes(x = AUC, y = PAS)) +
    geom_point(size = 5,
               col = "#3677a3",
               alpha = 0.8) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title = element_text(size = 18)) +
    stat_smooth(
      method = "lm",
      se = FALSE,
      color = "#f97750",
      size = 5
    )
  
  df <- data.frame(PAS = nscore(pas_auc$beataml$pas),
                   AUC = nscore(pas_auc$beataml$auc))
  plt_beataml <- ggplot(df, aes(x = AUC, y = PAS)) +
    geom_point(size = 5,
               col = "#3677a3",
               alpha = 0.8) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 30)) +
    theme(axis.text.y = element_text(size = 30)) +
    theme(axis.title = element_text(size = 18)) +
    stat_smooth(
      method = "lm",
      se = FALSE,
      color = "#f97750",
      size = 5
    )
  
  plt <- ggarrange(
    plt_gdsc,
    plt_beataml,
    nrow = 1,
    ncol = 2,
    common.legend = F,
    labels = c("A", "B"),
    font.label = list(size = 20)
  )
}
