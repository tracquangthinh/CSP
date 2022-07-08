z_transform <- function(cors) {
  sapply(cors, function(c) {
    0.5 * (log(1 + c) - log(1 - c))
  })
}

draw_rdr_figure <- function(cor_pas_auc, left_tail = TRUE) {
  pick <- c(
    which(is.na(cor_pas_auc$cor_gdsc)),
    which(is.na(cor_pas_auc$cor_beataml))
  )
  if(length(pick) > 0) {
    cor_pas_auc <- cor_pas_auc[-pick, ]
  }
  t_gdsc <- cor_pas_auc$cor_gdsc
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
    geom_point(aes(shape = group), size = 4, fill = "#153354") +
    labs(colour = "p-value", x = "Number of top CSPs") +
    theme_classic() +
    # theme(legend.position = "none") +
    theme(legend.position = c(0.8, 0.8)) +
    theme(axis.text.x = element_text(size = 22)) +
    theme(axis.text.y = element_text(size = 22)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title = element_text(size = 26)) +
    theme(legend.text = element_text(size = 22), legend.title = element_text(size = 22)) +
    ylim(c(0, 0.25)) +
    # scale_color_manual(values=c("#f97750", "#57b894"))+
    # scale_color_manual(values = c("#1f77b4", "#ed665d")) +
    scale_colour_manual("p-value", labels = c(0.01, 0.05), values = c("#1f77b4", "#ed665d")) +
    scale_shape_manual("p-value", labels =  c(0.01, 0.05), values = c(15, 17)) +

    geom_hline(yintercept = 0.01, linetype = "dashed", color = "#1f77b4", size = 1) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "#ed665d", size = 1) +
    scale_y_continuous(
      limits = c(0, 0.25), breaks = breaks, labels = labels,
      name = "RDR"
    ) +
    theme(legend.key.width = unit(3, "line"))

  plt1
}

