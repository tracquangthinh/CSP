library(ggplot2)
library(ggpubr)

source("functions.R")

gdsc_aml_cor_pas_auc <- readRDS("cor_pas_auc")
negative_cor_fig <- draw_rdr_figure(gdsc_aml_cor_pas_auc, left_tail = TRUE)
positive_cor_fig <- draw_rdr_figure(gdsc_aml_cor_pas_auc, left_tail = FALSE)

plt <- ggarrange(negative_cor_fig, positive_cor_fig,
  nrow = 1, ncol = 2, common.legend = F,
  labels = c("A", "B"), font.label = list(size = 20)
)

ggsave(filename = "rdr_fig.png", plt, width = 15, height = 7)
