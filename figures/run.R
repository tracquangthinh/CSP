library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggprism)
library(RColorBrewer)

source("functions.R")

# Figure 1B, C (Supplementary Figure S3)
gdsc_pas <- readRDS("data/marten_gdsc_pas")
tcga_pas <- readRDS("data/marten_tcga_pas")
figure_1bc <- draw_pas_figure(gdsc_pas, tcga_pas)
ggsave(figure_1bc, filename = "output/figure_1bc.png", width = 15, height = 15)

# Figure 1D, E
pas_auc <- readRDS("data/marten_pas_auc")
figure_1de <- draw_pas_auc_fig(pas_auc)
ggsave(figure_1de, file = "output/figure_1de.png", width = 15, height = 10)

# Figure 2A
gdsc_stat <- readRDS("data/gdsc_stat")
tcga_stat <- readRDS("data/tcga_stat")

figure_2a <- draw_n_csps_figure(gdsc_stat, tcga_stat)
ggsave(figure_2a, file = "output/figure_2a.png", width = 23, height = 10)


# Figure 2B, C
cor_pas_auc <- readRDS("data/cor_pas_auc")
figure_2b <- draw_rdr_figure(cor_pas_auc, TRUE)
figure_2c <- draw_rdr_figure(cor_pas_auc, FALSE)
figure_2bc <- ggarrange(figure_2b, figure_2c,
                 nrow = 1, ncol = 2, common.legend = F,
                 labels = c("A", "B"), font.label = list(size = 20)
)

ggsave(figure_2bc, file = "output/figure_2bc.png", width = 15, height = 10)

# Figure 2D
# this figure requires GDSC PAS
# you need to download the GDSC PAS from Zenodo reposititory
# figure_2d <- draw_permutation_figure("/link/to/gdsc_pas")
# ggsave(figure_2d, file = "output/figure_2d.png", width = 10, height = 10)


# Supplementary Figure S1
figure_s1 <- draw_tstat_chistat_figure(gdsc_stat)
ggsave(figure_s1, file = "output/figure_s1.png", width = 18, height = 10)


# Supplementary Figure S2
summary_table <- readRDS("data/summary_table")
# Out: output/GDSC.pdf; output/TCGA.pdf
draw_pie_chart_figure(summary_table)

# Supplementary Figure S4
negative_cor_pas_auc <- readRDS("data/negative_cor_pas_auc")
figure_s4a <- draw_rdr_figure(negative_cor_pas_auc, TRUE)
figure_s4b <- draw_rdr_figure(negative_cor_pas_auc, FALSE)
figure_s4 <- ggarrange(figure_s4a, figure_s4b,
                 nrow = 1, ncol = 2, common.legend = F,
                 labels = c("A", "B"), font.label = list(size = 20)
             )
ggsave(figure_s4, file = "output/figure_s4.png", width = 17, height = 10)
