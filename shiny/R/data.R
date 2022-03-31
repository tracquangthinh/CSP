lazyLoad("data/gdsc_pas_boxplot_env")
lazyLoad("data/tcga_pas_boxplot_env")

gdsc_stat <- readRDS("data/gdsc_stat")

cancer_abbreviations <- readRDS("data/cancer_abbreviations")

drugs <- unique(gdsc_stat$drug)
drugs <- drugs[order(drugs)]

cancers <- c(
   "BLCA", "BRCA", "CESC", "COAD/READ",
   "DLBC", "ESCA", "GBM", "HNSC", "KIRC", "LAML",
   "LGG", "LIHC", "LUAD", "LUSC", "MESO", "NB",
   "OV", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC"
)

top_csps <- readRDS("data/top_csps")
rownames(top_csps) <- NULL

beatAML_drug_description <- readRDS("data/beataml_drug_description")

cor_pas_auc <- readRDS("data/gdsc_aml_cor_pas_auc")

top_pas_auc <- readRDS("data/pas_auc_for_top_csps")
