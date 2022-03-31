source("generate_pas.R")

### z1: Association between target genes and pathway genes
### z2: Association between driver genes and target genes
### z3: Association between pathway genes and driver genes
gex <- readRDS("data/gex")
z1 <- readRDS("data/zscore_target_pathway")
z2 <- readRDS("data/zscore_driver_target")
z3 <- readRDS("data/zscore_pathway_driver")
geneset <- readRDS("data/geneset")

pas <- generate_pas(gex, z1, z2, z3, geneset, n_cores = 2)
