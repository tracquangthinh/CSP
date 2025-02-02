# Pathway Activation Score (PAS) Generation

## Requirements

- R 3.6 or newer

- R packages: Rfast, fastmatch, parallel

- Download `data_pas` from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6787033.svg)](https://doi.org/10.5281/zenodo.6787033), extract the `data` directory to the main directory. The compressed file contains gene expression, zscore obtained from network enrichment analysis (NEA) and up/down-stream gene sets. Note that gene expression contains 2 demo samples.

## Execution

- The pipeline is executed in parallel mode. You can change number of cores by modifying `n_cores`

```
Rscript run.R
```
or

```
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
```
