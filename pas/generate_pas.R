library(fastmatch)
library(Rfast)
library(parallel)

fn_gsum <- function(gene, xdat, genes, p1_, p2_, p3_) {
  pick <- gene %in% genes
  gene <- gene[pick]
  if (length(gene) == 0) return(rep(0, ncol(xdat)))
  pick <- fmatch(gene, genes)
  xdat <- xdat[pick, , drop = FALSE]
  pnew_ <- 1 + p1_ + p2_ + p3_
  xdat <- Rfast::eachrow(xdat, pnew_, "*")
  gsum <- colSums(xdat)
  return(gsum)
}

get_pas <- function(drug, exp, p1, p2, p3,
                    up_gene_set, down_gene_set) {
  cat(drug, " ")
  genes <- rownames(exp)

  ingene <- up_gene_set[[drug]]
  pathways <- names(ingene)
  up_pas <- sapply(pathways, function(pathway) {
    p1_ <- p1[[pathway]][drug]
    p2_ <- p2[[drug]]
    p3_ <- p3[[pathway]]
    gene <- ingene[[pathway]]
    fn_gsum(gene, exp, genes, p1_, p2_, p3_)
  })
  up_pas <- matrix(up_pas, ncol = ncol(exp), byrow = TRUE)

  ingene <- down_gene_set[[drug]]
  down_pas <- sapply(pathways, function(pathway) {
    p1_ <- p1[[pathway]][drug]
    p2_ <- p2[[drug]]
    p3_ <- p3[[pathway]]
    gene <- ingene[[pathway]]
    fn_gsum(gene, exp, genes, p1_, p2_, p3_)
  })
  down_pas <- matrix(down_pas, ncol = ncol(exp), byrow = TRUE)

  colnames(up_pas) <- colnames(down_pas) <- colnames(exp)
  rownames(up_pas) <- rownames(down_pas) <- pathways

  return(list(
    "up_pas" = up_pas,
    "down_pas" = down_pas
  ))
}

# main function to generate pathway activation score (PAS)
# INPUT
### gex: gene expressions
### z1, z2, z3: zscore obtained from network enrichment analysis (NEA)
### z1: Association between target genes and pathway genes
### z2: Association between driver genes and target genes
### z3: Association between pathway genes and driver genes
### geneset: contains up and down-stream genes
### For the structure of input, please check the data example
# OUTPUT
### pas$up: PAS for up-stream activity
### pas$down: PAS for down-stream activity

generate_pas <- function(gex, z1, z2, z3, geneset, n_cores = 2) {
  up_gene_set <- geneset$up
  down_gene_set <- geneset$down
  drugs <- names(up_gene_set)

  p1 <- lapply(z1, function(z) pnorm(z))
  p2 <- lapply(z2, function(z) pnorm(z))
  p3 <- lapply(z3, function(z) pnorm(z))

  pas_res <- mclapply(drugs,
    function(drug, exp, p1, p2, p3, up, down) {
      get_pas(drug, exp, p1, p2, p3, up, down)
    },
    exp = gex,
    p1 = p1,
    p2 = p2,
    p3 = p3,
    up = up_gene_set,
    down = down_gene_set,
    mc.preschedule = FALSE, mc.cores = n_cores
  )

  pas_up <- lapply(pas_res, function(l) l[[1]])
  pas_down <- lapply(pas_res, function(l) l[[2]])
  names(pas_down) <- names(pas_up) <- drugs
  pas <- list(
    "up" = pas_up,
    "down" = pas_down
  )

  return(pas)
}
