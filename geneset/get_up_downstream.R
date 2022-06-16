## expand drug targets
## to include upstream and downstream components
##

## directed networks
NET <- list()
NET[[1]] <- read.table("data/HTRIdb.TF-TG.txt", sep = "\t", as.is = TRUE)
NET[[2]] <- read.table("data/TF2genes.aliases.MSigDB_2013", header = TRUE, sep = "\t", as.is = TRUE)
NET[[2]][, 1] <- toupper(NET[[2]][, 1])
NET[[2]][, 2] <- toupper(NET[[2]][, 2])
NET[[3]] <- read.table("data/Kinase_Substrate.2015.human", sep = "\t", as.is = TRUE)
NET[[4]] <- read.table("data/TF_GenomeUCSC_z2.2_Factors2targets.5p100000.3p5000.human",
  sep = "\t", as.is = TRUE
)


## KEGG pathway

kegg <- readRDS("data/gene_path_link")
pick <- sapply(kegg, length) > 0
kegg <- kegg[pick]
a <- Map(cbind, names(kegg), kegg)
NET[[5]] <- Reduce(rbind, a)

load("data/target_objects.RData")


## synonyms
alltarget <- sapply(tgt_lst, function(x) unlist(humanSyno(x)))

## up-downstream of drug targets, using directed net
down <- function(gene, net) {
  pick <- net[, 1] %in% gene
  if (sum(pick) == 0) down <- NULL
  if (sum(pick) > 0) down <- net[pick, 2]
  return(unique(down))
}
up <- function(gene, net) {
  pick <- net[, 2] %in% gene
  if (sum(pick) == 0) up <- NULL
  if (sum(pick) > 0) up <- net[pick, 1]
  return(unique(up))
}

## Kinases i=3
i <- 3
tgt_down_kns <- sapply(alltarget, down, net = NET[[i]])
tgt_up_kns <- sapply(alltarget, up, net = NET[[i]])
## KEGG
i <- 5
tgt_down_keg <- sapply(alltarget, down, net = NET[[i]])
tgt_up_keg <- sapply(alltarget, up, net = NET[[i]])

## Trans Fac
tgt_down_tf <- tgt_up_tf <- sapply(alltarget, function(x) NULL)
for (i in c(1, 2, 4)) {
  tgt_down_tf <- Map(union, tgt_down_tf, sapply(alltarget, down, net = NET[[i]]))
  tgt_up_tf <- Map(union, tgt_up_tf, sapply(alltarget, up, net = NET[[i]]))
}
## ALL
tgt_down_all <- tgt_up_all <- sapply(alltarget, function(x) NULL)
for (i in 1:5) {
  tgt_down_all <- Map(union, tgt_down_all, sapply(alltarget, down, net = NET[[i]]))
  tgt_up_all <- Map(union, tgt_up_all, sapply(alltarget, up, net = NET[[i]]))
}

## including synonyms for the targets
target <- alltarget
## clear the overlaps with targets
tgt <- target
tgt_down_kns <- Map(setdiff, tgt_down_kns, tgt)
tgt_down_keg <- Map(setdiff, tgt_down_keg, tgt)
tgt_down_tf <- Map(setdiff, tgt_down_tf, tgt)
tgt_down_all <- Map(setdiff, tgt_down_all, tgt)
tgt_up_kns <- Map(setdiff, tgt_up_kns, tgt)
tgt_up_keg <- Map(setdiff, tgt_up_keg, tgt)
tgt_up_tf <- Map(setdiff, tgt_up_tf, tgt)
tgt_up_all <- Map(setdiff, tgt_up_all, tgt)

save(target, tgt_up_kns, tgt_down_kns,
  tgt_down_keg, tgt_up_keg,
  tgt_down_tf, tgt_up_tf,
  tgt_down_all, tgt_up_all,
  file = "target_updown_stream.RData"
)
