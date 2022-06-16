# library(devtools)
# install_github('oganm/geneSynonym')
# library(geneSynonym)

ll <- readRDS("data/drug_gene")$target
ll <- lapply(ll, function(x) gsub("^[[:space:]]", "", x))
ll <- lapply(ll, function(x) gsub("[[:space:]]$", "", x))
ll <- lapply(ll, function(x) gsub("^[[:space:]]", "", x))
llu <- sort(unique(unlist(ll)))
mat <- matrix(0, ncol = length(llu), nrow = length(ll))
rownames(mat) <- names(ll)
colnames(mat) <- llu

for (id in seq(along = ll)) {
  mat[names(ll)[id], ll[[id]]] <- 1
}
tgt_mat <- mat
tgt_lst <- ll
## Expand using synonyms
tgt_sy_lst <- lapply(tgt_lst, function(ii) {
  unique(unlist(humanSyno(ii)))
})

save(tgt_mat, tgt_lst, tgt_sy_lst, file = "data/target_objects.RData")
