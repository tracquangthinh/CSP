load("data/target_objects.RData")

rm(tgt_mat)

x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
rm(x, mapped_genes)


paths_lst <- lapply(tgt_lst, function(ii) {
  as.vector(na.omit(unique(unlist(mget(unique(unlist(xx[ii])), org.Hs.egPATH)))))
})

paths_sy_lst <- lapply(tgt_sy_lst, function(ii) {
  as.vector(na.omit(unique(unlist(mget(unique(unlist(xx[ii])), org.Hs.egPATH)))))
})


un_paths <- unique(unlist(paths_lst))

tmp <- tempfile()

paths_lst <- lapply(un_paths, function(pId) {
  cat(pId, "\n")
  retrieveKGML(pId, organism="hsa", destfile=tmp, method="wget", quiet=TRUE)
  # for windows
  # try(retrieveKGML(pId, organism = "hsa", destfile = tmp, method = "internal", quiet = TRUE))
  try(parseKGML2Graph(tmp, expandGenes = TRUE, genesOnly = TRUE))
})
names(paths_lst) <- un_paths


nodes_edges <- lapply(paths_lst[sapply(paths_lst, class) == "graphNEL"], edges)
names(nodes_edges) <- NULL
u_names <- unique(unlist(lapply(nodes_edges, names))) ## 6086

eg_list <- vector(length = length(u_names), mode = "list")
names(eg_list) <- u_names

for (i in seq(along = nodes_edges)) {
  lst <- nodes_edges[[i]]
  for (j in seq(along = lst)) {
    if (length(lst[[j]])) {
      eg_list[[names(lst)[j]]] <- unique(c(eg_list[[names(lst)[j]]], gsub("hsa:", "", lst[[j]])))
    }
  }
}

rm(i, j, lst)

names(eg_list) <- gsub("hsa:", "", names(eg_list))

eg2sym <- unlist(mget(names(eg_list), org.Hs.egSYMBOL, ifnotfound = as.list(NA)))

eg_list <- eg_list[!is.na(eg2sym)]
eg2sym <- eg2sym[!is.na(eg2sym)]

names(eg_list) <- eg2sym

xS <- org.Hs.egSYMBOL
xS <- as.list(xS[mappedkeys(xS)])


sym_list <- vector(length = length(eg_list), mode = "list")
names(sym_list) <- names(eg_list)

for (i in seq(along = eg_list))
{
  sym_list[[i]] <- unname(unlist(xS[eg_list[[i]]]))
}

saveRDS(sym_list, file = "data/gene_path_link")
