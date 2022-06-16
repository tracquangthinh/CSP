library(KEGGgraph)
library(org.Hs.eg.db)
# to install geneSynonym if needed
# library(devtools)
# install_github('oganm/geneSynonym')
library(geneSynonym)

print("Generating target matrix ...")
source("generate_target_matrix.R")
print("Getting directed pathways ...")
source("get_directed_pathway.R")
print("Getting up/downstream genes ...")
source("get_up_downstream.R")

# output: target_updown_stream.RData
