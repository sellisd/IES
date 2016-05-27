spEvents <- read.table("~/data/IES/analysis/tables/spEvents.dat", header = TRUE)
# Find the speciation events where gene tree nodes are inferred to have an IES (with p > 0.99)
cutoff <- 0.99
asrs <- read.table("~/data/IES/analysis/tables/avNodeProb.dat", header = TRUE) # nodes are in rb index notation
haveIES <- asrs[which(asrs$presence > cutoff), ]
keys <- paste0(haveIES$cluster, '.', haveIES$iesColumn)
cat("cluster", "iesColumn", "spEvents", "\n")
for(key in unique(keys)){
  # key <- unique(keys)[1]
  index <- which(keys == key)
  nodes <- haveIES$node[index]
  cluster <- haveIES$cluster[index]
  iesColumn <- haveIES$iesColumn[index]
  if(all(cluster == cluster[1]) & all(iesColumn == iesColumn[1])){
    E <- unique(spEvents[spEvents$cluster == cluster[1] & spEvents$nodeP %in% nodes, "spEvent"])
    if(length(E) == 0){
      # no IES present with p>cutoff in any speciation node
    }else{
      cat(cluster[1], iesColumn[1], E, "\n")
    }
  }else{
    stop()
  }
}
