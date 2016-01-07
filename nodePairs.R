# find pairs of nodes where we should calculate gain and loss events
library(phangorn)
library(ggtree)
load("~/data/IES_data/rdb/nodeDictionary")
source("~/projects/IES/src/sharedFunctions.R")
phyldogPath <- "~/data/IES_data/msas/phyldog/results/"
# pairs of speciation events r: regular, s: skipping node
branches <- data.frame(from = c(0, 0, 2, 2, 4, 4, 0, 0, 0, 0, 2, 2),
                       to   = c(1, 2, 3, 4, 6, 5, 3, 4, 5, 6, 5, 6),
                       type = c(rep("r", 6), rep("s", 6)))
clusters <- read.table("~/data/IES_data/msas/asr/geneFamilies.dat", stringsAsFactors = FALSE, as.is = TRUE)
counter <- 1
nodePairs <- data.frame(cluster   = numeric(0),
                        fromP     = numeric(0),
                        toP       = numeric(0),
                        fromEvent = numeric(0),
                        toEvent   = numeric(0),
                        fromR     = numeric(0),
                        toR       = numeric(0),
                        stringsAsFactors = FALSE)
for(cluster in clusters$V1){
  print(paste(counter, "/", length(clusters$V1)))
  #find from which to which nodes transitions are of interest
  ph <- read.nhx(paste0(phyldogPath, cluster, ".ReconciledTree"))
  #find node pairs in tree that correspond to events.
  nP <- data.frame(cluster = numeric(0), fromP = numeric(0), toP = numeric(0), fromEvent = numeric(0), toEvent = numeric(0))
  for(i in c(1:nrow(branches))){
    nP <- rbind(nP, getNodePairs(ph, cluster, branches[i, "from"], branches[i, "to"]))
  }
  # translate to R id notation
  nP <- data.frame(nP,
                   fromR = phyldog2r(nP$fromP, cluster = cluster),
                   toR = phyldog2r(nP$toP, cluster = cluster),
                   stringsAsFactors = FALSE)
  nodePairs <- rbind(nodePairs, nP)
  counter <- counter + 1
}
#save(nodePairs, file = "~/data/IES_data/rdb/nodePairs")
write.table(nodePairs, file = "~/data/IES_data/rdb/nodePairs.dat", sep = "\t", quote = FALSE, row.names = FALSE)
