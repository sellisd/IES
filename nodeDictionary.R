# make a node dictionary linking the node numbering in revBayes, phyldog and ape (R)
library(ape)
source("~/projects/IES/src/sharedFunctions.R")

for(asrRun in c(1,2,3)){
  rbPath <- paste0("~/data/IES/analysis/asr", asrRun, "/rbNodeIndexes/")
  phyldogPath <- paste0("~/data/IES/analysis/phyldogT", asrRun, "/results/")
  fileNames <- dir(path = rbPath, pattern = "^nodeIndex.\\d+.tre$")
  clusters <- substr(fileNames, 11, (nchar(fileNames) - 4))
  m <- matrix(nrow = 0, ncol = 4)
  counter <- 0
  for(cluster in clusters){
    cat(paste(counter , "/", length(clusters), "\r"))
    tr <- read.tree(paste("~/data/IES/analysis/phyldogT", asrRun, "/results/",cluster,".ReconciledTree",sep=""))
    ktr <- read.table(paste("~/data/IES/analysis/phyldogT", asrRun, "/results/",cluster,".ReconciledTree.key",sep=""))
    dictPhyldog <- linkNodes(tr, ktr)
    tr <- read.tree(paste("~/data/IES/analysis/asr", asrRun, "/rbNodeIndexes/nodeIndex.",cluster,".tre",sep=""))
    ktr <- read.table(paste("~/data/IES/analysis/asr", asrRun, "/rbNodeIndexes/nodeIndex.",cluster,".tre.key",sep=""))
    dictRB <- linkNodes(tr,ktr)
    if(!all(dictPhyldog[, 1] == dictRB[, 1])){
      stop("error with node dictionaries!")
    }
    m <- rbind(m, cbind(cluster, dictPhyldog[, 1], dictPhyldog[, 2], dictRB[, 2]))
    counter <- counter + 1
  }
  nodeDictionary <- data.frame(m, stringsAsFactors = FALSE)
  names(nodeDictionary) <- c("cluster", "r", "phyldog", "rb")
  #save(nodeDictionary, file = "~/data/IES_data/rdb/nodeDictionary")
  #also save it as text
  write.table(nodeDictionary, file = paste0("~/data/IES/analysis/tables/nodeDictionary", asrRun, ".dat"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}