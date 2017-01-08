# make a node dictionary linking the node numbering in revBayes, phyldog and ape (R)
library(ape)
setwd('.')
source("sharedFunctions.R")
opt <- loadUserOptions()
basePath <- opt["basePath", ]

for(asrRun in c(1,2,3)){
  rbPath <- file.path(basePath, 'analysis', paste0("asr", asrRun), 'rbNodeIndexes')
  phyldogPath <- file.path(basePath, 'analysis', paste0("phyldogT", asrRun), "results")
  fileNames <- dir(path = rbPath, pattern = "^nodeIndex.\\d+.tre$")
  clusters <- substr(fileNames, 11, (nchar(fileNames) - 4))
  m <- matrix(nrow = 0, ncol = 4)
  counter <- 0
  for(cluster in clusters){
    cat(paste(counter , "/", length(clusters), "\r"))
    tr  <- read.tree(file.path(phyldogPath, paste0(cluster,".ReconciledTree")))
    ktr <- read.table(file.path(phyldogPath, paste0(cluster,".ReconciledTree.key")))
    dictPhyldog <- linkNodes(tr, ktr)
    tr  <- read.tree(file.path(rbPath, paste0("nodeIndex.",cluster,".tre")))
    ktr <- read.table(file.path(rbPath, paste0("nodeIndex.",cluster,".tre.key",sep="")))
    dictRB <- linkNodes(tr,ktr)
    if(!all(dictPhyldog[, 1] == dictRB[, 1])){
      stop("error with node dictionaries!")
    }
    m <- rbind(m, cbind(cluster, dictPhyldog[, 1], dictPhyldog[, 2], dictRB[, 2]))
    counter <- counter + 1
  }
  nodeDictionary <- data.frame(m, stringsAsFactors = FALSE)
  names(nodeDictionary) <- c("cluster", "r", "phyldog", "rb")
  write.table(nodeDictionary, file = file.path(basePath, 'analysis', 'tables', paste0("nodeDictionary", asrRun, ".dat")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}