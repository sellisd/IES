# calculate statistics on trees, the total branch length, number of genes and whether they have a character matrix
library(ape)
treePath <- "~/data/IES_data/msas/phyldog/results/"
load("~/data/IES_data/rdb/charMats")
clusters <- dir(path = treePath, pattern = "*.ReconciledTree$")
clusters <- gsub(pattern = ".ReconciledTree", replacement = "", clusters, fixed = TRUE)
n <- length(clusters)
treeStatsDF <- data.frame(cluster = character(n), branchLength = numeric(n), ngenes = numeric(n), hasCharMat = logical(n), stringsAsFactors = FALSE)
counter <- 1
for(cluster in clusters){
  cat(counter, "/", length(clusters),"\r")
  tr <- read.tree(file = paste0(treePath, cluster, ".ReconciledTree"))
  ngenes <- length(tr$tip.label)
  totalBrLength <- sum(tr$edge.length)
  treeStatsDF[counter, "cluster"] <- cluster
  treeStatsDF[counter, "branchLength"] <- totalBrLength
  treeStatsDF[counter, "ngenes"] <- ngenes
  if(cluster %in% charMats$cluster[charMats$ies!=0]){
    # there is a character matrix with at least one IES
    treeStatsDF[counter, "hasCharMat"] <- TRUE
  }else{
    treeStatsDF[counter, "hasCharMat"] <- FALSE
  }
  counter <- counter + 1
}
save(treeStatsDF, file = "~/data/IES_data/rdb/treeStats")
