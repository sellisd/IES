# prepare data for revBayes
library(reshape2)
library(ape)
source("sharedFunctions.R")
load("~/data/IES_data/rdb/charMats")
clusters <- unique(charMats$cluster)
geneFamiliesProcessed <- character()
for(i in c(1:length(clusters))){
  cat(i, "/", length(clusters), "\r")
  nex <- save2nexus(clusters[i])
  #check if we have a tree file
  dataNexusOut <- paste0("~/data/IES_data/msas/asr/charMat",clusters[i],".nexus")
  treeFileIn <- paste0("~/data/IES_data/msas/phyldog/results/",clusters[i],".ReconciledTree")
  if(file.exists(treeFileIn)){
    # write a nexus data file for each gene family
    modif.write.nexus.data(nex, file = dataNexusOut, format = "Standard")
    # and strip the reconciled tree from NHX annotations
    geneTreeString <- readLines(con = treeFileIn) #read nhx tree as string, remove [] and read as newick formated string
    geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
    geneTree <- read.tree(text = geneTreeString)
    # we could nicely plot the trees also
    # source("~/projects/IES/src/sharedPlotFunctions.R")
    # cl <- colBySpec(geneTree)
    # plot(geneTree, tip.col = cl)
    write.nexus(geneTree, file = paste0("~/data/IES_data/msas/asr/tree",clusters[i],".nexus"))
    geneFamiliesProcessed <- append(geneFamiliesProcessed, clusters[i])
  }else{
    # skip gene family
  }
}

write.table(geneFamiliesProcessed, file = "~/data/IES_data/msas/asr/geneFamilies.dat", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
