# prepare data for revBayes
library(reshape2)
library(ape)

source("sharedFunctions.R")
load("~/data/IES_data/rdb/charMats")
clusters <- unique(charMats$cluster)
for(i in c(1:clusters)){
  cat(i, "/", length(clusters), "\r")
  nex <- save2nexus(clusters[i])
  #check if we have a tree file
  dataNexusOut <- paste0("~/data/IES_data/msas/asr/charMat",clusters[i],".nexus")
  treeFileIn <- paste0("~/data/IES_data/msas/phyldog/results/",clusters[i],".ReconciledTree")
  if(file.exists(treeFileIn)){
    # write a nexus data file for each gene family
    modif.write.nexus.data(nex, file = dataNexusOut, format = "Standard")
    # and strip the reconciled tree from NHX annotations
    geneTreeString <- readLines(con = treeFileIn)
    geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
    geneTree <- read.tree(text = geneTreeString)
    write.nexus(geneTree, file = paste0("~/data/IES_data/msas/asr/tree",clusters[i],".nexus"))
  }else{
    # skip gene family
  }
}


#source("~/projects/IES/src/sharedPlotFunctions.R")
#read nhx tree as string, remove [] and read as newick formated string
geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
geneTree <- read.tree(text = geneTreeString)
#cl <- colBySpec(geneTree)
#plot(geneTree, tip.col = cl)
write.nexus(geneTree, file = "~/projects/IES/src/trb/tree1.nexus")
cl <- colBySpec(geneTree)
plot(geneTree, tip.col = cl)
write.nexus(geneTree, file = "~/projects/IES/src/trb/tree2.nexus")
```
#nex10000 <- save2nexus(10000)
#ex10009 <- save2nexus(10009)
#list4nexus
#row.names(boolM) <- DFcasted$geneId
#write.table(matrix(c(DF$geneId, DF$ies), ncol = 2), file = "~/projects/IES/src/trb/charMat.tsv", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#knitr::kable(DFcasted)
#modif.write.nexus.data(nex10000, file = "~/projects/IES/src/trb/charMat1.nexus", format = "Standard")
#modif.write.nexus.data(nex10009, file = "~/projects/IES/src/trb/charMat2.nexus", format = "Standard")
