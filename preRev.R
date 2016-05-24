# prepare data for revBayes
library(reshape2)
library(ape)
source("~/projects/IES/src/sharedFunctions.R")
charMats <- read.table("~/data/IES/analysis/iesdb/charMats.tab", stringsAsFactors = FALSE, header = TRUE)
#load("~/data/IES_data/rdb/charMats")
clusters <- unique(charMats$cluster)
geneFamiliesProcessed <- character()
skippedFamilies <- character()
for(i in c(1:length(clusters))){
  cat(i, "/", length(clusters), "\r")
  nex <- save2nexus(clusters[i])
  spNames <- gene2species(names(nex))
  # if data matrix has T. thermophila exclude list elements
  #nex <- nex[which(spNames != "Tetrahymena_thermophila")]
 #check if we have a tree file
  dataNexusOut <- paste0("~/data/IES/analysis/asr/charMat",clusters[i],".nexus")
  treeFileIn <- paste0("~/data/IES/analysis/phyldog/results/",clusters[i],".ReconciledTree")
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
    if (!identical(sort(geneTree$tip.label), sort(gene2protName(names(nex))))){
      stop(paste("Genes do not match between tree and character matrix!: ", clusters[i], i))
    }
    write.nexus(geneTree, file = paste0("~/data/IES/analysis/asr/tree",clusters[i],".nexus"))
    geneFamiliesProcessed <- append(geneFamiliesProcessed, clusters[i])
  }else{
    skippedFamilies <- append(skippedFamilies, clusters[i])
  }
}

write.table(geneFamiliesProcessed, file = "~/data/IES/analysis/asr/geneFamilies.dat", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("skipped families")
print(skippedFamilies)