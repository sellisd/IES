# select gene trees with the same topology as the species tree
library(ape)
selectedGroups <- numeric(0)
treePath <- "~/data/IES/analysis/phyldog/results/"
#charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(treePath,pattern="[^.]*.ReconciledTree$")

spTreeString <- "(PCA:1,((PSE:1,PSO:1):1,(((PPE:1,PPR:1):1,(PBI:1,(POC:1,PTE:1):1):1):1,PTR:1):1):1):1;"

spTree <- read.tree(text=spTreeString)
spGT <- data.frame(stringsAsFactors = FALSE) # data frame with details of gene trees ~ species trees
counter <- 1
for(fileName in files){
  #fileName <- files[1]
  groupNo <- sub(".ReconciledTree","",fileName)
#  if(file.exists(paste0(charmatPath,"cluster.",groupNo,".dat"))){
  #read nhx tree as string, remove [] and read as newick formated string
  geneTreeString <- readLines(con = paste0(treePath, groupNo, ".ReconciledTree"))
  geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
  geneTree <- read.tree(text = geneTreeString)
  geneNames <- geneTree$tip.label
  geneTree$tip.label <- substr(geneNames, 0, 3)
  if(length(geneTree$tip.label) == length(spTree$tip.label)){
    if(identical(sort(geneTree$tip.label), sort(spTree$tip.label))){ # ignore trees of wrong size
      if(dist.topo(spTree, geneTree) == 0){
        selectedGroups <- append(selectedGroups,groupNo)
      }
    }
  }
  counter <- counter + 1
}

write.table(selectedGroups, file = "~/data/IES/analysis/tables/geneTreeSpeciesTree.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



