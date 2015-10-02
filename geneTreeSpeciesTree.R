# select gene trees with the same topology as the species tree
library(ape)
selectedGroups <- numeric(0)
treePath <- "~/data/IES_data/msas/phyldog/results/"
charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(path,pattern="[^.]*.ReconciledTree$")
spTreeString <- "(((PBIA:1,PTET:1):1,PSEX:1):1, PCAU:1);"
spTree <- read.tree(text=spTreeString)
for(fileName in files){
  #fileName <- files[1]
  groupNo <- sub(".ReconciledTree","",fileName)
  if(file.exists(paste0(charmatPath,"cluster.",groupNo,".dat"))){
    #read nhx tree as string, remove [] and read as newick formated string
    geneTreeString <- readLines(con = paste0(treePath, groupNo, ".ReconciledTree"))
    geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
    geneTree <- read.tree(text = geneTreeString)
    geneTree$tip.label <- substr(geneTree$tip.label,0,4)
    if(length(geneTree$tip.label) == 4 && sort(geneTree$tip.label) == sort(spTree$tip.label)){ # ignore trees of wrong size
      if(dist.topo(spTree, geneTree) == 0){
        selectedGroups <- append(selectedGroups,groupNo)
      }
    }
  }
}
save(selectedGroups, file = "~/data/IES_data/msas/geneTreesSpeciesTrees")

