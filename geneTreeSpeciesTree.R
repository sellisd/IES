# select gene trees with the same topology as the species tree
library(ape)
selectedGroupsA <- numeric(0)
selectedGroupsB <- numeric(0)
selectedGroupsC <- numeric(0)

treePath <- "~/data/IES/analysis/phyldog/results/"
#charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(treePath,pattern="[^.]*.ReconciledTree$")

spTreeStringA <- "((PSE:1,PSO:1):1,(((PPE:1,PPR:1):1,(PBI:1,(POC:1,PTE:1):1):1):1,PTR:1):1):1;"
spTreeStringB <- "(PCA:1,((PSE:1,PSO:1):1,(((PPE:1,PPR:1):1,(PBI:1,(POC:1,PTE:1):1):1):1,PTR:1):1):1);"
spTreeStringC <- "((PCA:1,((PSE:1,PSO:1):1,(((PPE:1,PPR:1):1,(PBI:1,(POC:1,PTE:1):1):1):1,PTR:1):1):1):1,TTH:1);"


spTreeA <- read.tree(text=spTreeStringA)
spTreeB <- read.tree(text=spTreeStringB)
spTreeC <- read.tree(text=spTreeStringC)

ta <- read.tree("~/data/IES/analysis/brlen/brlenA/part.nexus.treefile")
tb <- read.tree("~/data/IES/analysis/brlen/brlenB/part.nexus.treefile")
tc <- read.tree("~/data/IES/analysis/brlen/brlenC/part.nexus.treefile")
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
  counter <- 1
  for(spTree in c(spTreeA, spTreeB, spTreeC)){
    if(length(geneTree$tip.label) == length(spTree$tip.label)){
      if(identical(sort(geneTree$tip.label), sort(spTree$tip.label))){ # ignore trees of wrong size
        if(dist.topo(spTree, geneTree) == 0){
          if(counter == 1) {
            selectedGroupsA <- append(selectedGroupsA, groupNo)
          }
          if(counter == 2){
            selectedGroupsB <- append(selectedGroupsB, groupNo)
          }
          if(counter == 3){
            selectedGroupsC <- append(selectedGroupsC, groupNo)
          }
        }
      }
    }
    counter <- counter + 1
  }
}

length(selectedGroupsA)
length(selectedGroupsB)
length(selectedGroupsC)

write.table(selectedGroupsA, file = "~/data/IES/analysis/tables/geneTreeSpeciesTreeA.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(selectedGroupsB, file = "~/data/IES/analysis/tables/geneTreeSpeciesTreeB.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(selectedGroupsC, file = "~/data/IES/analysis/tables/geneTreeSpeciesTreeC.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



