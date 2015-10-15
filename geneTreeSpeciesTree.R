# select gene trees with the same topology as the species tree
library(ape)
selectedGroups <- numeric(0)
treePath <- "~/data/IES_data/msas/phyldog/results/"
charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
files <- dir(treePath,pattern="[^.]*.ReconciledTree$")
spTreeString <- "(((PBIA:1,PTET:1):1,PSEX:1):1, PCAU:1);"
spTree <- read.tree(text=spTreeString)
spGT <- data.frame(stringsAsFactors = FALSE) # data frame with details of gene trees ~ species trees
counter <- 1
for(fileName in files){
  #fileName <- files[8]
  groupNo <- sub(".ReconciledTree","",fileName)
  if(file.exists(paste0(charmatPath,"cluster.",groupNo,".dat"))){
    #read nhx tree as string, remove [] and read as newick formated string
    geneTreeString <- readLines(con = paste0(treePath, groupNo, ".ReconciledTree"))
    geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
    geneTree <- read.tree(text = geneTreeString)
    geneNames <- geneTree$tip.label
    geneTree$tip.label <- substr(geneNames, 0, 4)
    if(length(geneTree$tip.label) == 4 && sort(geneTree$tip.label) == sort(spTree$tip.label)){ # ignore trees of wrong size
     # if(dist.topo(spTree, geneTree) == 0){
        selectedGroups <- append(selectedGroups,groupNo)
        spGT <- rbind(spGT,
                      data.frame(subtree = counter,
                                 cluster = as.integer(groupNo),
                                 geneId = geneNames,
                                 group = "ingroup",
                                 species = gene2species(geneNames),
                                 stringsAsFactors = FALSE
                                 ),
                      make.row.names = FALSE
                      )
     # }
    }
  }
  counter <- counter + 1
}
save(selectedGroups, spGT, file = "~/data/IES_data/rdb/geneTreeSpeciesTree")


