# prepare data for revBayes
library(reshape2)
library(ape)
source("~/projects/IES/src/sharedFunctions.R")

for(run in c(1,2,3)){
  charMats <- read.table("~/data/IES/analysis/iesdb/charMats.tab", stringsAsFactors = FALSE, header = TRUE)
  clusters <- unique(charMats$cluster)
  geneFamiliesProcessed <- character()
  skippedFamilies <- character()
  # build arrays with homologous IES group and column numbering correspondance
  TgeneFamily <- numeric(0)
  Tid <- numeric(0)
  Tcolumn <- numeric(0)
  
  for(i in c(1:length(clusters))){
    cat("run ", run, ":", i, "/", length(clusters), "\r")
    l <- save2nexus(clusters[i])
    nex <- l[[1]]
    spNames <- gene2species(names(nex))
    #check if we have a tree file
    dataNexusOut <- paste0("~/data/IES/analysis/asr", run, "/charMat",clusters[i],".nexus")
    treeFileIn <- paste0("~/data/IES/analysis/phyldogT", run, "/results/",clusters[i],".ReconciledTree")
    if(file.exists(treeFileIn)){
      # write a nexus data file for each gene family
      modif.write.nexus.data(nex, file = dataNexusOut, format = "Standard")
      # and strip the reconciled tree from NHX annotations
      geneTreeString <- readLines(con = treeFileIn) #read nhx tree as string, remove [] and read as newick formated string
      geneTreeString <- gsub("\\[&&NHX:Ev=[SD]:S=\\d+:ND=\\d+\\]","",geneTreeString, perl = TRUE)
      geneTree <- read.tree(text = geneTreeString)
      if (!identical(sort(geneTree$tip.label), sort(gene2protName(names(nex))))){
        stop(paste("Genes do not match between tree and character matrix!: ", clusters[i], i))
      }
      write.nexus(geneTree, file = paste0("~/data/IES/analysis/asr", run, "/tree",clusters[i],".nexus"))
      geneFamiliesProcessed <- append(geneFamiliesProcessed, clusters[i])
      # gather information linking homologous IES group id and order of columns in character matrix
      homCol <- l[[2]]
      columnsIGF  <- length(homCol) # homologous IES columns in gene family
      TgeneFamily <- append(TgeneFamily, rep(clusters[i], columnsIGF))
      Tcolumn     <- append(Tcolumn, c(0:(columnsIGF-1))) # 0-based
      Tid         <- append(Tid, homCol)
    }else{
      skippedFamilies <- append(skippedFamilies, clusters[i])
    }
  }
  
  write.table(geneFamiliesProcessed, file = paste0("~/data/IES/analysis/asr", run, "/geneFamilies.dat"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  print("skipped families")
  print(skippedFamilies)
  write.table(data.frame(geneFamily = TgeneFamily, homIES = Tid, column = Tcolumn, stringsAsFactors = FALSE), 
              file = paste0("~/data/IES/analysis/tables/homIES", run, ".columns.link"),
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}
