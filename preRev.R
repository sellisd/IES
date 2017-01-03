# prepare data for revBayes
library(reshape2)
library(ape)
setwd('.')
source("sharedFunctions.R")
opt <- loadUserOptions()
basePath <- opt["basePath", ]
for(run in c(1,2,3)){
  outdir <- file.path(basePath, 'analysis', 'asr', run)
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  charMats <- read.table(file.path(basePath, 'analysis', 'iesdb', 'charMats.tab'), stringsAsFactors = FALSE, header = TRUE)
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
    dataNexusOut <- file.path(basePath, 'analysis', paste0('asr', run), paste0("charMat",clusters[i],".nexus"))
    treeFileIn   <- file.path(basePath, 'analysis', paste0('phyldogT', run), "results", paste0(clusters[i],".ReconciledTree"))
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
      write.nexus(geneTree, file = file.path(basePath, 'analysis', paste0("asr", run), paste0("tree", clusters[i], ".nexus")))
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
  
  write.table(geneFamiliesProcessed, file = file.path(basePath, 'analysis', paste0('asr', run), "geneFamilies.dat"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  print("skipped families")
  print(skippedFamilies)
  write.table(data.frame(geneFamily = TgeneFamily, homIES = Tid, column = Tcolumn, stringsAsFactors = FALSE), 
              file = file.path(basePath, "analysis", "tables", paste0("homIES", run, ".columns.link")),
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}
