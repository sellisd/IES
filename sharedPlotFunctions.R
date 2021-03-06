# useful functions for plotting
source("~/projects/fgmo/colors.R")

# required libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(binom))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(RColorBrewer))

abr = c( 'ppr' = 'P. primaurelia',
         'pbi' = 'P. biaurelia',
         'pte' = 'P. tetraurelia',  
         'ppe' = 'P. pentaurelia',
         'pse' = 'P. sexaurelia',
         'poc' = 'P. octaurelia',
         'ptr' = 'P. tredecaurelia',
         'pso' = 'P. sonneborni',
         'pca' = 'P. caudatum',
         'tth' = 'T. thermophila')

# from http://stackoverflow.com/questions/29943251/displaying-values-from-a-character-vector-as-italic-labels-in-boxplot-in-r
# by http://stackoverflow.com/users/516548/g-grothendieck
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

printIESalign <- function(geneFamily, iesColumn, pad, filtered, introns){
  # plot a multiple sequence alignment around an IES
  # print the nucleotide alignment around iesColumn with TAs annotated
  DF <- charMats[charMats$cluster == geneFamily & charMats$column == iesColumn, c("geneId", "ies", "begin", "end")]
  #printRegion(geneFamily, DF$begin[1] - pad, DF$end[1] + pad)
  alignF <- paste0("/home/dsellis/data/IES_data/msas/alignments/filtered/cluster.", geneFamily, ".nucl.fa")
  aln <- read.alignment(file = alignF, format = "fasta")
  alnM <- as.matrix.alignment(aln)
  # find in which genes IES are present
  geneIds <- DF$geneId[which(DF$ies != "0" & !is.na(DF$ies))]
  I <- which(row.names(alnM) %in% geneIds)
  IES <- DF[which(DF$ies != "0" & !is.na(DF$ies)), ]
  for(i in I){
    alnM[i, DF$begin[1]] <- toupper(alnM[i, DF$begin[1]])
    alnM[i, DF$end[1]] <- toupper(alnM[i, DF$end[1]])
  }
  #print(alnM[, c((DF$begin[1] - pad):(DF$end[1] + pad))])
  m <- alnM[, c((DF$begin[1] - pad):(DF$end[1] + pad))]
  plotMSA(m, IES, filtered, introns)
}

dnacol <- function(a){
  # a color pallette for DNA sequence printing
#  a <- c("a","c","t","g")
  p1 <- c("A" = dred,
          "T" = dblue,
          "C" = dgreen,
          "G" = dorange)
  p1[toupper(a)]
}


plotMSA <- function(m, IES, filtered, introns){
  # plot a MSA, shade IES insertion locations and mark arbitrary regions
  namesWidth <- max(sapply(row.names(m),nchar))
  y <- nrow(m)
  x <- as.numeric(dimnames(m)[[2]])
  plot.new()
  #leave some space between name and sequence
  sp <- 2
  xmin <- min(x) - namesWidth - sp
  plot.window(xlim = c(xmin, max(x)), ylim = c(-y, 0))
  greytone <- c("#c1c9c8", "#a5b5ab")
  for(i in x){
    rect(i - 0.5, -y - 1, i + 0.5, 0, border = NA, col = greytone[findTriplet(i)%%2+1])
  }
  for(yi in c(1:y)){
    if(row.names(m)[yi] %in% IES$geneId){ # if it has an IES draw background
      be <- IES[IES$geneId==row.names(m)[yi], c("begin","end")]
      points(be[1], y = - yi, col = "grey60",pch = 19, cex = 3)
      points(be[2], y = - yi, col = "grey60",pch = 19, cex = 3)
    }
    filtered <- filtered[filtered$eventType=="loss"] #only mark lost IES
    if(row.names(m)[yi] %in% filtered$gene){
      be <- filtered[which(filtered$gene==row.names(m)[yi]), c("begin","end")]
      points(be[1], y = -yi, col = "grey60", cex = 3)
      points(be[2], y = -yi, col = "grey60", cex = 3)
    }
    if(row.names(m)[yi] %in% introns$gene){
      be <- introns[which(introns$gene==row.names(m)[yi]), c("begin","end")]
      points(be[, 1], y = rep(-yi, nrow(be)), col = "black", cex = 3)
      points(be[, 2], y = rep(-yi, nrow(be)), col = "black", cex = 3)
    }
    text(x = x, y = - yi, m[yi, ], col = dnacol(m[yi, ]))
    geneName <- s2c(row.names(m)[yi])
    text(x = c(xmin:(xmin + length(geneName) - 1)) , y = - yi, geneName)
  }
  par(las = 3)
  axis(3, at = x, labels = x)
}

plotASRScenario <- function(cluster){
  #plot scenarios of ancestral state presence using the output of ancestralStates.R which summarizes the revBayes output into an tidyR format.
  cluster <- cluster
  name <- paste0("~/data/IES_data/rdb/ancestralStates/",cluster)
  if (!file.exists(name)){
    return("File does not exist")
  }
  load(name)
  nodePresence <- aggregate(ancestralStates$presenceAbsence, by = list(iesColumn = ancestralStates$iesColumn, nodeR = ancestralStates$r), FUN = sum)
  totalIterations <- length(unique(ancestralStates$Iteration))
  nodePresence <- data.frame(nodePresence, mean = nodePresence$x/totalIterations, stringsAsFactors = FALSE)
  ci <- binom.confint(nodePresence$x, rep(totalIterations, nrow(nodePresence)), methods = "exact", conf.level = 0.95)
  DF <- data.frame(ci[nodePresence$iesColumn == 2, c("lower", "upper")], nodePresence[nodePresence$iesColumn == 2, ])
  
  wide <- ancestralStates[ancestralStates$iesColumn == 2, c("Iteration", "presenceAbsence", "r")]
  wide <- spread(wide, r, presenceAbsence)
  key <- wide[, c(2:ncol(wide))]
  nodesRids <- names(key)
  key <- apply(key, 1, paste, collapse = "")
  simCounts <- table(key)
  total <- sum(simCounts)
  
  Pscenario <- numeric(0)
  for(i in c(1:length(simCounts))){
    pattern <- names(simCounts)[i]
    DF <- pattern2nodes(pattern, nodesRids)
    plotAsrInt(geneTree@phylo, title = paste("P:", round(simCounts[i]/total, 4)), DF)
    Pscenario <- append(Pscenario, simCounts[i]/total)
  }
  unname(Pscenario)
}

plotAsrInt <- function(phylo, title, DF){
  # simple ape-based plot of ancestral states with presence absence states only
  cl <- c(0,"grey60")
  plot(phylo, main = title, adj = 0.5)
  tipI <- which(DF$r <= length(phylo$tip.label))
  intI <- which(DF$r > length(phylo$tip.label))
  painternal <- DF$pa[intI]
  patips <- DF$pa[tipI]
  nodelabels(text = painternal, node = DF$r[intI], bg = cl[painternal + 1])
  tiplabels(text = patips, tip = DF$r[tipI], bg = cl[patips + 1])
}

plotAsrFloat <- function(phylo, title, DF){
  # simple ape-based plot of ancestral states with presence absence states only
  # str(DF) data.frame r: num, pa: num
  cl <- c(1, 0)
  plot(phylo, main = title, adj = 0.5)
  tipI <- which(DF$r <= length(phylo$tip.label))
  intI <- which(DF$r > length(phylo$tip.label))
  painternal <- DF$pa[intI]
  patips <- DF$pa[tipI]
  nodelabels(thermo = painternal, node = DF$r[intI], piecol = cl)
  tiplabels(thermo = patips, tip = DF$r[tipI], piecol = cl)
}

recttext <- function(xl, yb, xr, yt, text, rectArgs = NULL, textArgs = NULL) {
  # from http://stackoverflow.com/questions/31371296/how-to-write-text-inside-a-rectangle-in-r
  center <- c(mean(c(xl, xr)), mean(c(yb, yt)))
  do.call('rect', c(list(xleft = xl, ybottom = yb, xright = xr, ytop = yt), rectArgs))
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

plotCompl <- function(matBool, matBoolR, flankLength, speciesName, conf.level = 0.95){
  # make complementarity plot
  windowSize <- ncol(matBool) - flankLength
  ci <- binom.confint(colSums(matBool), rep(nrow(matBool), ncol(matBool)), methods = "exact", conf.level = conf.level)
  ciR <- binom.confint(colSums(matBoolR), rep(nrow(matBoolR), ncol(matBoolR)), methods = "exact", conf.level = conf.level)
  bpCoo <- barplot(ci$mean, 
                   ylim = c(0, 1), 
                   names.arg = c(-flankLength:-1, 1:windowSize),
                   col=c(rep("grey", flankLength),
                         dblue, dblue, 
                         rep("grey", (windowSize - 2))),
                   xlab = "",
                   ylab = "identity",
                   main = bquote(italic(.(speciesName))))
  arrows(bpCoo, ci$lower, bpCoo, ci$upper, code = 0, angle = 90)
  points(bpCoo, ciR$mean, type = "l", lwd = 2)
  points(bpCoo, ciR$lower, type = "l", lty = 3)
  points(bpCoo, ciR$upper, type = "l", lty = 3)
  prevParXPD <- par()$xpd
  par(xpd = NA)
  recttext(bpCoo[flankLength]+.5, 0, bpCoo[length(bpCoo)]+.5, -.05, "IES", rectArgs = list(col = "black"), textArgs = list(col = "grey90"))
  par(xpd = prevParXPD)
  
#   mb <- colMeans(matBool, na.rm = TRUE)
#   sb <- apply(matBool, 2, sd, na.rm = TRUE)
#   mbr <- colMeans(matBoolR, na.rm = TRUE)
#   sbr <- apply(matBoolR, 2, sd, na.rm = TRUE)
#   bpCoo <- barplot(mb, ylim = c(0, 1.2), names.arg = c(-flankLength:-1, 1:windowSize),col=c(rep("grey", flankLength),     dblue,dblue,rep("grey", (windowSize - 2))), xlab = "", ylab = "identity", main = bquote(italic(.(speciesName))))
#   arrows(bpCoo, mb + sb, bpCoo, mb - sb, code = 0, angle = 90)
#   points(bpCoo,mbr, type = "l", lwd = 2)
#   points(bpCoo,mbr + sbr, type = "l", lty = 3)
#   points(bpCoo,mbr - sbr, type = "l", lty = 3)
#   prevParXPD <- par()$xpd
#   par(xpd = NA)
#   recttext(bpCoo[flankLength]+.5, 0, bpCoo[length(bpCoo)]+.5, -.05, "IES", rectArgs = list(col = "black"), textArgs = list(col = "grey90"))
#   par(xpd = prevParXPD)
}

asrPlot <- function(ph, orderedMM){
  # function that plots the ancestral states on a ggtree
  # Input: ggree::read.nhx object and an orderedMM data.frame
  spEventsV <- extractEvents(ph)
  spN <- data.frame(tipLabels = ph@phylo$tip.label, speciesName = unname(gene2species(ph@phylo$tip.label)))
  data4tree <- cbind(orderedMM[orderedMM$iesColumn == 1, c("r", "mean", "sd")], spEventsV)
  iesColumns <- levels(orderedMM$iesColumn)
  p <- ggplot(ph@phylo, aes(x, y)) + geom_tree() + theme_tree() + xlab("") + ylab("")
  for(i in c(1:length(iesColumns))){
    print(p %<+% spN + geom_text(aes(color = speciesName, label=label, adj = -.1)) + scale_size(range = c(0,10)) +
            geom_point(aes(size = orderedMM$mean[orderedMM$iesColumn == iesColumns[i]])) +
            ggtitle(paste("gene family:", cluster, "IES:", iesColumns[i])))
  }
}

colBySpec <- function(tree){
  # color code tips of a tree based on species
  names <- tree$tip.label
  speciesAbr <- substr(names,0,4)
  colV <- speciesAbr
  cls <- brewer.pal(10, "Set3")
  colV[which(speciesAbr == "PPRI")] <- cls[1]
  colV[which(speciesAbr == "PBIA")] <- cls[3]
  colV[which(speciesAbr == "PTET")] <- cls[4]
  colV[which(speciesAbr == "PPEN")] <- cls[5]
  colV[which(speciesAbr == "PSEX")] <- cls[6]
  colV[which(speciesAbr == "POCT")] <- cls[7]
  colV[which(speciesAbr == "PTRE")] <- cls[8]
  colV[which(speciesAbr == "PSON")] <- cls[9]
  colV[which(speciesAbr == "PCAU")] <- cls[10]
  colV[which(speciesAbr == "TTHE")] <- "black"
  names(colV) <- names
  colV
}

shortNames <- function(tipNames){
  # abbreviate species names, for plotting with larger fonts
  shortNames <- sub("Paramecium", "P", tipNames, fixed = TRUE)
  shortNames <- sub("Tetrahymena", "T", shortNames, fixed = TRUE)
  shortNames
}


plotConsPast <- function(pattern){
  # function to plot conservation barplots
  # INPUT: a matrix with patterns of conservation for each species
  #  skip if at least one ingroup genes are missing, or no outgroup
  skipBV <- c(pattern[, 1] == -1 | pattern[, 2] == -1 | pattern[, 3] == -1 | pattern[, 4] == -1)
  pattern <- pattern[!skipBV, ]
  stringPat <- apply(pattern, 1, paste, collapse = "")
  inPbi <- pattern[, 1] == 1
  inPte <- pattern[, 2] == 1
  inPse <- pattern[, 3] == 1
  inPca <- pattern[, 4] == 1
  
  Pbi <- sum( inPbi & !inPte & !inPse & !inPca) # Number of IES only in P. biaurelia
  Pte <- sum(!inPbi &  inPte & !inPse & !inPca) # Number of IES only in P. tetraurelia
  Pse <- sum(!inPbi & !inPte &  inPse & !inPca) # Number of IES only in P. sexaurelia
  Pca <- sum(!inPbi & !inPte & !inPse &  inPca) # Number of IES only in P. caudatum
  PbiPte <- sum( inPbi &  inPte &  !inPse & !inPca)
  PbiPtePse <- sum( inPbi &  inPte &  inPse & !inPca)
  PbiPtePsePca <- sum( inPbi &  inPte &  inPse & inPca)
  
  totalPbi <- sum(inPbi) # Number of IES in P. biaurelia
  totalPte <- sum(inPte) # Number of IES in P. biaurelia
  totalPse <- sum(inPse) # Number of IES in P. biaurelia
  totalPca <- sum(inPca) # Number of homologous IES with an element present in P. caudatum
  
  par(las = 2, bty = "l", oma = c(1, 0, 0, 0))
  ybi <- c(Pbi, PbiPte, PbiPtePse, PbiPtePsePca)
  pbiBP <- barplot(100 * ybi/totalPbi, ylim = c(0, 100), names.arg = c("P. biaurelia", "P. biaurelia\nP. tetraurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. biaurelia IES (%)")))
  rybi <- round(100*ybi/totalPbi,2)
  text(pbiBP, 10+100*ybi/totalPbi, paste(ybi, "\n(", rybi,"%)"))
  
  yte <- c(Pte, PbiPte, PbiPtePse, PbiPtePsePca)
  pteBP <- barplot(100 * yte/totalPte, ylim = c(0, 100), names.arg = c("P. tetraurelia", "P. biaurelia\nP. tetraurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. tetraurelia IES (%)")))
  ryte <- round(100*yte/totalPte,2)
  text(pteBP, 10+100*yte/totalPte, paste(yte, "\n(", ryte,"%)"))
  
  yse <- c(Pse, PbiPtePse, PbiPtePsePca)
  pseBP <- barplot(100 * yse/totalPse, ylim = c(0, 100), names.arg = c("P. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. sexaurelia IES (%)")))
  ryse <- round(100*yse/totalPse,2)
  text(pseBP, 10+100*yse/totalPse, paste(yse, "\n(", ryse,"%)"))
  
  yca <- c(Pca, PbiPtePsePca)
  pcaBP <- barplot(100 * yca/totalPca, ylim = c(0, 100), names.arg = c("P. caudatum", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. caudatum IES(%)")))
  ryca <- round(100*yca/totalPca,2)
  text(pcaBP, 10+100*yca/totalPca, paste(yca, "\n(", ryca,"%)"))
  
  numbers <- c("Pbi" = Pbi, "Pte" = Pte, "Pse" = Pse, "Pca" = Pca,
               "PbiPte" = PbiPte, "PbiPtePse" = PbiPtePse, "PbiPtePsePca" = PbiPtePsePca,
               "totalPbi" = totalPbi, "totalPte" =  totalPte, "totalPse" = totalPse, "totalPca" = totalPca,
               "rybi" = rybi, "ryte" = ryte, "ryse" = ryse, "ryca" = ryca)
}

