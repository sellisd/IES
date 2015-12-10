# useful functions for plotting
source("~/projects/fgmo/colors.R")

# required libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(binom))

recttext <- function(xl, yb, xr, yt, text, rectArgs = NULL, textArgs = NULL) {
  # from http://stackoverflow.com/questions/31371296/how-to-write-text-inside-a-rectangle-in-r
  center <- c(mean(c(xl, xr)), mean(c(yb, yt)))
  do.call('rect', c(list(xleft = xl, ybottom = yb, xright = xr, ytop = yt), rectArgs))
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

plotCompl <- function(matBool, matBoolR, flankLength, speciesName){
  # make complementarity plot
  windowSize <- ncol(matBool) - flankLength
  ci <- binom.confint(colSums(matBool), rep(nrow(matBool), ncol(matBool)), methods = "exact")
  ciR <- binom.confint(colSums(matBoolR), rep(nrow(matBoolR), ncol(matBoolR)), methods = "exact")
  bpCoo <- barplot(ci$mean, 
                   ylim = c(0, 1.2), 
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
  colV[which(speciesAbr == "PBIA")] <- cgreen1
  colV[which(speciesAbr == "PTET")] <- cred1
  colV[which(speciesAbr == "PSEX")] <- cblue1
  colV[which(speciesAbr == "PCAU")] <- "grey60"
  colV[which(speciesAbr == "TTHE")] <- cbrown1
  names(colV) <- names
  colV
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