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
  y <- c(Pbi, PbiPte, PbiPtePse, PbiPtePsePca)
  pbiBP <- barplot(100 * y/totalPbi, ylim = c(0, 100), names.arg = c("P. biaurelia", "P. biaurelia\nP. tetraurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. biaurelia IES (%)")))
  text(pbiBP, 10+100*y/totalPbi, paste(y, "\n(",round(100*y/totalPbi,2),"%)"))
  
  y <- c(Pte, PbiPte, PbiPtePse, PbiPtePsePca)
  pteBP <- barplot(100 * y/totalPte, ylim = c(0, 100), names.arg = c("P. tetraurelia", "P. biaurelia\nP. tetraurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. tetraurelia IES (%)")))
  text(pteBP, 10+100*y/totalPte, paste(y, "\n(",round(100*y/totalPte,2),"%)"))
  
  y <- c(Pse, PbiPtePse, PbiPtePsePca)
  pseBP <- barplot(100 * y/totalPse, ylim = c(0, 100), names.arg = c("P. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. sexaurelia IES (%)")))
  text(pseBP, 10+100*y/totalPse, paste(y, "\n(",round(100*y/totalPse,2),"%)"))
  
  y <- c(Pca, PbiPtePsePca)
  pcaBP <- barplot(100 * y/totalPca, ylim = c(0, 100), names.arg = c("P. caudatum", "P. biaurelia\nP. tetraurelia\nP. sexaurelia\nP. caudatum"), ylab = expression(italic("P. caudatum IES(%)")))
  text(pcaBP, 10+100*y/totalPca, paste(y, "\n(",round(100*y/totalPca,2),"%)"))
}