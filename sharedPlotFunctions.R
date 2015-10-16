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