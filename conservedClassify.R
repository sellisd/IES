# classify homologous IES columns at different patterns of conservation
load("~/data/IES_data/msas/geneTreeSpeciesTree")
charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
cons <- data.frame(stringsAsFactors = FALSE)
for(group in selectedGroups){
  #group <- selectedGroups[1]
  iesNames <- read.table(paste0(charmatPath,"cluster.",group,".dat"), header = TRUE, as.is = TRUE, row.names = 1)
  # remove rows with T. thermophila (that could be multiple)
  spNames <- substr(row.names(iesNames),0,4)
  keepRows <- which(spNames == "PBIA" | spNames == "PTET" | spNames == "PSEX" | spNames == "PCAU")
  a <- iesNames # binary copy
  a[a != 0] <- 1
  a[a == 0] <- 0
  a <- a[keepRows, , drop = FALSE]
  row.names(a) <- spNames[keepRows]
  for(i in c(1:ncol(a))){
    #i <- 1
    cons <- rbind(cons, 
                  data.frame(group = group,
                             col = i,
                             pattern = paste(a[c("PBIA", "PTET", "PSEX", "PCAU"), i], sep="", collapse = ""),
                             ies = iesNames[a[,i] != 0, i],
                             stringsAsFactors = FALSE)
                  )    
  }
}
save(cons, file = "~/data/IES_data/rdb/conservationPatterns")
