# classify homologous IES columns at different patterns of conservation
load("~/data/IES_data/msas/geneTreesSpeciesTrees")
charmatPath <- "~/data/IES_data/msas/alignments/charMat/"
cons <- data.frame(stringsAsFactors = FALSE)
for(group in selectedGroups){
  #group <- selectedGroups[1]
  a <- read.table(paste0(charmatPath,"cluster.",group,".dat"), header = TRUE, as.is = TRUE, row.names = 1)
  # remove rows with T. thermophila (that could be multiple)
  spNames <- substr(row.names(a),0,4)
  keepRows <- which(spNames == "PBIA" | spNames == "PTET" | spNames == "PSEX" | spNames == "PCAU")
  a[a != 0] <- 1
  a[a == 0] <- 0
  a <- a[keepRows, , drop = FALSE]
  row.names(a) <- spNames[keepRows]
  for(i in ncol(a)){
    #  cons <- rbind(cons, data.frame(group = group, col = i, inpbi = a["PBIA", i], inptet = a["PTET", i], inpsex = a["PSEX", i], inpca = a["PCAU", i], stringsAsFactors = FALSE))
    cons <- rbind(cons, data.frame(group = group, col = i, pattern = paste(a[c("PBIA", "PTET", "PSEX", "PCAU"), i], sep="", collapse = "")))
  }
}
save(cons, file = "~/data/IES_data/rdb/conservationPatterns")
