# create a tidy data.frame with all character matrices
charMatPath <- "~/data/IES_data/msas/alignments/charMat/"
fileNames <- dir(path = charMatPath, pattern = "cluster.\\d+.dat")
clusters <- as.numeric(substr(fileNames, 9, (nchar(fileNames) - 4)))
charMats <- data.frame()
for(cluster in clusters){
  print(cluster)
  charMat <- read.table(paste(charMatPath, "cluster.", cluster, ".dat", sep=""), header = TRUE, as.is = TRUE)
  melted <- melt(charMat, id.vars = "geneName", variable.name = "begin", value.names = "IES", factorsAsStrings = TRUE)
  columnNames <- match(melted$begin, names(charMat)[-1])
  #remove X
  noX <- substr(as.character(melted$begin), 2, nchar(as.character(melted$begin)))
  # split start and end coordinates
  l <- strsplit(noX, ".", fixed = TRUE)
  # make character to numeric
  m <- matrix(as.numeric(unlist(l)), ncol = 2, byrow = TRUE)
  # add back to data.frame
  charMats <- rbind(charMats, data.frame(cluster = cluster, column = columnNames, geneName = melted$geneName, begin = m[,1], end = m[, 2], ies = melted$value))
}

save(charMats, file = "~/data/IES_data/rdb/charMats")