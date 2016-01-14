# calculate length of aligned blocks per cluster
path <- "~/data/IES_data/msas/alignments/filtered/"
files <- dir(path = path, pattern = "*.aln.fasta.gblocks")
l <- length(files)
blockLengths <- data.frame(cluster = character(l), nucleotides = numeric(l), stringsAsFactors = FALSE)
counter <- 1
for(fileN in files){
  # fileN <- files[1]
  print(paste0(counter, "/", l))
  cluster <- sub("^cluster\\.(\\d+)\\.aln\\.fasta\\.gblocks$", "\\1", fileN, perl = TRUE)
  pathFile <- paste0(path, fileN)
  if(file.size(pathFile) == 0){
    next # skip empty files
  }
  DF <- read.table(pathFile)
  blockLengths[counter, 1] <- cluster
  blockLengths[counter, 2] <- sum(DF[,3] - DF[2] + 1)
  counter <- counter + 1
}
save(blockLengths, file = "~/data/IES_data/rdb/blockLengths")