#iesInGenes
#read IESin file and for each line find transcript coordinates
#then print per gene
# read exon bed files and print corresponding intron bed files
source("~/projects/IES/src/sharedFunctions.R")
args<-commandArgs(TRUE);
inFile <- args[1]
cdsdbF <- args[2]
outFile <- args[3]
#Rscript --vanilla iesInGenes.R inFile cdsdbF
# inFile <- "~/data/IES/analysis/bed/ppr.IESin.be"
# cdsdbF <- "/home/dsellis/data/IES/analysis/iesdb/ppr.cdsdb"

iesin <- read.table(inFile, stringsAsFactors = FALSE)
cdsdb <- read.table(cdsdbF, stringsAsFactors = FALSE, header = TRUE)


iesin <- iesin[which(iesin$V6 == "cds"), ]
matchIndex <- match(iesin$V10, cdsdb$id)
genes <- cdsdb$geneName[matchIndex]
strands <- cdsdb$strand[matchIndex]
outFileCon <- file(outFile, open = "at")

for(i in unique(genes)){
  index <- which(genes == i)
  cds <- data.frame(name = iesin[index, 10], start = iesin[index, 8], end = iesin[index, 9], stringsAsFactors = FALSE)
  strand <- strands[index][1]
  for(loc in index){
    inCDS <- iesin[loc, 10]
    qrange <- c(iesin[loc, 2], iesin[loc, 3])
    inT <- inTranscript(cds, strand, qrange, inCDS)
    cat(file = outFileCon, i,inT[1], inT[2], iesin[loc,4], "\n", sep = "\t", append = TRUE)
  }
}


#write.table(intergenic, file = outFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)