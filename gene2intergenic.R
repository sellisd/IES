# read exon bed files and print corresponding intron bed files
source("~/projects/IES/src/sharedFunctions.R")
args<-commandArgs(TRUE);
inFile <- args[1]
outFile <- args[2]
lengthsFile <- args[3]
#Rscript --vanilla gene2intergenic.R inFile outFile scaffoldLengthsFile
#"~/data/IES/analysis/pbi.scaf"
scaffoldLengths <- read.table(lengthsFile, stringsAsFactors = FALSE, header = TRUE)
genes <- read.table(inFile, stringsAsFactors = FALSE)
names(genes) <- c("scaffold", "start", "end", "name")
intergenic <- gene2intergenicBed(genes, scaffoldLengths)
write.table(intergenic, file = outFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)