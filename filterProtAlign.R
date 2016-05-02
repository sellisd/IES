library(seqinr)
path <- "/home/dsellis/data/IES/analysis/msas/aln/"
filteredPath <- "/home/dsellis/data/IES/analysis/msas/filtered/"
if(!dir.exists(filteredPath)){
  dir.create(filteredPath)
}
files <- dir(path= path,pattern="*.aln.fa$")
nb <- numeric()
id <- numeric()
for(fileName in files){
  print(fileName)
  a <- read.alignment(paste(path,fileName,sep=""),format="fasta")
  nb <- append(nb,a$nb) #number of sequeces
  distM <- dist.alignment(a,matrix="identity")
  id <- append(id,mean(1-distM*distM)) #average pairwise identity
}
write.table(data.frame(file=files,sequenceNo=nb,avPairId=id),file="~/data/IES/analysis/msas/protAlignStats.dat", quote=FALSE, row.names = FALSE)
highId <- files[which(id>0.5)]
enoughSeq <- files[which(nb>3)]
filtered <- intersect(enoughSeq,highId)
from <- paste(path,filtered,sep="")
to <- paste(filteredPath,filtered,sep="")
r <- file.copy(from, to)
print(all(r)==TRUE)
