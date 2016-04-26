#source("/home/dsellis/projects/fgmo/colors.R")
library(seqinr)
library(ape)
#colorCodes<-c("PCAU"=cbrown1,"PBIA"=cred1,"PTET"=cgreen1,"PSEX"=cblue1,"TTHE"= "black")
#pathT <- "/home/dsellis/data/IES_data/msas/alignments/dnd/"
pathA <- "/home/dsellis/data/IES/analysis/msas/filtered/"
filesF <- dir(path=pathA,pattern="*.nucl.fa")
m <- data.frame(matrix(nrow=length(filesF),ncol=4))
#pdf("/home/dsellis/projects/IES/allTrees.pdf")
#par(mfcol=c(2,2))
counter <- 1
for(f in filesF){
  print(paste(counter,length(filesF)))
  #alnF <- gsub(".dnd",".aln",f)
  b <- seqinr::read.alignment(paste(pathA,f,sep=""),format="fasta")
  cns <- seqinr::consensus(b,method="threshold",threshold=1)
  seqL <- length(cns) # length of consensus
  seqN <- b$nb #number of sequences
  ident <- 100*(1-sum(is.na(cns))/length(cns))
  if(sum(is.na(cns)) == length(cns)){
    stop("identical")
  }
  m[counter,] <- c(f, seqN, seqL, ident)
  counter <- counter + 1
  #a<-read.tree(paste(pathT,f,sep=""))  
  #if(length(a$tip.label)>10){
  #   plot.phylo(a,tip.color=colorCodes[substr(a$tip.label,0,4)],cex=0.5)
  #   title(f)
  # }
}
#dev.off()
#write.table(m,file="/home/dsellis/projects/IES/reports/pipeline/msaStats.dat",row.names=FALSE,col.names=FALSE,quote=FALSE)

