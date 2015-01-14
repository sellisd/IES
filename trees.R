args<-commandArgs(TRUE);
alnFile <- args[1]
group <- args[2]
output <- args[2]
# trees.R: estimate tree with maximum likelihood
# load libraries
library(ape)
library(phangorn)
library(methods)
# set input output
aln <- read.phyDat(alnFile,type="DNA", format="fasta")
outputF <- paste(output,sep="")

dM <- dist.dna(as.DNAbin(aln)) # calculate distance matrix
njT <- njs(dM) # calculate NJ tree
fit <- pml(njT,aln) # calculate likelihood of tree
fitJC <- optim.pml(fit,potNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optEdge=TRUE,control = pml.control(trace=0)) # optimize tree
spNames  <- substr(fitJC$tree$tip.label,0,4)
outgroup <- which( spNames != "PBIA" & spNames != "PTET" & spNames != "PSEX")
rerooted  <- tryCatch({root(fitJC$tre,fitJC$tree$tip.label[outgroup],resolve.root=TRUE)},
                      error = function(e){"lathos"})
# draw trees in file
png(paste("tree.",group,".png",sep=""))
if(class(rerooted) == "phylo"){
    is.rooted(rerooted)
    plot(rerooted,main=group) # plot tree
    write.tree(rerooted,file=outputF) # save tree
}else{
    plot(fitJC$tre,main=paste(group,"unrooted")) # plot tree
    write.tree(fitJC$tre,file=outputF) # save tree
}
dev.off()
