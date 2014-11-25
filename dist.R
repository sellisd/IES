library("seqinr")
library("ape")
# read fasta file of largest protein silix cluster
a<-read.alignment(file="/Users//diamantis/data//IES_data/working//102Pclusters/cluster.1.fa",format="fasta")
# calculate the matrix of pairwise identity
m<-dist.alignment(a,matrix="identity")
#tr<-nj(m)
#colorCodes<-c("Pbi"=2,"Pte"=4,"Pse"=6)
#plot.phylo(tr,tip.color=colorCodes[substr(tr$tip.label,0,3)],cex=0.1)
#write.tree(tr,file="/Users/diamantis/data/IES_data//working/sp102.tre")
#geneNames<-sub("P.._scaffold.*T(\\d+).*","\\1",tr$tip.label,perl=TRUE)
#cr<-rainbow(length(geneNames))
#plot.phylo(tr,tip.color=cr,cex=0.1)
#subet distance matrix to calculate hist and average pairwise identities
# calulate the identity of 102pb exons within weach species
# and within each gene
l<-length(labels(m))
mm <- as.matrix(m)
#make vector of species
#make vector of gene numbers
speciesCodes <- c("Pbi" = 2, "Pte" = 4, "Pse" = 6)
speciesV <- speciesCodes[substr(labels(m),0,3)]
genesV <- sub("(P..)_scaffold.*T(\\d+).*","\\1\\2",labels(m),perl=TRUE)
PbiL <- length(which(speciesV == 2))
PteL <- length(which(speciesV == 4))
PseL <- length(which(speciesV == 6))
PbiPairs <- numeric(PbiL)
PtePairs <- numeric(PteL)
PsePairs <- numeric(PseL)
geneNames <- unique(genesV) # unique gene Names
geneDistSum <- numeric(length(geneNames)) # array of sum of pairwise distances 
geneDistCount <- numeric(length(geneNames)) # array of counts of pairwise distances calculated
PbiIndex <- 1
PteIndex <- 1
PseIndex <- 1
for(i in c(1:l)){
  print(i)
  for(j in c(i:l)){
#  print(paste(i,j,speciesV[i],speciesV[j], genesV[i], genesV[j]))
    if(speciesV[i] == 2 & speciesV[j] == 2){
      PbiPairs[PbiIndex] <- mm[i,j]
      PbiIndex <- PbiIndex + 1
    }
    if(speciesV[i] == 4 & speciesV[j] == 4){
      PtePairs[PteIndex] <- mm[i,j]
      PteIndex <- PteIndex + 1
    }
    if(speciesV[i] == 6 & speciesV[j] == 6){
      PsePairs[PseIndex] <- mm[i,j]
      PseIndex <- PseIndex + 1
    }
    if(genesV[i] == genesV[j]){
      counter <- 1
      for(n in geneNames){
          if(genesV[i] == n){
            geneDistSum[counter] <- geneDistSum[counter] + mm[i,j]
            geneDistCount[counter] <- geneDistCount[counter] + 1          
            }
          counter <- counter + 1
      }
    }
      #loop through unique names
    #if speciesV[i] == uniqueName and speciesV[j] also push in list genes[[name][index]
  }
}

write.table(PbiPairs,file="/Users/diamantis/data/IES_data/working/Pbi.dist",row.names=F,col.names=F)
write.table(PtePairs,file="/Users/diamantis/data/IES_data/working/Pte.dist",row.names=F,col.names=F)
write.table(PsePairs,file="/Users/diamantis/data/IES_data/working/Pse.dist",row.names=F,col.names=F)
write.table(geneDistSum/geneDistCount,file="/Users/diamantis/data/IES_data/working/gene.dist",row.names=F,col.names=F)