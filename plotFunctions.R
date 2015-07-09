# useful functions for plotting
source("~/projects/fgmo/colors.R")

colBySpec <- function(tree){
  # color code tips of a tree based on species
  names <- tree$tip.label
  speciesAbr <- substr(names,0,4)
  colV <- speciesAbr
  colV[which(speciesAbr == "PBIA")] <- cgreen1
  colV[which(speciesAbr == "PTET")] <- cred1
  colV[which(speciesAbr == "PSEX")] <- cblue1
  colV[which(speciesAbr == "PCAU")] <- "grey60"
  colV[which(speciesAbr == "TTHE")] <- cbrown1
  names(colV) <- names
  colV
}

linkNodes <- function(tr,ktr){
  # read a phyldog outpout file and the corresponding key file
  # and create a dictionary linking the two node numbering schemes
  library(ape)
  library(phangorn)
  nodes <- unique(c(tr$edge))
  dict <- matrix(nrow=length(nodes),ncol=2)
  #find root
  #get Descendents create node correspondance
  i <- 1
  for(node in nodes){
    offspring <-Descendants(tr,node,type="tips")
    names <- tr$tip.label[offspring[[1]]]
    keyString <- paste(sort(names),collapse="")
    dict[i,] <- c(node,ktr[which(ktr[,2]==keyString),1])
    i <- i+1
  }
  dict # return matrix with node correspondance
}

sumAsr <- function(asr,burnIn,iesNo,dict){
  # calculate probabilities of presence/absence from treelist files (multibiphy output)
  # INTPUT: asr      file as table
  #         burnIn   cycles to skip
  #         iesNo    which IES to use
  #         dict     is the node correspondance
  ps <- matrix(ncol=2,nrow=length(dict[,1])) # P[as=1] for each node
  i <- 1
  for(node in dict[,2]){ # foreach node find ancestral state of ies 'iesNo'
    columnI <- which(names(asr)==paste("X",node,sep="")) # column with asr of node 'node'
    vector <- as.numeric(substr(asr[,columnI],iesNo,iesNo))
    p <- mean(vector[burnIn:length(vector)]) # probability of state 1
    ps[i,] <- c(1-p,p) # Prob absence, Prob presence
    i <- i+1
  }
  # nodes is the R ape numbering system
  # dict[,2] is the phyldog numebring system
  ps # return vector of probabilities
}
