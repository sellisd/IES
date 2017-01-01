# Shared functions

# required libraries
suppressPackageStartupMessages(library(seqinr)) # make s2c availabel
suppressPackageStartupMessages(library(seqLogo))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(dplyr))

# shared variables

prefixes = c( 'PPRIM.AZ9-3.1.' = 'Paramecium primaurelia',
              'PBIA.V1_4.1.'   = 'Paramecium biaurelia',
              'PTET.51.1.'     = 'Paramecium tetraurelia',  
              'PPENT.87.1.'    = 'Paramecium pentaurelia',
              'PSEX.AZ8_4.1.'  = 'Paramecium sexaurelia',
              'POCTA.138.1.'   = 'Paramecium octaurelia',
              'PTRED.209.2.'   = 'Paramecium tredecaurelia',
              'PSON.ATCC_30995.1.'     = 'Paramecium sonneborni',
              'PCAU.43c3d.1.'  = 'Paramecium caudatum')

loadUserOptions <- function(){
  # load user options from python options file
  # example usage:
  # opt <- loadUserOptions()
  # opt["basePath",]
  setwd("~/projects/IES/src")
  opt <- read.table("./pyies/userOptions.py", stringsAsFactors = FALSE, sep = "=", strip.white = TRUE, header = FALSE, row.names = 1)
}

withinRange <- function(q, r){
  # if a be range is partially within another be range
  if((q[1]>=r[1]) & (q[1]<r[2])      # if query beginning is within range
     | ((q[2]>r[1]) & (q[2]<=r[2]))){ # or query end is within range
    return(1)
  }else{
    return(0)
  }
}

inTranscript <- function(cds, strand, qrange, inCDS){
  #data frame should be sorted
  if(!all(cds$start == sort(cds$start))){
    stop("unsorted input")
  }
  cumLength <- 0
  found <- 0
  if(strand == 1){
    for(i in c(1:nrow(cds))){
      cdsL <- cds$end[i] - cds$start[i]
      if(inCDS == cds$name[i]){
        # check if at least partially within
        if(!(withinRange(qrange, c(cds$start[i], cds$end[i])))){  #  or end is within CDS
          stop("IES is not in ", inCDS)
        }
        # partial sum
        startT <- qrange[1] - cds$start[i]
        endT <- qrange[2] - cds$start[i]
        if(startT < 0){ # if partial match
          startT <- 0
        }
        if(endT > cdsL){
          endT <- cdsL
        }
        startT <- startT + cumLength
        endT <- endT + cumLength
        found <- 1
        break
      }else{
        cumLength <- cumLength + cdsL
      } 
    }
  }else if(strand == -1){
    for(i in c(nrow(cds):1)){ # count backwards
      #cat(cumLength,"\n")
      cdsL <- cds$end[i] - cds$start[i]
      if(inCDS == cds$name[i]){
        # check if at least partially within
        if(!(withinRange(qrange, c(cds$start[i], cds$end[i])))){  #  or end is within CDS
          stop(stop("IES is not in ", inCDS))
        }
        startT <- cds$end[i] - qrange[2]
        endT <- cds$end[i] - qrange[1]
        if(startT < 0){
          startT <- 0
        }
        if(endT > cdsL){
          endT <- cdsL
        }
       # cat("   ",startT, endT,"\n")
        startT <- startT + cumLength
        endT <- endT + cumLength
        found <- 1
        break
      }else{
       cumLength <- cumLength + cdsL 
      }
    }
  }else{
    stop(paste("unknown strand", strand))
  }
  if(found == 0){
    stop("query not in cds")
  }
  c(startT, endT)
}

exon2intronsBed <- function(exons){
  #  read a data.frame from a bed file and return the complement keeping fixed a group
  # INPUT: exons
  # scaffold  (start   end]   name  gene
  # sc1         2       5     exon1 gene1
  # sc1         7      10     exon2 gene1
  # sc1         12     20     exon3 gene1
  # OUTPUT: introns

  genes <- unique(exons$gene)
  l <- nrow(exons)
  introns <- data.frame(scaffold = character(l), intronStart = numeric(l), intronEnd = numeric(l), gene = character(l), stringsAsFactors = FALSE)
  rowCounter <- 1
  for(gene in genes){
    # gene <- genes[1]
    geneI <- which(exons$gene == gene)
    if(length(geneI) < 2){
      next # no introns
    }
    intronStart <- exons[geneI[-length(geneI)], 3]
    intronEnd <- exons[geneI[-1], 2]
    if(!all(intronEnd - intronStart > 1)){
      stop("zero sized intron???")
    }
    for(i in c(1:length(intronStart))){
      introns[rowCounter, "scaffold"] <- exons[geneI[1], 1]
      introns[rowCounter, "intronStart"]    <- intronStart[i]
      introns[rowCounter, "intronEnd"]      <- intronEnd[i]
      introns[rowCounter, "gene"]     <- gene
      rowCounter <- rowCounter + 1
    }
    cat(rowCounter,"/",l,"\r")
  }
  introns[1:(rowCounter-1),]
}

gene2intergenicBed <- function(genes, scaffoldLengths){
  #  read a data.frame from a bed file and return the complement keeping fixed a groups
  scaffolds <- unique(genes$scaffold)
  l <- nrow(genes)
  intergenic <- data.frame(scaffold = character(l), interStart = numeric(l), interEnd = numeric(l), name = character(l), stringsAsFactors = FALSE)
  rowCounter <- 1
  for(scaffold in scaffolds){
    nameId <- 1
    # scaffold <- scaffolds[1]
    scafI <- which(genes$scaffold == scaffold)
    firstStart <- 0
    lastEnd <- scaffoldLengths$length[which(scaffoldLengths$scaffold == scaffold)]
    if(length(scafI) > 0){ # if at least one gene in scaffold
      firstEnd <- genes[scafI[1], 2]
      lastStart <- genes[scafI[length(scafI)], 3]
      if(length(scafI) > 1){ # there is at least one intergenic region between genes
        interStart <- c(firstStart, genes[scafI[-length(scafI)], 3], lastStart)
        interEnd <- c(firstEnd, genes[scafI[-1], 2], lastEnd)
      }else{
        interStart <- c(firstStart, lastStart)
        interEnd <- c(firstEnd, lastEnd)
      }
    }else{ # if no genes
      interStart <- 0
      interEnd <- lastEnd
    }
    for(i in c(1:length(interStart))){
      if(interStart[i] >= interEnd[i]){
        # Exclude cases where genes are overlapping or have no intergenic space between them
      }else{
        intergenic[rowCounter, "scaffold"]   <- scaffold
        intergenic[rowCounter, "interStart"] <- interStart[i]
        intergenic[rowCounter, "interEnd"]   <- interEnd[i]
        intergenic[rowCounter, "name"]       <- paste0(scaffold, '.', nameId)
        nameId <- nameId + 1
        rowCounter <- rowCounter + 1
        cat(rowCounter,"/",l,"\r")
      }
    }
  }
  intergenic[1:(rowCounter-1),]
}

table2dataFrame <- function(l){
  # convert a table to an x,y data.frame
  h <- table(l)
  x <- as.numeric(dimnames(h)[[1]])
  y <- as.vector(h)
  data.frame(x = x, y = y, stringsAsFactors = FALSE)
}

gene2protName <- function(geneNames){
  protNames <- character(length(geneNames))
  for(prefix in names(prefixes)){
    index <- which(substr(geneNames, 0, nchar(prefix)) == prefix) 
    protNames[index] <- paste0(prefix, "P", substr(geneNames[index], nchar(prefix) + 2, nchar(geneNames[index])))
  }
  tindex <- which(substr(geneNames, 0, 7) == "TTHERM_")
  if(length(tindex) !=0 ){
    protNames[tindex] <- paste0("TTHERM_", substr(geneNames[tindex], 8, nchar(geneNames[tindex])))
  }
  if(all(protNames != "")){
    return(protNames)
  }else{
    stop("error matching gene/protein names: ", geneNames[which(protNames == "")])
  }
}

gene2protCDS <- function(cds){
  # function that translates genomic to protein coordinates for CDS  
  # input is a data.frame with cdsId geneId geneStart, geneEnd and strand columns
  l <- nrow(cds)
  protcds <- data.frame(cdsId = character(l), geneId = character(l), geneStart = character(l), geneEnd = character(l), protStart = character(l), protEnd = character(l), stringsAsFactors = FALSE)
  geneCounter <- 1
  uniqueGenes <- unique(cds$geneId)
  #gene <- uniqueGenes[1]
  for(gene in uniqueGenes){
    cat(paste(geneCounter, "/", length(uniqueGenes),"\r"))
    rowIndex <- which(cds$geneId==gene)
    DF <- cds[rowIndex, ]
    # make sure all on the same strand
    if(all(DF$strand == "-1")){
      strand <-  -1
    }else if(all(DF$strand == "1")){
      strand <-  1
    }else{
      stop(paste("mixed strands within the same gene, or unknown strand: ", DF$strand))
    }
    # order by gene start
    #DF <- DF[order(DF$geneStart), ]
    l <- nrow(DF)
    CDSLengths <- DF$geneEnd - DF$geneStart + 1
    if(CDSLengths[l] == 3 & strand == 1){
        CDSLengths <- CDSLengths[-l]
        DF <- DF[-l, ]
        rowIndex <- rowIndex[-l]
        l <- l - 1
    }else if(CDSLengths[1] == 3 & strand == -1){
        CDSLengths <- CDSLengths[-1]
        DF <- DF[-1, ]
        rowIndex <- rowIndex[-1]
        l <- l - 1
    }
    if(l == 1){
        pStart <- 1
        pEnd <- CDSLengths[1] - 3 # remove termination codon
    }else{
      if(strand == 1){
        CDSLengths[l] <- CDSLengths[l] - 3 # remove termination codon
        cumulLen <- cumsum(CDSLengths)
        pEnd <- cumulLen
        pStart <- c(1, cumulLen[1:(l-1)]+1)
      }else if(strand == -1){
        CDSLengths[1] <- CDSLengths[1] - 3 # remove termination codon
        cumulLen <- cumsum(rev(CDSLengths))
        pStart <- rev(c(1, cumulLen[1:(l-1)]+1))
        pEnd <- rev(cumulLen)
      }
    }
    # remove the last codon (TGA)
    protcds[rowIndex, "cdsId"]  <- DF$cdsId
    protcds[rowIndex, "geneId"] <- DF$geneId
    protcds[rowIndex, "geneStart"] <- DF$geneStart
    protcds[rowIndex, "geneEnd"] <- DF$geneEnd
    protcds[rowIndex, "protStart"] <- pStart
    protcds[rowIndex, "protEnd"] <- pEnd
    geneCounter <- geneCounter + 1
  }
  protcds
  # remove empty rows, if one column is empty all should be
  protcds <- protcds[which(protcds$protEnd != ""), ]
}

introns <- function(cds){
  # from a dataframe of geneId and CDS Ids extract some information for introns
  l <- nrow(cds)
  intronsPerGene <- data.frame(geneId = character(l), introns = numeric(l), stringsAsFactors = FALSE)
  intronLengths <- numeric()
  geneCounter <- 1
  uniqueGenes <- unique(cds$geneId)
  for(gene in uniqueGenes){
    cat(paste(geneCounter, "/", length(uniqueGenes)),"\r")
    rowIndex <- which(cds$geneId==gene)
    DF <- cds[rowIndex, ]
    # order by gene start
    DF <- DF[order(DF$geneStart), ]
    l <- nrow(DF)
    if(l == 1){
      # no introns
      intronLength <- numeric(0)
    }else{
      intronLength <- DF$geneStart[2:l] - DF$geneEnd[1:(l-1)] - 1 # -1 because (geneStart, geneEnd)
      # cInLen <- c(1, cumsum(intronLength)) # cumulative length of introns
    }
    intronsPerGene[geneCounter, "geneId"] <- gene
    intronsPerGene[geneCounter, "introns"] <- length(intronLength)
    intronLengths <- append(intronLengths, intronLength)
    geneCounter <- geneCounter + 1
  }
  # remove extra rows
  intronsPerGene <- intronsPerGene[which(intronsPerGene$geneId != ""), ]
  list(intronsPerGene, intronLengths)
}

geneInScaffold <- function(geneNames){
  # given a gene find in which scaffold it is (old gene names)
  scaffolds <- character(length(geneNames))
  counter <- 1
  for(x in geneNames){
    if(substr(x,0,4) == "PTET"){
      scaffold <- paste0("scaffold51_", as.numeric(substr(x,12,14)))
    }else if (substr(x,0,4) == "PBIA"){
      scaffold <- paste0("scaffold_", substr(x,14,17)) # pbiaurelia scaffolds retain padding 0
    }else{
      warning(paste("unknown gene name: ", x))
    }
    scaffolds[counter] <- scaffold
    counter <- counter + 1
  }
  scaffolds
}

findFrame <- function(x){
  # find the frame of a nucleotide from its coordinates
  if(x<1){
    stop("coordinates should be positive")
  }
  (x-1)%%3+1
}

findTriplet <- function(x){
  # find the aa triplet of a nucleotide from its coordinates
  if(x<1){
    stop("coordinates should be positive")
  }
  trunc(x/3)
}

gainLossOnPath <- function(L){
  # from a list of probability of presence in ancestral nodes along a path calculate sum of insertion and loss probability
  l <- length(L)
  if(!all(L >= 0) & all(L <= 1)){
    warning(L)
  }
  gain <- L[2:l] - L[1:(l-1)]
  gain[gain<0] <- 0
  loss <- L[1:(l-1)] - L[2:l]
  loss[loss<0] <- 0
  c(sum(gain), sum(loss))
}

pattern2nodes <- function(pattern, nodesRids){
  # from a string pattern of presence absence generate a data.frame with node ids and presence absence
  data.frame(pa = as.numeric(s2c(pattern)), r = as.numeric(nodesRids), stringsAsFactors = FALSE)
}

eventType <- function(a, b){
  #classify events as gains, loss, presence or absence from two vectors of states
  l <- length(a)
  if(l != length(b)){
    stop("vectors have unequal lengths!")
  }
  absence <- 0
  presence <- 0
  gain <- 0
  loss <- 0
  for(i in c(1:l)){
    if(a[i] == 0 & b[i] == 0){
      absence <- absence + 1
    }else if(a[i] == 0 & b[i] == 1){
      gain <- gain + 1
    }else if(a[i] == 1 & b[i] == 0){
      loss <- loss + 1
    }else if(a[i] == 1 & b[i] == 1){
      presence <- presence + 1
    }else{
      stop(paste("presence absence vectors not interpretable at:", i, a[i], b[i]))
    }
  }
  c(presence, absence, gain, loss)
}

boundaryCompl <- function(pabD, pabLengthBins){
  # calculate IES boundary complementarity from a data frame of IES info and lenght bin class
  #extend IES info table with length bin and conservation pattern
  if(length(pabLengthBins) != nrow(pabD)){
    stop()
  }
  DE <- cbind(pabD, lengthBin = unname(pabLengthBins[row.names(pabD)]))
  # exclude floating
  DE <- DE[DE$isFloating == 0, ]
  # permute back sequences(and downstream) reserving lengthBin
  DERI <- permuteBy(DE, c("back", "downstreamFlank"), "lengthBin", ret = "index")
  # calculate complementarity of all
  flankLength <- nchar(DE$upstreamFlank[1])
  windowSize <- nchar(pabD$front)
  matBool <- matrix(ncol = windowSize + flankLength , nrow = nrow(DE))
  matBoolR <- matrix(ncol = windowSize + flankLength , nrow = nrow(DE))
  for(i in c(1:nrow(DE))){
    matBool[i, ] <- complementarity(s2c(paste0(DE$upstreamFlank[i], DE$front[i])),
                                    s2c(paste0(DE$back[i], DE$downstreamFlank[i])))
    matBoolR[i, ] <- complementarity(s2c(paste0(DE$upstreamFlank[i], DE$front[i])),
                                     s2c(paste0(DE$back[DERI[i]], DE$downstreamFlank[DERI[i]])))
  }
  list(matBool, matBoolR, flankLength)
}

getNodePairs <- function(ph, cluster, fromEvent, toEvent){
# get all pairs of parental descendent nodes between speciation events A and B
  DF <- data.frame(cluster = character(0), fromP = character(0), toP = character(0), fromEvent = character(0), toEvent = character(0))
  fromPs <- getEvents(ph, fromEvent)
  toPs <- getEvents(ph, toEvent)
  for(fromP in fromPs){
    for(toP in toPs){
      # check if toP is among Descendants of gromP
      offspringR <- Descendants(ph@phylo, as.numeric(phyldog2r(phyldogNodeId = fromP, cluster = cluster)), type = c("all"))
      if(phyldog2r(phyldogNodeId = toP, cluster = cluster) %in% offspringR){
        DF <- rbind(DF, data.frame(cluster = cluster, fromP =  fromP, toP = toP, fromEvent = fromEvent, toEvent = toEvent, stringsAsFactors = FALSE))
      }
    }
  }
  DF
}

transitions <- function(paV){
  if(length(paV)<2){
    stop("Not enough observations")
  }
  # calculate from a presence absence vector the number of transitions
  presence <- 0
  absence <- 0
  gain <- 0
  loss <- 0
  curV <- paV[1]
  for(i in c(2:length(paV))){
    if(curV == 1 & paV[i] == 1){
      presence <- presence + 1  
    }else if(curV == 0 & paV[i] == 0){
      absence <- absence + 1
    }else if (curV == 0 & paV[i] == 1){
      gain <- gain + 1
    }
    else if (curV == 1 & paV[i] == 0){
      loss <- loss + 1
    }else{
      stop("unknown state")
    }
    curV <- paV[i]
  }
  data.frame(presence = presence, absence = absence, gain = gain, loss = loss, stringsAsFactors = FALSE)
}

#transitions(c(0,1,0,1,1,1,0)) 
# gain 2
# loss 2
# absence 0
# presence 2

rol <- function(bigDF, index){
  # calculate rate of IES loss per IES for each branch type (index)
  sum(pas[index] & !pae[index])/ sum(pas[index]) #loss/present at start
}

tTypeOnTree<- function(ph, orderedMM, fromSpEv, toSpEvs, cluster){
  # calculate transition type from speciation 'fromEv' to speciation event 'toEv' on tree 'ph' (nhx object)
  # given a table 'orderedMM' with node properties (probability of presence)
  bigDF <- data.frame(cluster    = character(0),
                      iesColumn  = character(0),
                      nodeSTartR = character(0),
                      nodeStartPhyldog = character(0),
                      nodeStopR  = character(0),
                      spEventStart = character(0),
                      spEventEnd = character(0),
                      meanStart  = numeric(0),
                      sdStart    = numeric(0),
                      meanStop   = numeric(0),
                      sdStop     = numeric(0),
                      branchType = character(0),
                      transitionTypes = character(0)
  )
  fromNodesP <- getEvents(ph, fromSpEv)
  # find descendents of each node
  for(fromNodeP in fromNodesP){
    #fromNodeP <- fromNodeP[1]
    # translate to R numbering
    fromNodeR <- phyldog2r(phyldogNodeId = fromNodeP, cluster = cluster)
    descendantNodesR <- Descendants(ph@phylo, as.numeric(fromNodeR), type = "all")
    # find type of events
    evS <- getEventsR(RNodeId = descendantNodesR, cluster = cluster, phyldogTree = ph)
    # keep nodes of appropriate speciation event
    for(toSpEv in toSpEvs){
      toNodesR <- descendantNodesR[which(evS$Ev == "S" & evS$S == toSpEv)]
      # next if there are no toNodesR
      # some events cannot be classified as belonging to one type of branch because of gene loss
      if(length(toNodesR) == 0){
        next
      }
      toNodesP <- r2phyldog(RNodeId = toNodesR, cluster)
      fromDF <- orderedMM[orderedMM$phyldog == fromNodeP, c("iesColumn", "mean", "sd", "r", "phyldog")]
      toDF   <- orderedMM[orderedMM$phyldog %in% toNodesP, c("iesColumn", "mean", "sd", "r", "phyldog")]
      DF <- dplyr::full_join(fromDF, toDF, by = "iesColumn")
      DF <- cbind(DF, cluster = cluster, spEventStart = fromSpEv, spEventEnd = toSpEv,stringsAsFactors = FALSE)
      names(DF) <- c("iesColumn", "meanStart", "sdStart", "nodeStartR", "nodeStartPhyldog", "meanStop", "sdStop", "nodeStopR", "nodeStopPhyldog", "cluster", "spEventStart", "spEventEnd") # give nice names
      DF <- DF[,c("cluster", "iesColumn", "nodeStartR", "nodeStartPhyldog", "nodeStopR", "nodeStopPhyldog", "spEventStart", "spEventEnd", "meanStart", "sdStart", "meanStop", "sdStop")] # reodrer
      branchType <- paste0(DF$spEventStart, DF$spEventEnd)
      bigDF <- rbind(bigDF, cbind(DF, branchType, transitionTypes = transitionTypes(DF)))
    }
  }
  bigDF
}
  

transitionTypes <- function(DF){
  # calculate transition type from columns of a data.frame
  tt <- character(nrow(DF))
  for(i in c(1:nrow(DF))){
    tt[i] <- transitionType(DF$meanStart[i], DF$meanStop[i], DF$sdStart[i], DF$sdStop[i])
  }
  tt
}

transitionType <- function(meanFrom, meanTo, sdFrom, sdTo){
  # function that classifies the transition type based on probabilities of the two nodes
  # possibly a better approach would be to do a stochastic simulation on the tree and get probabilities of events on branches
  # possible values are absence, presence, gain, loss, notSignificant
  # for absence both means should be < 0 + interval
  # for presence both means should be > 1 - interval
  interval <- 0.05
  from <- almostP(meanFrom, interval)
  to <- almostP(meanTo, interval)
  if(from == 0 & to == 0){
    return("absent")
  }else if(from == 0 & to == 1){
    return("gain")
  }else if(from == 1 & to == 1){
    return("presence")
  }else if(from == 1 & to == 0){
    return("loss")
  }else if(from == -1 | to == -1){
    return("notSignificant")
  }else{
    stop()
  }
}

almostP <- function(prob, interval){
  # is the probability almost 1 or almost 0
  if(prob < interval){
    return(0)
  }else if(prob > 1 - interval){
    return(1)
  }else{
    return(-1)
  }
}

r2phyldog <- function(RNodeId, cluster){
  # translate node id from R to phyldog
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(RNodeId, nodeDictionary$r[clusterIndex])
  (nodeDictionary[clusterIndex, "phyldog"])[index]
}

phyldog2r <- function(phyldogNodeId, cluster){
  # translate node id from phyldog to R
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(phyldogNodeId, nodeDictionary$phyldog[clusterIndex])
  (nodeDictionary[clusterIndex, "r"])[index]
}

phyldog2rb <- function(phyldogNodeId, cluster){
  # translate node id from phyldog to revBayes
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(phyldogNodeId, nodeDictionary$phyldog[clusterIndex])
  (nodeDictionary[clusterIndex, "rb"])[index]
}

rb2r <- function(rbNodeId, cluster){
  # translate node id from rb (revBayes) to R
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(rbNodeId, nodeDictionary$rb[clusterIndex])
  (nodeDictionary[clusterIndex, "r"])[index]
}

rb2phyldog <- function(rbNodeId, cluster){
  # translate node id from rb (revBayes) to PHYLDOG
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(rbNodeId, nodeDictionary$rb[clusterIndex])
  (nodeDictionary[clusterIndex, "phyldog"])[index]
}

r2rb <- function(RNodeId, cluster){
  # translate node id from R to rb (revBayes)
  clusterIndex <- which(nodeDictionary$cluster %in% as.character(cluster))
  index <- match(RNodeId, nodeDictionary$r[clusterIndex])
  (nodeDictionary[clusterIndex, "rb"])[index]
}

getEventsR <- function(RNodeId, cluster, phyldogTree){
  # get type of event from an R node
  phyldogNodeId <- r2phyldog(RNodeId, cluster)  
  evS <- ph@nhx_tags[match(phyldogNodeId, ph@nhx_tags$ND), c("Ev","S")]
  evS
}

getEvents <- function(phyldogTree, eventType){
  # find which node(s) corresponds to event of type eventType {0,1,2,3..} for Ev = S, and {D} for duplications
  # returns phyldogNodeId
  spEventsV <- extractEvents(phyldogTree)
  phyldogTree@nhx_tags$ND[which(spEventsV == as.character(eventType))]
}

extractEvents <- function(ph){
  # Input an nhx object
  # Returns a vector of events
  spEventsV <- character(nrow(ph@nhx_tags))
  for(i in 1:nrow(ph@nhx_tags)){
    if(ph@nhx_tags$Ev[i] == "S"){
      # speciation event
      spEventsV[i] <- ph@nhx_tags$S[i]
    }else if(ph@nhx_tags$Ev[i] == "D"){
      # duplication event
      spEventsV[i] <- "D"
    }else{
      stop("unknown event")
    }
  }
  spEventsV
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

pairwiseIdentity <- function(a,b){
  # function that returnes a boolean vector with identities of the character matrix with two input strings
  if(length(a) != length(b)){
    stop()
  }
  ac <- toupper(a)
  bc <- toupper(b)
  ac==bc
}

complementarity <- function(a, b){
  # calculate per locus complementarity of two sequences, by reverse complementing the second and calculating identity
  if(length(a) != length(b)){
    stop()
  }
  b <- rev(comp(b))
  pairwiseIdentity(a, b)
}

combS <- function(a,b){
  # function that finds all combinations of vectors a and b
  la <- length(a)
  lb <- length(b)
  m <- min(la, lb)
  paired <- matrix(nrow=la*lb,ncol=2)
  counter <- 1
  for(i in c(1:la)){
    for(j in c(1:lb)){
      paired[counter,] <- c(a[i],b[j])
      counter <- counter + 1
    }
  }
  paired
}

prepareLOGO <- function(l){
  aln <- seqinr::as.alignment(nb = length(l),nam=as.character(c(1:length(l))),seq = c(l))
  m <- as.matrix.alignment(aln)
  as <- numeric(0)
  ts <- numeric(0)
  cs <- numeric(0)
  gs <- numeric(0)
  for(i in c(1:length(m[1,]))){
    mt <- table(m[,i])
    as <- append(as,mt["A"])
    cs <- append(cs,mt["C"])
    ts <- append(ts,mt["T"])
    gs <- append(gs,mt["G"])
  }
  #NA to 0
  as[is.na(as)]<-0
  ts[is.na(ts)]<-0
  cs[is.na(cs)]<-0
  gs[is.na(gs)]<-0
  df <- data.frame("A"=as,"C"=cs,"G"=gs,"T"=ts)
  seqLogo(makePWM(t(df/rowSums(df))))
}

#remove NAs and names from vectors and
nadel <- function(x){
  as.vector(x[!is.na(x)])
}
# 
# x <- data.frame(letters = letters[1:10], numbers1 = c(1:10), numbers2 = c(1,1,1,1,1,1,1,2,2,2), ignored = letters[1:10], numbers3 = c(NA,NA,NA,2,3,3,4,4,5,5))
# col2permute <- c("letters", "numbers1")
# col2keepFixed <- c("numbers2", "numbers3")

permuteBy <- function(x, col2permute, col2keepFixed, ret = "index"){
  # function that permutes values in data.frame column but keeping values in another column fixed
  # col2permute can be a column name of the x data.frame  or a character vector with multiple such names
  # col2keepFixed can be a column name of the x data.frame
  if(!is.data.frame(x)){
    stop("x must be a data.frame")
  }
  permutedI <- numeric(nrow(x))
  #make a key from all columns that we need to keep
  if(length(col2keepFixed) >  1){
    key <-  apply(x[ , col2keepFixed], 1, paste0, collapse =" ")
  }else{
    key <- x[, col2keepFixed]
  }
  #find unique combinations
  lev <- unique(key)
  #find index of unique combinations
  for(l in lev){
    index <- which(key == l)
    if(length(index) == 1){
      newIndex = index
    }else{
      newIndex <- sample(index, length(index))    
    }
    permutedI[index] <- newIndex
  }
  for(i in col2permute){
    x[i] <- x[permutedI, i]
  }
  if(ret == "index"){
    return(permutedI)
  }else if(ret == "data.frame"){
    return(x)
  }else{
    stop("return error")
  }
}

compl <- function(gc){
  # The expected complementarity of random sequences with gc GC%
  at <- 1 - gc  
  at*at/2+gc*gc/2
}

appa <- function(x){
  # approximate presence-absence
  # The function determines the state of an IES in a species from character vector containing information from all paralogs
  # INPUT: a character string with the following options
  #       0:  absence
  #       NA: unknown IES state
  #       !=0: presence
  # OUTPUT: NA: no IES information
  #          0: all absent
  #          1: at least one present
  #         -1: empty input
  if(length(x) == 0){
    return(-1)
  }
  if(all(is.na(x))){
    return(NA)
  }
  if(all(x == 0, na.rm = TRUE)){
    return(0)
  }
  if(any(x != 0)){
    return(1)
  }  
}
#sp <- spGT
presenceAbsence <- function(sp){
  # Function to assign patterns of presence / absence of IES in subtrees.
  # In the case of multiple paralogs, presence is counted as at least one paralog having an IES
  load("~/data/IES_data/rdb/charMats")
  # add to charMats a column with in-/outgroup/NA and species
  extMat <- join(charMats, sp, by = "geneId")
  biI <- which(extMat[, "species"] == "Paramecium_biaurelia")
  teI <- which(extMat[, "species"] == "Paramecium_tetraurelia")
  seI <- which(extMat[, "species"] == "Paramecium_sexaurelia")
  caI <- which(extMat[, "species"] == "Paramecium_caudatum")
  uniqueSubtrees <- unique(extMat$subtree)
  # the total number of patterns will be the product of the number of subtrees by the number of columns in each one
  countRows <- extMat[,c("column","subtree")]
  countRows <- countRows[!is.na(countRows$subtree), ]
  countRows <- unique(paste(countRows[,1], countRows[,2]))
  patternsM <- matrix(ncol = 4, nrow = length(countRows)) # minimum length of matrix
#  patternsV <- character(length(uniqueSubtrees))
  counter <- 1
  #for(clustcol in unique(homIES)){
  for(subtree in uniqueSubtrees){
    #subtree <- uniqueSubtrees[2]
    if(is.na(subtree)){
      next
    }
    indexCM <- which(extMat$subtree == subtree) # create and index for the current subtree
    homologIES <- unique(extMat[indexCM, "column"])
    for(column in homologIES){
      # column <- homologIES[1]                 # and one for each column
      indexCMC <- which(extMat$subtree == subtree & extMat$column == column)
      # for each species find if at least one IES
      patternsM[counter, ] <- c(appa(extMat[intersect(biI, indexCMC), "ies"]),
                                appa(extMat[intersect(teI, indexCMC), "ies"]),
                                appa(extMat[intersect(seI, indexCMC), "ies"]),
                                appa(extMat[intersect(caI, indexCMC), "ies"]))
      counter <- counter + 1
    }
  }
  patternsM
}


gene2species <- function(string){
  # get a character vector of species names from 
  # a character vector of gene names
  knownSpeciesNames <- c("Paramecium_primaurelia", "Paramecium_biaurelia", "Paramecium_tetraurelia", "Paramecium_pentaurelia", "Paramecium_sexaurelia", "Paramecium_octaurelia", "Paramecium_tredecaurelia", "Paramecium_sonneborni",  "Paramecium_caudatum", "Tetrahymena_thermophila")
  names(knownSpeciesNames) <-  c("PPRI", "PBIA", "PTET", "PPEN", "PSEX", "POCT", "PSON", "PTRE", "PCAU", "TTHE")
  speciesNames <- knownSpeciesNames[substr(string, 0, 4)]
  if(any(is.na(speciesNames))){
    stop(paste("unknown name", string[which(is.na(speciesNames))]))
  }
  speciesNames
}

extractGenes <- function(ggff){
  # extract information for genes from a gff3 file
  geneIndex <- which(ggff$V3 == "gene")
  geneIds <- ggff[geneIndex, 9]
  matched <- regexpr("ID=[^;]*;", geneIds, perl = TRUE)
  # clean names fron ID=...;
  geneIds <- gsub("(ID=)|;", "", regmatches(geneIds, matched), perl = TRUE)
  # G to P 
  geneIds <- sub("G","P",geneIds)
  scaffold <- ggff[geneIndex, 1]
  strand <- ggff[geneIndex, 7]
  begin <- ggff[geneIndex, 4]
  end <- ggff[geneIndex, 5]
  df <- data.frame(id = geneIds, begin = begin, end = end, strand = strand, scaffold = scaffold, stringsAsFactors = FALSE)
  #order by scaffold, then by start and then by end
  df[order(df$scaffold, df$begin, df$end), ]
}


overlapping <- function(a, b){
  # function that tests if two elements a and b are overlapping
  # strand is ignored as we care for both strands
  if(b$begin >= a$begin & b$begin <= a$end){  # B starts within A
    return(1) 
  }
  if(b$end >= a$begin & b$end <= a$end){      # B ends within A
    return(1)
  }
  if(b$end >= a$end & b$begin <= a$begin){
    return(1)                             # A is within B
  }
  return(0)
}

closest <- function(DF, index, direction){
  # look up or down in sorted data.frame for closest gene to the one at row index
  # return index of closest gene, distance and relative orientation
  # if no closest gene present (first or last gene) return -1
  # if focus gene is the only in scaffold also return -1
  # if there is overlap with closest genes return -1
  scaffold <- DF[index, "scaffold"]
  if(direction == "upstream"){
    ua <- -1
  }else if(direction == "downstream"){
    ua <- +1
  }else{
    stop("unknown direction")
  }
  candidate <- index + ua
  if(length(DF[candidate, "scaffold"]) == 0){
    # candidate gene does not exist (focus gene is first in scaffold)
    return(-1)
  }
  if(is.na(DF[candidate, "scaffold"])){
    # last gene reached
    return(-1)
  }
  if(DF[index, "scaffold"] != DF[candidate, "scaffold"]){
    # genes are on different scaffolds (focus gene is the only one in scaffold)
    return(-1)
  }
  if(overlapping(DF[index, ], DF[candidate, ])){
    return(-1) # overlapping gene
  }
  if(direction == "downstream"){
    return(c(index = candidate, geneDist(DF[index, ], DF[candidate, ])))
  }else if(direction == "upstream"){
    return(c(index = candidate, geneDist(DF[candidate, ], DF[index, ]))) # distance is always positive
  }else{
    stop("unknown direction")
  }
}

geneDist <- function(a, b){
  # The genomic coordinates (begin, end) are not dependent on orientation, begin is always more upstream than end
  # genomic distance and relative orientation calculation assuming b is downstream a.
  if(a$scaffold != b$scaffold){
    print(a)
    stop("Can't compute distance of genes from different chromosomes")
  }
  if(a$strand == "+" & b$strand == "+"){
    #d <- b$begin - a$end - 1 # -> ->
    relOri <- "same"
  }else if(a$strand == "-" & b$strand == "-"){ # <- <-
    #d <- b$end - a$begin - 1
    relOri <- "same"
  }else if(a$strand == "+" & b$strand == "-"){
    #d <- b$end - a$end - 1
    relOri <- "in" # heads in ->  <-
  }else if(a$strand == "-" & b$strand == "+"){
    #d <- b$begin - a$begin 
    relOri <- "out" # heads out <- ->
  }else{
    stop("unknown orientation(s)")
  }
  d <- b$begin - a$end
  return(c(distance = d, relOri = relOri, begin = a$end, end = b$begin))
}


genePairs <- function(DF){
  # from the output of extractGenes function find all pairs of consecutive genes, their distances and relative orientation, excluding overlapping genes
  l <- nrow(DF)
  genePairsM <- matrix(ncol = 8, nrow = 2 * l)
  counter <- 1
  for(i in c(1:l)){
    scaffold <- DF[i, "scaffold"]
    if(length(which(DF[, "scaffold"] == scaffold)) == 1){
      next # skip scaffolds with only one gene
    }
    upNeighbor <- closest(DF, i, "upstream")
    downNeighbor <- closest(DF, i, "downstream")
    if(length(upNeighbor) == 1){
      # either first gene or overlapping another gene
    }else{
      genePairsM[counter, ] <- c(DF[i, "id"], DF[as.numeric(unname(upNeighbor["index"])), "id"], "upstream", unname(upNeighbor[c("distance", "relOri", "begin", "end")]), scaffold)
      counter <- counter + 1
    }
    if(length(downNeighbor) == 1){
      # either last gene in scaffold or overlapping another one
    }else{
      genePairsM[counter, ] <- c(DF[i, "id"], DF[as.numeric(unname(downNeighbor["index"])), "id"], "downstream", unname(downNeighbor[c("distance", "relOri", "begin", "end")]), scaffold)
      counter <- counter + 1
    }
  }
  genePairsM[apply(genePairsM, 1, notAllNA), ] # return matrix without the extra rows of NAs
}

# small functions useful for use in apply
notAllNA <- function(a){
  !all(is.na(a))
}

naSum <- function(a){
  # the sun of non NA values
  sum(is.na(a))
}

iesInIntergenic <- function(synblockDF, pabD, speciesName){
  # function that finds if IES are within an intergenic region of a given species
  # as input is given the synblockDF data.frame of blocks of intergenic regions, the pabD data.frame of iesInfo
  # and the species Name for which to perform the search
  # The value returned is a boolean vector with the same length as the rows of synblockDF
  if(!speciesName %in% c("Paramecium_biaurelia", "Paramecium_tetraurelia", "Paramecium_sexaurelia", "Paramecium_caudatum")){
    stop(paste("error in speciesNane", speciesName))
  }
  hasIES <- rep(NA, nrow(synblockDF))
  spV <- unname(gene2species(synblockDF$upanchor))
  pabI <- which(spV == speciesName)
  # all scaffolds that have IES or genes
  scaffolds <- unique(c(pabD$scaffold, unique(synblockDF$scaffold[pabI])))
  for(i in c(1:length(scaffolds))){
    # index of intergenic regions in the current scaffold
    intergenicScaffoldI <- which(synblockDF$scaffold == scaffolds[i])
    # build IRanges objects for intergenic regions in each species
    spscI <- intersect(pabI, intergenicScaffoldI)
    intergenicIR <- IRanges(start = as.numeric(synblockDF[spscI, "begin"]),
                            end   = as.numeric(synblockDF[spscI, "end" ])
    )
    # for each species buildindex of IES in the current scaffold
    iesI <- which(pabD$scaffold == scaffolds[i])
    # build corresponding IRanges object
    iesIR <- IRanges(start = pabD$start[iesI],
                     end   = pabD$end[iesI]
    )
    # find if intergenic there are IES
    hasIES[spscI] <- intergenicIR %over% iesIR
  }
  hasIES
}


modif.write.nexus.data <- function (x, file, format = "dna", datablock = TRUE, interleaved = TRUE, charsperline = NULL, gap = NULL, missing = NULL){
  # minor modification in ape::write.nexus.data to save standard and restriction type data  
  format <- match.arg(toupper(format), c("DNA", "PROTEIN", "STANDARD", "RESTRICTION"))
  indent <- "  "
  maxtax <- 5
  defcharsperline <- 80
  defgap <- "-"
  defmissing <- "?"
  if (is.matrix(x)) {
    if (inherits(x, "DNAbin")) 
      x <- as.list(x)
    else {
      xbak <- x
      x <- vector("list", nrow(xbak))
      for (i in seq_along(x)) x[[i]] <- xbak[i, ]
      names(x) <- rownames(xbak)
      rm(xbak)
    }
  }
  ntax <- length(x)
  nchars <- length(x[[1]])
  zz <- file(file, "w")
  if (is.null(names(x))) 
    names(x) <- as.character(1:ntax)
  fcat <- function(..., file = zz) cat(..., file = file, sep = "", 
                                       append = TRUE)
  find.max.length <- function(x) max(nchar(x))
  print.matrix <- function(x, dindent = "    ") {
    Names <- names(x)
    printlength <- find.max.length(Names) + 2
    if (!interleaved) {
      for (i in seq_along(x)) {
        sequence <- paste(x[[i]], collapse = "")
        taxon <- Names[i]
        thestring <- sprintf("%-*s%s%s", printlength, 
                             taxon, dindent, sequence)
        fcat(indent, indent, thestring, "\n")
      }
    }
    else {
      ntimes <- ceiling(nchars/charsperline)
      start <- 1
      end <- charsperline
      for (j in seq_len(ntimes)) {
        for (i in seq_along(x)) {
          sequence <- paste(x[[i]][start:end], collapse = "")
          taxon <- Names[i]
          thestring <- sprintf("%-*s%s%s", printlength, 
                               taxon, dindent, sequence)
          fcat(indent, indent, thestring, "\n")
        }
        if (j < ntimes) 
          fcat("\n")
        start <- start + charsperline
        end <- end + charsperline
        if (end > nchars) 
          end <- nchars
      }
    }
  }
  fcat("#NEXUS\n[Data written by write.nexus.data.R, ", date(), 
       "]\n")
  NCHAR <- paste("NCHAR=", nchars, sep = "")
  NTAX <- paste0("NTAX=", ntax)
  DATATYPE <- paste0("DATATYPE=", format)
  if (is.null(charsperline)) {
    if (nchars <= defcharsperline) {
      charsperline <- nchars
      interleaved <- FALSE
    }
    else charsperline <- defcharsperline
  }
  if (is.null(missing)) 
    missing <- defmissing
  MISSING <- paste0("MISSING=", missing)
  if (is.null(gap)) 
    gap <- defgap
  GAP <- paste0("GAP=", gap)
  INTERLEAVE <- if (interleaved) 
    "INTERLEAVE=YES"
  else "INTERLEAVE=NO"
  if (datablock) {
    fcat("BEGIN DATA;\n")
    fcat(indent, "DIMENSIONS ", NTAX, " ", NCHAR, ";\n")
    fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", 
         GAP, " ", INTERLEAVE, ";\n")
    fcat(indent, "MATRIX\n")
    print.matrix(x)
    fcat(indent, ";\nEND;\n\n")
  }
  else {
    fcat("BEGIN TAXA;\n")
    fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
    fcat(indent, "TAXLABELS\n")
    fcat(indent, indent)
    j <- 0
    for (i in seq_len(ntax)) {
      fcat(names(x[i]), " ")
      j <- j + 1
      if (j == maxtax) {
        fcat("\n", indent, indent)
        j <- 0
      }
    }
    fcat("\n", indent, ";\n")
    fcat("END;\n\nBEGIN CHARACTERS;\n")
    fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")
    fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", DATATYPE, 
         " ", INTERLEAVE, ";\n")
    fcat(indent, "MATRIX\n")
    print.matrix(x)
    fcat(indent, ";\nEND;\n\n")
  }
  close(zz)
}

save2nexus <- function(geneFamily){
  # prepare IES character matrices to be writen in a nexus file by
  DF <- data.frame(charMats[charMats$cluster == geneFamily, ], stringsAsFactors = FALSE)
  DFcasted <- dcast(DF, geneId ~ column, value.var = "ies")
  list4nexus <- list()
  for(i in c(1:nrow(DFcasted))){
    sequenceS <- DFcasted[i, -1]
    sequenceS[sequenceS != "0"] <- "1"
    sequenceS[is.na(sequenceS)] <- "?"
    sequenceV <- as.vector(unlist((sequenceS)))
    sequenceL <- list(sequenceV)
    names(sequenceL)<-DFcasted$geneId[i]
    list4nexus <- append(list4nexus, sequenceL)
  }
  homIESids <- names(DFcasted[2:ncol(DFcasted)])
  list(list4nexus, homIESids)
}

brfilt <- function(brl){
  # for each branch filter-out gene families with extremely long branches (above the 90th quantile of the log branch length)
  outliers <- character()
  for(i in c(1:ncol(brl))){
    # skip all NA
    if(all(is.na(brl[, i]))){
      next
    }
    r <- quantile(log10(brl[,i]), probs = 0.9, names = FALSE)
    outliers <- append(outliers, row.names(brl)[log10(brl[, i]) > r[1]])
  }
  unique(outliers)
}

getCoreCluster <- function(clusterId, iess, iesdb){
  # select the core of a cluster by filtering out floating, merged and not in peak elements
  cl <- iess$V2[iess$V1 == clusterId]
  notFloatingI <- which(iesdb$isFloating !=1)
  notMergedI <- which(iesdb$merged !=1)
  inCl <- which(iesdb$id%in%cl)
  ind <- which(d$cluster == clusterId)
  if (d$peak[ind]){
    minL <- d$lengthL[ind]
    maxL <- d$lengthU[ind]
    peakI <- which(iesdb$length >= minL & iesdb$length < maxL)
    includeIES <- intersect(intersect(notFloatingI, notMergedI), intersect(inCl, peakI))
  }else{
    includeIES <- intersect(intersect(notFloatingI, notMergedI), inCl)
  }
  iesdb[includeIES,]
}

goodHSP<- function(w, b, e, l){
  # a good HSP has begining and end of query matching with begin and end of subject (or vice versa) for a window of size w
  (b <= w & e >= l - w + 1) | (b >= l - w + 1 & e <= w)
}

HSPtype <- function(w, b, e, l){
  # classify HSPs as
  retV <- matrix(nrow = length(b), ncol = (4))
  aB <- b <= w
  retV[aB, 1] <- TRUE
  retV[!aB, 1] <- FALSE
  bB <- e >= l - w + 1 
  retV[bB, 2] <- TRUE
  retV[!bB, 2] <- FALSE
  cB <- b >= l - w + 1
  retV[cB, 3] <- TRUE
  retV[!cB, 3] <- FALSE
  dB <- e <= w 
  retV[dB, 4] <- TRUE
  retV[!dB, 4] <- FALSE
  return(retV)
  # 1. HSP with beginning matching beginning
  # 2. HSP with end matching end
  # 3. both of the above
  # 4. HSP with beginning matching end
  # 5. HSP with end matching beginning
  # 6, both of the above
  #IES with 1 & 2 | 3 & 4
}

source("~/projects/IES/src/testSharedFunctions.R")
