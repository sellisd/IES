# calculate complementarity of IES boundaries and flanking regions
load("~/data/IES_data/rdb/iesInfo")
load("~/data/IES_data/rdb/lengthBins")
source("~/projects/IES/src/sharedFunctions.R")
pbiL <- boundaryCompl(pbiD, pbiLengthBins)
pteL <- boundaryCompl(pteD, pteLengthBins)
pseL <- boundaryCompl(pseD, pseLengthBins)
pcaL <- boundaryCompl(pcaD, pcaLengthBins)
save(pbiL, pteL, pseL, pcaL, file = "~/data/IES_data/rdb/complementarity")
