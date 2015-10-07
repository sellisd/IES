# classify IES to length classes
# define length bins based on observations.
library(IRanges)
maxLength <- 10000
peakX <- c(27,38,45,seq(56, maxLength, by = 10))
startX <- c(18,34,41,seq(52, maxLength, by = 10))
endX <- startX[-1]
lengthBins <- IRanges(start = startX[-length(startX)], end = endX-1)
# read ies lengths
load("~/data/IES_data/rdb/iesInfo")
pbiIesLengthsIR <- IRanges(start = pbiD$length, end = pbiD$length)
pbiBins <- findOverlaps(pbiIesLengthsIR,lengthBins)
pbiLengthBins <- pbiBins@subjectHits
names(pbiLengthBins) <- row.names(pbiD)

pteIesLengthsIR <- IRanges(start = pteD$length, end = pteD$length)
pteBins <- findOverlaps(pteIesLengthsIR,lengthBins)
pteLengthBins <- pteBins@subjectHits
names(pteLengthBins) <- row.names(pteD)

pseIesLengthsIR <- IRanges(start = pseD$length, end = pseD$length)
pseBins <- findOverlaps(pseIesLengthsIR,lengthBins)
pseLengthBins <- pseBins@subjectHits
names(pseLengthBins) <- row.names(pseD)

pcaIesLengthsIR <- IRanges(start = pcaD$length, end = pcaD$length)
pcaBins <- findOverlaps(pcaIesLengthsIR,lengthBins)
pcaLengthBins <- pcaBins@subjectHits
names(pcaLengthBins) <- row.names(pcaD)

save(pbiLengthBins, pteLengthBins, pseLengthBins, pcaLengthBins, file = "~/data/IES_data/rdb/lengthBins")
