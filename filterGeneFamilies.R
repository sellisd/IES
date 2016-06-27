#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# filter gene families and create plots of branch length distributions

if (length(args) != 3) {
  stop("Missing argument(s), input output")
}

inputF <- args[1] # input "~/data/IES/analysis/tables/normBrLensC.tab"
pdfOUT <- args[2] # pdf output "~/data/IES/analysis/figures/brLenDistr.pdf"
brlenPath <- args[3] # filtered path output "/home/dsellis/data/IES/analysis/brlen/"

d <- read.table(inputF, header = TRUE)
d <- d[,!is.na(d[1,])] # drop columns with NA in the first entry (should be full of NAs)

pdf(pdfOUT)
par(mfcol=c(4,4), las = 1, bty = "l")

# for each branch filter-out gene families with extremely long branches (above the 90th quantile of the log branch length)
# plot distribution of log lengths for each branch and manually determine outliers (keep the smallest peak that is not almost zero)
outliers <- character()
boeN <- ncol(d)
mins <- numeric(boeN)
maxs <- numeric(boeN)
for(i in c(1:boeN)){
  r <- quantile(log10(d[, i]), probs = c(0.9), names = FALSE)
  outliers <- append(outliers, row.names(d)[log10(d[, i]) > r[1]])
  h <- hist(log10(d[,i]), breaks = 50, plot = FALSE)
  plot(h$mids, h$counts, main = names(d)[i], log = "", type = "h", xlab = "log(branch length)", ylab = "counts")
  abline(v = r, col = "red")
  #l <- locator(2)
  #abline(v = l$x, col = "red")
  #mins[i] <- l$x[1]
  #maxs[i] <- l$x[2]
  #outliers <- append(outliers, row.names(d)[(log10(d[, i]) < mins[i]) | (log10(d[, i]) > maxs[i])])
  #print(l$x)
  #  abline(v = log10(r), col = "red")
}
brlgf <- paste0("cluster.", setdiff(row.names(d), outliers), ".aln.fa")

if(!dir.exists(brlenPath)){
  dir.create(brlenPath)
}
filteredPath <- "/home/dsellis/data/IES/analysis/msas/filtered/"
#files <- dir(path= filteredPath, pattern="*.aln.fa$")

from <- paste(filteredPath, brlgf, sep= "")
to <- paste(brlenPath, brlgf, sep = "")
r <- file.copy(from, to)
print(all(r) == TRUE)
dev.off()