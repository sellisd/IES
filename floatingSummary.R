# Count percentage of floating IES for each species

PbiIES <- read.table("~/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.gff3")
PbiIESfl <- read.table("~/data/IES_data/pbiaurelia/internal_eliminated_sequence_MIC_biaurelia.pb_V1-4.fl.gff3")
PteIES <- read.table("~/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.gff3")
PteIESfl <- read.table("~/data/IES_data/ptetraurelia/internal_eliminated_sequence_PGM_IES51.pt_51.fl.gff3")  
PseIES <- read.table("~/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.gff3")
PseIESfl <- read.table("~/data/IES_data/psexaurelia/internal_eliminated_sequence_MIC_sexaurelia.ps_AZ8-4.fl.gff3")  
PcaIES <- read.table("~/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.gff3")
PcaIESfl <- read.table("~/data/IES_data/pcaudatum_43c3d_annotation_v2.0/PCAUD_MIC10_IES.fl.gff3")  

PbiIESNo <- length(PbiIES$V1)
PteIESNo <- length(PteIES$V1)
PseIESNo <- length(PseIES$V1)
PcaIESNo <- length(PcaIES$V1)
PbiIESNoFl <- sum(paste(PbiIES$V4,PbiIES$V5) != paste(PbiIESfl$V4,PbiIESfl$V5))
PteIESNoFl <- sum(paste(PteIES$V4,PteIES$V5) != paste(PteIESfl$V4,PteIESfl$V5))
PseIESNoFl <- sum(paste(PseIES$V4,PseIES$V5) != paste(PseIESfl$V4,PseIESfl$V5))
PcaIESNoFl <- sum(paste(PcaIES$V4,PcaIES$V5) != paste(PcaIESfl$V4,PcaIESfl$V5))

total <- c(PbiIESNo, PteIESNo, PseIESNo, PcaIESNo)
floating <- c(PbiIESNoFl, PteIESNoFl, PseIESNoFl, PcaIESNoFl)
floatingRatio <- floating/total

floatingT <- data.frame(total = total, floating = floating, floatingP = format(100*floatingRatio,digits = 3), stringsAsFactors = FALSE)
save(floatingT, file = "~/data/IES_data/rdb/floatingSummary")
