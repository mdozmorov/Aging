## load the library
library(FDb.InfiniumMethylation.hg19)
FDb.InfiniumMethylation.hg19
## extract features for use in constructing SummarizedExperiments
## or comparing chip features against other data (e.g. ChIP-seq)
InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
InfiniumMethylation["cg20822990"] # Check random CpG ID
# Create BED-like data frame
Illumina450K.bed <- data.frame(chrom=as.character(seqnames(InfiniumMethylation)),
                                  chromStart=start(InfiniumMethylation),
                                  chromEnd=end(InfiniumMethylation),
                                  name=names(InfiniumMethylation),
                                  score=0,
                                  strand=as.character(strand(InfiniumMethylation)))
# Save to file
write.table(Illumina450K.bed, "data/Illumina450K.bed", sep="\t", quote=F, col.names=F, row.names=F)

## Working with 27K array
hm27.hg19 <- get27k()
summary(hm27.hg19)
# Create BED-like data frame
Illumina27K.bed <- data.frame(chrom=as.character(seqnames(hm27.hg19)),
                              chromStart=start(hm27.hg19),
                              chromEnd=end(hm27.hg19),
                              name=names(hm27.hg19),
                              score=0,
                              strand=as.character(strand(hm27.hg19)))
# Save to file
write.table(Illumina27K.bed, "data/Illumina27K.bed", sep="\t", quote=F, col.names=F, row.names=F)

# Function to read in ID file and write out BED subset of Illumina 450K annotation
subsetIllumina450 <- function(fnamein, fnameout) {
  datain <- read.table(fnamein, sep="\t", header=T, row.names=1)
  #if (setdiff(rownames(datain), Illumina450K.bed$name) <> 0) return(NA)
  write.table(Illumina450K.bed[Illumina450K.bed$name %in% rownames(datain), ], fnameout, sep="\t", quote=F, row.names=F, col.names=F)
  #return(0)
}

# Function to read in ID file and write out BED subset of Illumina 450K annotation
subsetIllumina27 <- function(fnamein, fnameout) {
  datain <- read.table(fnamein, sep="\t", header=T, row.names=1)
  write.table(Illumina27K.bed[Illumina27K.bed$name %in% rownames(datain), ], fnameout, sep="\t", quote=F, row.names=F, col.names=F)
}

subsetIllumina450("data/Hannum//mmc2_Model_Primary_data.txt", "data/Hannum/mmc2_Model_Primary_data.bed")
subsetIllumina450("data/Hannum//mmc2_Model_All_data.txt", "data/Hannum/mmc2_Model_All_data.bed")
subsetIllumina450("data/Hannum//mmc3_Model_Breast.txt", "data/Hannum/mmc3_Model_Breast.bed")
subsetIllumina450("data/Hannum//mmc3_Model_Kidney.txt", "data/Hannum/mmc3_Model_Kidney.bed")
subsetIllumina450("data/Hannum//mmc3_Model_Lung.txt", "data/Hannum/mmc3_Model_Lung.bed")

subsetIllumina27("data/Rakyan//Supp_Table_3.txt", "data//Rakyan/Supp_Table_3.bed")

subsetIllumina450("data/Alisch//Alisch_et_al_Sup_Table2.txt", "data/Alisch//Alisch_et_al_Sup_Table2.bed")
subsetIllumina450("data/Alisch//Alisch_et_al_Sup_Table2_pos.txt", "data/Alisch//Alisch_et_al_Sup_Table2_pos.bed")
subsetIllumina450("data/Alisch//Alisch_et_al_Sup_Table2_neg.txt", "data/Alisch//Alisch_et_al_Sup_Table2_neg.bed")

