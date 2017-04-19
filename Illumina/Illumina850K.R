library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Annotations only
annot.cg <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19@data$Other

z <- gzfile("Illumina850K.csv.gz")
write.csv(annot.cg, z, quote = FALSE)

# hg19 coordinates
BED.cg <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19@data$Locations
BED.cg <- data.frame(chr = BED.cg$chr, start = BED.cg$pos + 1, end = BED.cg$pos + 1, name = rownames(BED.cg), score = 1000, strand = BED.cg$strand)

write.table(BED.cg, "Illumina850K.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
