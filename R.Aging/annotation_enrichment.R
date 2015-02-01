# Extract and append multiple values embedded in rows
#
# data: data.frame
# col: column name containing embedded values
# sep: regular expression to split column by
#
# df <- data.frame(key = c("a", "a;b", "a;b;c"), val = 1:3)
# unembed(df, "key", ";")

unembed <- function(dat, col, sep, ...) {
  
  stopifnot(is.data.frame(dat))
  col_i <- which(names(dat) == col)
  
  dat[[col]] <- as.character(dat[[col]])
  pieces <- strsplit(dat[[col]], sep, fixed=T)
  ns <- vapply(pieces, length, integer(1))
  
  #   structure(data.frame(unlist(pieces), 
  #                        data[rep(seq_along(ns), ns), -col_i]), 
  #                        names = c(col, names(data)[-col_i]))
  data.unembed <- cbind(unlist(pieces), dat[rep(seq_along(ns), ns), -col_i])
  names(data.unembed) <- c(col, names(dat)[-col_i])
  return(data.unembed)
}

## =========================================
# ChromStates
#  Preprocess annotation matrix
mtx.chrom <- read.table("../GenomeRunner/Annotation_BroadHmm_H1hesc/all_annot_chromstates.txt", sep="\t", header=T, stringsAsFactors=F)
mtx.chrom <- mtx.chrom[, c(-1, -2, -3)]
rownames(mtx.chrom) <- make.names(mtx.chrom[, 1], unique=T)
mtx.chrom <- as.matrix(mtx.chrom[, -1])
# Flooring
mtx.chrom[mtx.chrom > 1] <- 1
# Get enrichment summary matrix
mtx.chrom.stats <- as.matrix(read.table("../GenomeRunner/Annotation_BroadHmm_H1hesc/summary_all_embryo_chrom_stats.txt", sep="\t", row.names=1, header=T, stringsAsFactors=F))

# Histones
# Preprocess annotation matrix
mtx.histone <- read.table("../GenomeRunner/Annotation_Histone_H1hesc/all_annot_histone.txt", sep="\t", header=T, stringsAsFactors=F)
mtx.histone <- mtx.histone[, c(-1, -2, -3)]
rownames(mtx.histone) <- make.names(mtx.histone[, 1], unique=T)
mtx.histone <- as.matrix(mtx.histone[, -1])
# Flooring
mtx.histone[mtx.histone > 1] <- 1
# Get enrichment summary matrix
mtx.histone.stats <- as.matrix(read.table("../GenomeRunner/Annotation_Histone_H1hesc/summary_all_embryo_hist_stats.txt", sep="\t", row.names=1, header=T, stringsAsFactors=F))

# TFBSs
#  Preprocess annotation matrix
mtx.tfbs <- read.table("../GenomeRunner/Annotation_Tfbs_H1hesc/all_annot_tfbs.txt", sep="\t", header=T, stringsAsFactors=F)
mtx.tfbs <- mtx.tfbs[, c(-1, -2, -3)]
rownames(mtx.tfbs) <- make.names(mtx.tfbs[, 1], unique=T)
mtx.tfbs <- as.matrix(mtx.tfbs[, -1])
# Special processing of TFBS matrix
tmp <- t(mtx.tfbs)
tmp <- data.frame(tmp, table=rownames(tmp))
tmp1 <- left_join(tmp, gfAnnot[, c(1, 3, 5, 2)], by=c("table" = "V1"))
tmp1.1 <- (tmp1[, c( -27635, -27637, -27636)])
tmp1.1 <- as.matrix(tmp1.1[, -1])
tmp1.5 <- tmp1$V5
tmp2 <- aggregate(data.frame(tmp1.1), by=list(as.character(tmp1$V5)), sum)
tmp2 <- t(tmp2)
colnames(tmp2) <- tmp2[1, ]
tmp2 <- tmp2[-1, ]
mtx.tfbs <- as.matrix(tmp2)
class(mtx.tfbs) <- "numeric"
# colnames(mtx.tfbs) <- make.names(colnames(mtx.tfbs))
# Flooring
mtx.tfbs[mtx.tfbs > 1] <- 1
# Get enrichment summary matrix
mtx.tfbs.stats <- as.matrix(read.table("../GenomeRunner/Annotation_Tfbs_H1hesc/summary_all_embryo_tfbs_stats.txt", sep="\t", row.names=1, header=T, stringsAsFactors=F))

###==================
# Sanity check before merging annotation matrixes
head(cbind(rownames(mtx.chrom), rownames(mtx.histone), rownames(mtx.tfbs)))
tail(cbind(rownames(mtx.chrom), rownames(mtx.histone), rownames(mtx.tfbs)))
mtx <- cbind(mtx.chrom, mtx.histone[rownames(mtx.chrom), ], mtx.tfbs[rownames(mtx.chrom), ])
# Unembedding multiple genes
mtx.un <- unembed(data.frame(mtx, genes=rownames(mtx)), "genes", sep=".")
rownames(mtx.un) <- NULL
# Collapling duplicates
mtx.un.agg <- aggregate(data.frame(mtx.un[, !grepl("genes", colnames(mtx.un))]), by=list(mtx.un$genes), min)
mtx.un.agg <- mtx.un.agg[ !grepl("^[0-9]", mtx.un.agg$Group.1, perl=T), ] # Remove numerical IDs
length(mtx.un.agg$Group.1) # Check that all names are unique
length(unique(mtx.un.agg$Group.1))
rownames(mtx.un.agg) <- mtx.un.agg$Group.1
mtx.un.agg$Group.1 <- NULL
# Merging stats matrix, check that the column names match
mtx.stats <- data.frame(cbind(mtx.chrom.stats, mtx.histone.stats, mtx.tfbs.stats))
setdiff(colnames(mtx.stats), colnames(mtx.un.agg))
setdiff(colnames(mtx.un.agg), colnames(mtx.stats))
mtx.stats <- mtx.stats[, colnames(mtx.un.agg)]
# Save merged matrixes
write.table(mtx.un.agg, "results/all_annot_merged.txt", sep="\t", quote=F, col.names=NA)
write.table(as.matrix(mtx.stats), "results/all_stats_merged.txt", sep="\t", quote=F, col.names=NA)
# Define a vector to store summary statsitics
sum.stat <- numeric(length=nrow(mtx.un.agg))
names(sum.stat) <- rownames(mtx.un.agg)
# Calculate summary statistics
# Set non-annotated AND negatively significant cells to "-1", to emphasize the significance of the depleted and absent annotations 
for (i in (1:nrow(mtx.un.agg))) {
  tmp <- mtx.un.agg[i, ]
  tmp[((mtx.un.agg[i,] == 0) & (mtx.stats < 0))] <- -1
  sum.stat[i] <- sum(tmp*mtx.stats)
}
# Check the results
head(sum.stat[order(sum.stat, decreasing=T)])
# Save them
write.table(data.frame(sum.stat[order(sum.stat, decreasing=T)]), "results/annot.txt", sep="\t", quote=F, col.names=NA)
