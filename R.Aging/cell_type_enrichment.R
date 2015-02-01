# cell.sig <- function(fname) {
#   mtx.count <- read.table(fname, sep="\t", header=F, row.names=1, stringsAsFactors=F)
#   sig <- sum(mtx.count$V2) # Count of all significant cell lines
#   tot <- sum(mtx.count$V3) # Count of all tested cell linesi
#   p <- list() # List to store p-values
#   for (i in 1:nrow(mtx.count)) { 
#     cell.2x2 <- matrix(c(mtx.count$V2[i], mtx.count$V3[i] - mtx.count$V2[i],
#                          sig - mtx.count$V2[i], tot - sig - (mtx.count$V3[i] - mtx.count$V2[i])),
#                        ncol=2, dimnames=list(c("sig", "not sig"), c("selected", "not selected"))) 
#     p <- c(p, list(fisher.test(cell.2x2)$p.value))
#   }
#   
#   tmp <- cbind(mtx.count, formatC(unlist(p), digits=2, format="e"))
#   colnames(tmp) <- c(sig", "tot", "pval")
#   tmp <- tmp[order(as.numeric(tmp[, "pval"])), ]
#   unlink(fname)
#   write.table(tmp, fname, sep="\t", col.names=NA, quote=F)
# }
# 
# # Run once, will overwrite the original matrix
# cell.sig("../GenomeRunner/Hannum/markers_vs_encBroadHmm//cell_counts.txt")
# cell.sig("../GenomeRunner/Hannum/markers_vs_encHistone//cell_counts.txt")
# cell.sig("../GenomeRunner/Hannum/markers_vs_encTfbs//cell_counts.txt")
# cell.sig("../GenomeRunner/Hannum/tissue_vs_encBroadHmm//cell_counts.txt")
# cell.sig("../GenomeRunner/Hannum/tissue_vs_encHistone//cell_counts.txt")
# cell.sig("../GenomeRunner/Hannum/tissue_vs_encTfbs//cell_counts.txt")

# All counts
fname <- "/Users/mikhail/Documents/Work/GenomeRunner/Aging/GenomeRunner/cell_counts_all_10.txt"
mtx.count <- read.table(fname, sep="\t", header=F, stringsAsFactors=F)
mtx.count <- ddply(mtx.count, "V1", summarise, sum(V2), sum(V3))
colnames(mtx.count) <- c("V1", "V2", "V3")
sig <- sum(mtx.count$V2) # Count of all significant cell lines
tot <- sum(mtx.count$V3) # Count of all tested cell lines
cell.2x2 <- matrix(c(mtx.count$V2[i], mtx.count$V3[i] - mtx.count$V2[i],
sig - mtx.count$V2[i], tot - sig - (mtx.count$V3[i] - mtx.count$V2[i])),
ncol=2, dimnames=list(c("sig", "not sig"), c("selected", "not selected")))
p <- c(p, list(fisher.test(cell.2x2)$p.value))
sig <- sum(mtx.count$V2) # Count of all significant cell lines
tot <- sum(mtx.count$V3) # Count of all tested cell linesi
p <- list() # List to store p-values
for (i in 1:nrow(mtx.count)) {
cell.2x2 <- matrix(c(mtx.count$V2[i], mtx.count$V3[i] - mtx.count$V2[i],
sig - mtx.count$V2[i], tot - sig - (mtx.count$V3[i] - mtx.count$V2[i])),
ncol=2, dimnames=list(c("sig", "not sig"), c("selected", "not selected")))
p <- c(p, list(fisher.test(cell.2x2)$p.value))
}
tmp <- cbind(mtx.count, formatC(unlist(p), digits=2, format="e"))
colnames(tmp) <- c("cell", "sig", "tot", "pval")
tmp <- tmp[order(as.numeric(tmp[, "pval"])), ]
tmp$cell <- toupper(tmp$cell)
cellAnnot$cell <- toupper(cellAnnot$cell)
tmp <- left_join(tmp, cellAnnot, by= c("cell" = "cell"))
write.table(tmp[tmp$cell != "", ], "/Users/mikhail/Documents/Work/GenomeRunner/Aging/GenomeRunner/cell_counts_all_pval_10.txt", sep="\t", row.names=F, quote=F)

