fileNameIn <- commandArgs(TRUE)[1]
mtx <- openxlsx::read.xlsx(fileNameIn, sheet = 1, colNames = FALSE)
fileNameOut <- sub("xlsx", "txt", fileNameIn)
write.table(mtx, fileNameOut, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
