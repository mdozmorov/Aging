library(dplyr)
library(methylumi)
library(lumi)
library(xlsx)
library(ggplot2)
library(ggrepel)
library(wateRmelon)
library(betareg)
library(pander)
library(sva)
# # https://github.com/mdozmorov/450K_DataProcessing
# source("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/450K_DataProcessing/SRC/detectionPval.filter.R")
source("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/450K_DataProcessing/SRC/filterXY.R")

# biocLite(c('IlluminaHumanMethylation450kmanifest'))
# biocLite('wateRmelon')
                 
# # A modification of the normalization function
# normalizeMethyLumiSet.modified <- function (x, beta.cuts = c(0.2, 0.8), mapfun = c("atan", "ratio")) 
# {
#   if (length(annotation(x)) > 0) {
#     if (annotation(x) == "IlluminaHumanMethylation450k") {
#       message("Normalizing via Illumina controls...")
#       return(normalizeViaControls(x))
#     }
#     if (annotation(x) == "IlluminaHumanMethylation27k") {
#       message("HumanMethylation27 data encountered, skipping...")
#       return(x)
#     }
#   }
#   mapfun = match.arg(mapfun)
#   history.submitted <- as.character(Sys.time())
#   good <- rep(TRUE, ncol(x))
#   cy3 <- unmethylated(x)
#   cy3[cy3 < 0] <- NA
#   cy5 <- methylated(x)
#   cy5[cy5 < 0] <- NA
#   for (i in 1:ncol(cy5)) {
#     cy3inc <- (!is.na(betas(x)[, i]) & !is.na(cy3[, i]))
#     cy5inc <- (!is.na(betas(x)[, i]) & !is.na(cy5[, i]))
#     cy3vec <- cy3[cy3inc, i]
#     cy5vec <- cy5[cy5inc, i]
#     cy3h <- median(cy3vec[betas(x)[cy3inc, i] < beta.cuts[1]])
#     cy3l <- median(cy3vec[betas(x)[cy3inc, i] > beta.cuts[2]])
#     cy5l <- median(cy5vec[betas(x)[cy5inc, i] < beta.cuts[1]])
#     cy5h <- median(cy5vec[betas(x)[cy5inc, i] > beta.cuts[2]])
#     corfactor <- (cy3h - cy3l)/(cy5h - cy5l)
#     cy5[, i] <- cy5[, i] * (corfactor)
#     cy5vec <- cy5[cy5inc, i]
#     newcy5l <- median(cy5vec[betas(x)[cy5inc, i] < beta.cuts[1]])
#     if (newcy5l < cy3l) {
#       cy5[, i] <- cy5[, i] + (cy3l - newcy5l)
#     }
#     else {
#       cy3[, i] <- cy3[, i] + (newcy5l - cy3l)
#     }
#     if (corfactor < 0) {
#       good[i] <- FALSE
#       warning(sprintf("Sample %d has medians that do not make sense for a normal sample\n(cy3l=%f ,cy5l=%f ,cy3h=%f ,cy5h=%f)\nRemoving sample!  Check quality control.", 
#                       i, cy3l, cy5l, cy3h, cy5h))
#       cy5[, i] <- NA
#       cy3[, i] <- NA
#     }
#   }
#   newbeta <- 0
#   if (mapfun == "atan") {
#     newbeta <- atan((cy5)/(cy3))/(pi/2)
#   }
#   else {
#     newbeta <- cy5/(cy5 + cy3 + 100)
#   }
#   assaydata <- new.env(hash = TRUE, parent = emptyenv())
#   assaydata[["unmethylated"]] <- cy3
#   assaydata[["methylated"]] <- cy5
#   assaydata[["betas"]] <- newbeta
#   history.finished <- as.character(Sys.time())
#   history.command <- capture.output(print(match.call(normalizeMethyLumiSet)))
#   ret <- new("MethyLumiSet", phenoData = phenoData(x), featureData = featureData(x), 
#              assayData = assaydata, annotation = annotation(x))
#   # QCdata(ret) <- QCdata(x)
#   ret <- ret[, good]
#   ret@history <- rbind(getHistory(x), data.frame(submitted = history.submitted, 
#                                                  finished = history.finished, command = history.command))
#   return(ret)
# }
# 
# A function to pull out p-value of LM. https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# # Read in 450K annotation
# ANNOT <- "/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/GPL13534-11288.txt.gz"
# annot <- read.table("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/GPL13534-11288.txt.gz", sep = "\t", skip=37, header=T, row.names=1, as.is = T, fill = T)
# colnames(annot)
# annot.slim <- annot[, c("CHR", "RANGE_START", "RANGE_END", "Name", "UCSC_RefGene_Name")]

# Read in methylation values
MethyFileName <- "data/subset_1000.csv.gz"
# MethyFileName <- "/Users/mikhail/Documents/Data/Work_OMRF/Amr/Control Naive CD4+ T cell Methylation/Control Naive CD4 Tcell Methylation Values.csv.gz"
mldat <- methylumiR(filename = MethyFileName)
# save(list = "mldat", file = "data/mldat.rda")
# load("/Users/mikhail/Documents/Data/Work_OMRF/Amr/Control Naive CD4+ T cell Methylation/mldat.rda")
# mldat

# Read in pheno data, and Match order of samples
PhenoFileName <- "data/Female SLE Naive CD4 T cell Control Info.xlsx"
samp <- read.xlsx2(file = PhenoFileName, sheetName = "Sheet1", stringsAsFactors = F)
samp$Age <- as.numeric(samp$Age) # Make age numeric
samp <- mutate(samp, names = paste(BeadChip.ID, Chip.Placement, sep = "_"))
samp <- samp[match(sampleNames(mldat), samp$names), ]
# Check names
all.equal(sampleNames(mldat), samp$names) # Should be true
rownames(samp) <- samp$names
# Add pheno data
pData(mldat) <- samp
# mldat

print(paste("Total number of CpG probes to start with:", nrow(mldat)))
# Batches
batch.BeadChip.ID <- table(pData(mldat)$BeadChip.ID) %>% sort # Summary how many samples per batch
batch.BeadChip.ID.one <- names(batch.BeadChip.ID)[ batch.BeadChip.ID == 1 ] # Sample names associated with one-sample batch
batch.BeadChip.ID.three <- names(batch.BeadChip.ID)[ batch.BeadChip.ID == 3 ] # Sample names associated with three-sample batch
batch.BeadChip.ID.four <- names(batch.BeadChip.ID)[ batch.BeadChip.ID == 4 ] # Etc.
batch.BeadChip.ID.six <- names(batch.BeadChip.ID)[ batch.BeadChip.ID == 6 ]
batch.BeadChip.ID.twelve <- names(batch.BeadChip.ID)[ batch.BeadChip.ID == 12 ]

# Removing CpG cites having NAs in batches
ind.one <- pData(mldat)$BeadChip.ID %in% batch.BeadChip.ID.one # TRUE/FALSE vector of samples associated with one-sample batch
ind.three <- pData(mldat)$BeadChip.ID %in% batch.BeadChip.ID.three # TRUE/FALSE vector of samples associated with three-sample batch
ind.four <- pData(mldat)$BeadChip.ID %in% batch.BeadChip.ID.four # Etc
ind.six <- pData(mldat)$BeadChip.ID %in% batch.BeadChip.ID.six
ind.twelve <- pData(mldat)$BeadChip.ID %in% batch.BeadChip.ID.twelve
ind.exclude <- vector(mode = "numeric") # Row numbers to exclude
for (i in 1:nrow(mldat)) {
  # If a row (CpG probe) contains NAs in samples associated with one- or three- or ... -sample batches, keep its index 
  if (any(sapply(batch.BeadChip.ID.one, function(x) sum(is.na(betas(mldat)[i, ind.one])) == 1 )) |
      any(sapply(batch.BeadChip.ID.three, function(x) sum(is.na(betas(mldat)[i, ind.three])) == 3 )) |
      any(sapply(batch.BeadChip.ID.three, function(x) sum(is.na(betas(mldat)[i, ind.four])) == 4 )) |
      any(sapply(batch.BeadChip.ID.three, function(x) sum(is.na(betas(mldat)[i, ind.six])) == 6 )) |
      any(sapply(batch.BeadChip.ID.three, function(x) sum(is.na(betas(mldat)[i, ind.twelve])) == 12 ))
      ) {
    ind.exclude <- c(ind.exclude, i)
  }
}
mldat <- mldat[-ind.exclude, ] # Remove offenting rows
print(paste("Total number of CpG probes after removing CpG probes with NAs confounding batch:", nrow(mldat)))

# # Check average p-values, should be similar
# avgPval <- colMeans(pvals(mldat))
# par(las=2, oma=c(5,0,0,0))
# barplot(avgPval,ylab="Average P-Value")

# # Cluster data as-is
# IAC <- cor(betas(mldat)[ complete.cases(betas(mldat)) & apply(betas(mldat), 1, sd) != 0, ]) # Remove rows with NAs & Precaution against zero SD
# cluster1 <- hclust(as.dist(1-IAC), method="ward.D") # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
# plot(cluster1, cex=0.7)#, xlab = mRNA.annot$Class)

# PCA: Check for batch effects. Select one batch, to color points by its assignment
# batch <- "Race"
batch <- "BeadChip.ID"
# batch <- "Chip.Placement"
pca <- betas(mldat)[ complete.cases(betas(mldat)) & apply(betas(mldat), 1, sd) != 0, ] %>% scale %>% t %>% prcomp
data.frame(summary(pca)$importance)[, 1:5] %>% pander # Percent of variance explained
# What correlates with the first component
pca.lm <- lm( PC1 ~ factor(BeadChip.ID), data = cbind(samp, pca$x))
print(paste("BeadChip.ID accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the first principle component, p-value", signif(lmp(pca.lm), 5)))
pca.lm <- lm( PC1 ~ factor(Chip.Placement), data = cbind(samp, pca$x))
print(paste("Chip.Placement accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the first principle component, p-value", signif(lmp(pca.lm), 5)))
pca.lm <- lm( PC1 ~ factor(Race), data = cbind(samp, pca$x))
print(paste("Race accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the first principle component, p-value", signif(lmp(pca.lm), 5)))

pt <- ggplot(data = data.frame(pca$x, pData(mldat), samples = pData(mldat)$Sample.ID, stringsAsFactors = F), 
             aes(x = PC1, y = PC2, label = eval(parse(text = batch)))) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("Batch coloring:", batch)) +
  geom_point(aes(color = eval(parse(text = batch))), size = 3) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text_repel(colour = "black", size = 3) +
  labs(color = batch)
plot(pt)


# # Filtering
print(paste("Number of probes before filtering:", nrow(mldat)))
# # Exclude CpG probes with SNPs
# snpsites <- read.table("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/snpsites.txt.gz", sep="\t", stringsAsFactors = FALSE)
# ind <- which(is.element(featureNames(mldat), snpsites$V1))
# mldat <- mldat[-ind, ]
# Exclude non-specific probes
annot.nonspecific <- openxlsx::read.xlsx("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/48639-non-specific-probes-Illumina450k.xlsx", sheet = 1, cols = 1)
ind <- which(is.element(featureNames(mldat), annot.nonspecific$TargetID))
mldat <- mldat[-ind, ]
print(paste("After removing", length(ind), "non-specific probes,", nrow(mldat), "probes remain."))
# Exclude polymorphic probes
annot.polymorphic <- openxlsx::read.xlsx("/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/48640-polymorphic-CpGs-Illumina450k.xlsx", sheet = 1, cols = 1)
ind <- which(is.element(featureNames(mldat), annot.polymorphic$PROBE))
mldat <- mldat[-ind, ]
print(paste("After removing", length(ind), " polymorphic probes,", nrow(mldat), "probes remain."))
# # Filter probes showing detection p-value > 0.05 in at least 10% of the samples
# mldat <- detectionPval.filter2(mldat, detectionPval.threshold=0.05, detectionPval.perc.threshold=10, projectName = NULL, PATH="./")
# print(paste("After removing probes showing detection p-value > 0.05 in at least 10% of the samples,", nrow(mldat), "probes remain."))
# Filter chrY probes
mldat <- filterXY(mldat, chr_remove = c("Y"))
print(paste("After removing chrY probes,", nrow(mldat), "probes remain."))
# Remove probes with median betas = 0% or 100% across all samples
ind <- rowMedians(betas(mldat), na.rm = T)
mldat <- mldat[ ind != 1 | ind != 0, ]
print(paste("After removing probes with median betas = 0% or 100% across all samples,", nrow(mldat), "probes remain."))

# Normalization 
# quality filter using default thresholds
# melon.pf <- pfilter(mldat)
# # preprocess using our best method
# melon.dasen.pf <- dasen(melon.pf)
# melon.dasen.pf <- BMIQ(melon.dasen.pf)
# save(list=c("melon.dasen.pf"), file = "data/melon_dasen_pf_BMIQ.rda")
# load(file = "/Users/mikhail/Documents/Data/Work_OMRF/Amr/Control Naive CD4+ T cell Methylation/melon_dasen_pf_BMIQ.rda")

# Convert to lumi object
# mldat.lumi <- as(mldat, "MethyLumiM")
# # Diagnostic plots from the tutorial
# plotSampleRelation(mldat.lumi, method='mds', cv.Th=0)
# plotSampleRelation(mldat.lumi, method='cluster', cv.Th=0)
# density(mldat.lumi, xlab="M-value")
# density(estimateIntensity(mldat.lumi), xlab="log2(CpG-site Intensity)")
# boxplot(mldat.lumi)
# boxplot(estimateIntensity(mldat.lumi))
# colorBiasSummary(mldat.lumi[,1:8], channel='methy')
# plotColorBias1D(mldat.lumi)
# plotColorBias1D(mldat.lumi, channel='methy')
# plotColorBias1D(mldat.lumi, channel='unmethy')
# plotColorBias1D(mldat.lumi, channel='sum')
# boxplotColorBias(mldat.lumi, channel='methy')
# boxplotColorBias(mldat.lumi, channel='sum')

# mldat.lumi.coloradj <- lumiMethyC(mldat.lumi, method = 'quantile')
# mldat.lumi.quantile <- lumiMethyN(mldat.lumi.coloradj, method='quantile')
# OR
# mldat.lumi.quantile <- lumiMethyB(mldat.lumi, method="bgAdjust2C", separateColor=TRUE)

# # Adjust for batch effect, using ComBat, select one
# ind <- c("9611519008", "9611519009", "9974366123") # Batches with single samples, should be removed
# melon.dasen.pf <- melon.dasen.pf[, !(pData(melon.dasen.pf)$BeadChip.ID %in% ind)] # Remove them
# batch <- pData(melon.dasen.pf)$BeadChip.ID # Batch
# modcombat <- model.matrix(~1 , data=pData(melon.dasen.pf)) # No effect
# combat_edata <- ComBat(betas(melon.dasen.pf), batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
# save(list = c("combat_edata"), file = "data/combat_edata_parpriorTrue.rda")
load(file = "data/combat_edata_parpriorTrue.rda")

# Adjusting for batch using swamp::combat
# combat_edata <- combat(betas(mldat), data.frame(Factor1 = factor(samp$BeadChip.ID), Factor2 = factor(samp$Chip.Placement), Factor3 = factor(samp$Race), row.names = samp$names), batchcolumn = 3)

# # Adjust for three batch effects. Modified after limma::removeBatchEffect
# # Batch 1
# batch <- as.factor(pData(melon.dasen.pf)$Race)
# contrasts(batch) <- contr.sum(levels(batch))
# batch <- model.matrix(~batch)[, -1, drop = FALSE]
# # Batch 2
# batch2 <- as.factor(pData(melon.dasen.pf)$BeadChip.ID)
# contrasts(batch2) <- contr.sum(levels(batch2))
# batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
# # Batch 3
# batch3 <- as.factor(pData(melon.dasen.pf)$Chip.Placement)
# contrasts(batch3) <- contr.sum(levels(batch3))
# batch3 <- model.matrix(~batch3)[, -1, drop = FALSE]
# # Combine them
# X.batch <-cbind(batch, batch2, batch3)
# design = matrix(1, ncol(betas(melon.dasen.pf)), 1)
# # Estimate predictions using batch effects, and subtract them
# fit <- limma::lmFit(betas(melon.dasen.pf), cbind(design, X.batch), na.action=na.omit)
# beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
# beta[is.na(beta)] <- 0
# combat_edata <- as.matrix(betas(melon.dasen.pf)) - beta %*% t(X.batch)

# Sanity check: Comparing correlations before/after adjustments
cor.raw <- cor(pData(melon.dasen.pf)$Age, t(betas(melon.dasen.pf))) # Correlation with age before any adjustments
cor.adj <- cor(pData(melon.dasen.pf)$Age, t(combat_edata)) # Correlation with age after batch adjustment
cor(as.numeric(cor.raw), as.numeric(cor.adj), use = "complete.obs") # How raw and adjusted vectors of correlation coefficients compare?

# Correlation analysis
mldat.final = new(Class = "MethyLumiSet")
betas(mldat.final) <- combat_edata
pData(mldat.final) <- pData(melon.dasen.pf)

# mldat.final <- melon.dasen.pf

# Allocate storage
coeffs <- vector(mode = "list", length = nrow(mldat.final)) # Store coefficients for regular regression
adj.r.squared <- vector(mode = "list", length = nrow(mldat.final)) # R2 for regular regression
lmpval <- vector(mode = "list", length = nrow(mldat.final)) # p-value for regular regression
coeffs.M <- vector(mode = "list", length = nrow(mldat.final)) # Store coefficients for regular regression on M values
adj.r.squared.M <- vector(mode = "list", length = nrow(mldat.final)) # R2 for regular regression on M values
lmpval.M <- vector(mode = "list", length = nrow(mldat.final)) # p-value for regular regression on M values
betareg.coeffs <- vector(mode = "list", length = nrow(mldat.final)) # Store coefficients for betareg regression
betareg.adj.r.squared <- vector(mode = "list", length = nrow(mldat.final)) # R2 for betareg regression
correl <- vector(mode = "list", length = nrow(mldat.final)) # Pearson correlation coefficients

# Test each row (CpG expression profile) for association with Age
for (i in 1:nrow(mldat.final)) {
  print(i)
  y <- betas(mldat.final)[i, ] # Vector of beta values to test
  # Regular regression
  fit <- lm(y ~ Age, data = pData(mldat.final), na.action=na.omit) # Age only
  #  fit2 <- lm(betas(mldat.lumi.quantile.final)[i, ] ~ Age + Race + BeadChip.ID + Chip.Placement, data = pData(mldat.lumi.quantile.final)) # Age plus batches
  coeffs[[i]] <- c(summary(fit)$coefficients[2, "Estimate"], summary(fit)$coefficients[2, "Pr(>|t|)"]) #, summary(fit2)$coefficients[2, "Estimate"], summary(fit2)$coefficients[2, "Pr(>|t|)"]) # The first coefficient, ignore intercept
  adj.r.squared[[i]] <- c(summary(fit)$adj.r.squared) #, summary(fit2)$adj.r.squared)
  lmpval[[i]] <- c(lmp(fit)) #, lmp(fit2))
  # Regular regression on M values
  fit.M <- lm(beta2m(y) ~ Age, data = pData(mldat.final), na.action=na.omit) # Age only
  #  fit2 <- lm(betas(mldat.lumi.quantile.final)[i, ] ~ Age + Race + BeadChip.ID + Chip.Placement, data = pData(mldat.lumi.quantile.final)) # Age plus batches
  coeffs.M[[i]] <- c(summary(fit.M)$coefficients[2, "Estimate"], summary(fit.M)$coefficients[2, "Pr(>|t|)"]) #, summary(fit2)$coefficients[2, "Estimate"], summary(fit2)$coefficients[2, "Pr(>|t|)"]) # The first coefficient, ignore intercept
  adj.r.squared.M[[i]] <- c(summary(fit.M)$adj.r.squared) #, summary(fit2)$adj.r.squared)
  lmpval.M[[i]] <- c(lmp(fit.M)) #, lmp(fit2))
  # betareg regression
  y[y > 1] <- 0.9999999999999999 # Set large outliers to maximum of (0, 1) range
  y[y < 0] <- 0.0000000000000001 # Set small outliers to minimum of (0, 1) range
  betareg.fit <- betareg(y ~ Age, data = pData(mldat.final), na.action=na.omit)
  betareg.coeffs[[i]] <- c(summary(betareg.fit)$coefficients$mean[2, "Estimate"], summary(betareg.fit)$coefficients$mean[2, "Pr(>|z|)"]) #, summary(fit2)$coefficients[2, "Estimate"], summary(fit2)$coefficients[2, "Pr(>|t|)"]) # The first coefficient, ignore intercept
  betareg.adj.r.squared[[i]] <- c(summary(betareg.fit)$pseudo.r.squared) #, summary(fit2)$adj.r.squared)
  # Correlation
  corr <- Hmisc::rcorr(betas(mldat.final)[i, ], Biobase::pData(mldat.final)$Age)
  correl[[i]] <- c(corr$r[1, 2], corr$P[1, 2])
  if (i %% 100000 == 0) { save(list = c("fit", "coeffs", "adj.r.squared", "lmpval", "correl"), file = paste0("results/lm_results_", i, ".rda")) }
}

# Combine results
results <- cbind(do.call(cbind, lapply(coeffs, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(adj.r.squared, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(lmpval, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(coeffs.M, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(adj.r.squared.M, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(lmpval.M, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(betareg.coeffs, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(betareg.adj.r.squared, data.frame, stringsAsFactors=F)) %>% t,
                 do.call(cbind, lapply(correl, data.frame, stringsAsFactors=F)) %>% t)
rownames(results) <- rownames(betas(mldat.final))
colnames(results) <- c("coeff.lm", "pval.lm", "r2.lm", "lmpval", "coeff.lm.M", "pval.lm.M", "r2.lm.M", "lmpval.M", "coeff.br", "pval.br", "r2.br", "corr.rcorr", "pval.rcorr")
write.table(results, "results/correlation_results.txt", sep = "\t", quote = FALSE, col.names = NA)
# 
# # # Diagnostic plots, comparing histograms of the two distributions. Select one factor set
# factor1 <- "coeff.lm"; factor2 <- "coeff.br"
# factor1 <- "r2.lm"; factor2 <- "r2.br"
# factor1 <- "pval.lm"; factor2 <- "pval.br"
# # Prepare data, -log10-transform p-values only
# mtx.plot <- melt(results[, c(factor1, factor2)])
# mtx.plot$value <- -log10(mtx.plot$value)
# # Plot histograms
# mtx.plot %>% ggplot(aes(x=value)) +
#   geom_histogram(data=subset(mtx.plot, Var2 == factor1), fill = "red", alpha = 0.2) +
#   geom_histogram(data=subset(mtx.plot, Var2 == factor2), fill = "blue", alpha = 0.2)
# 
