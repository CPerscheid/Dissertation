#! /usr/local/bin/Rscript
library(genefilter)
library(GEOquery)
library(zFPKM)
library(PCAtools)


home.folder <- "/path/to/working/directory/"

file.prefix <- "BRCA/SCANBII"
accession.number <- "GSE81538"

data <- getGEO(accession.number, GSEMatrix=TRUE)
#get expression data
filePaths = getGEOSuppFiles(accession.number)
expr.filename <- row.names(filePaths)[2]
expr <- read.csv(expr.filename, row.names = 1)
#add sample IDs as colnames (we assume that the order remains the same)
colnames(expr) <- colnames(data[[1]]@assayData[["exprs"]])
  
#metadata
metadata <- data[[1]]@phenoData@data
  
#remove NA samples
expr <- expr[, colSums(is.na(expr)) != nrow(expr)]
  
#remove unclassified samples
conditioned <- row.names(metadata[!is.na(metadata[["pam50 subtype:ch1"]]),])
intersection <- intersect(colnames(expr), conditioned)
expr <- expr[,intersection]
metadata <- metadata[intersection,]
  
#reverse-transform from log2-FPKM and subtract 0.1 to have original FPKM values
#excerpt from study: "...adding to each expression measurement 0.1 FPKM, performing a log2 transformation."
exp.fpkm <- 2^expr
exp.fpkm.original <- exp.fpkm - 0.1
exp.zfpkm <- zFPKM(exp.fpkm.original)

#pp <- zFPKMPlot(exp.fpkm.original)
p <- pca(exp.zfpkm, metadata = metadata)
pairsplot(p, colby="pam50 subtype:ch1")
exp.zfpkm.nooutliers <- exp.zfpkm[,(p$rotated$PC1 > -800)]
m2 <- metadata[(p$rotated$PC1 > -800),]
p2 <- pca(exp.zfpkm.nooutliers, metadata = m2)
pairsplot(p2, colby="pam50 subtype:ch1")

#match sample ids for expression and metadata again
intersection <- intersect(colnames(exp.zfpkm.nooutliers), conditioned)
exp.zfpkm.nooutliers <- exp.zfpkm.nooutliers[,intersection]
metadata <- metadata[intersection,]


p2 <- pca(exp.zfpkm.nooutliers, metadata = metadata)
pairsplot(p2, colby="pam50 subtype:ch1")

  
#filter out lowly expressed genes via genefilter
thres <- (ncol(exp.zfpkm.nooutliers) * 30) / 100
#filter all expression values that have zfpkm score above 3.0 in more than 30% of the samples
yy <- rowMedians(as.matrix(exp.zfpkm.nooutliers))
expr.zfpkm.nooutliers.filtered <- exp.zfpkm.nooutliers[(yy > -3.0), ]


#put genes into columns
expr.zfpkm.nooutliers.filtered <- t(expr.zfpkm.nooutliers.filtered)
  
expr.prep.filename <- paste0(home.folder, file.prefix, "_normalized_nonormals_expressions.csv")
  
#save with genes in columns
metadata <- metadata[,c("title","geo_accession", "pam50 subtype:ch1")]
colnames(metadata) <- c("title", "geo_accession", "pam50_subtype:ch1")
write.table(expr.zfpkm.nooutliers.filtered, file= expr.prep.filename, sep = ";")

#write metadata with samples in columns
write.table(t(metadata), file=paste0(home.folder, file.prefix, "_metadata.csv"), sep = ";")



