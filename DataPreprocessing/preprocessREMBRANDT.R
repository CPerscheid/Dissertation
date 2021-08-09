#! /usr/local/bin/Rscript
library(genefilter)
library(GEOquery)

home.folder <- "/Users/Cindy/Studium/PhD/ResearchProjects/Diss/AnalysisFramework/Case_Studies/"

data <- getGEO("GSE108474", GSEMatrix=TRUE)

#metadata
metadata <- data[[1]]@phenoData@data
#expression data
expr <- data[[1]]@assayData[["exprs"]]

#remove NA samples
expr <- expr[, colSums(is.na(expr)) != nrow(expr)]

#remove normal samples
conditioned <- row.names(metadata[!metadata[["disease:ch1"]] %in% c("normal", "unknown", "unclassified", "mixed"),])
intersection <- intersect(colnames(expr), conditioned)
expr <- expr[,intersection]
metadata <- metadata[intersection,]

#plot
p <- pca(expr, metadata = metadata)
pairsplot(p, colby="disease:ch1")

#filter out lowly expressed genes via genefilter
thres <- (ncol(expr) * 30) / 100
#filter all expression values that are below 0 (neg. log2 values mean a count betw. 0 and 1)
#TODO: set threshold accordingly
expr <- expr[(rowSums(expr > 7.5)) > thres, ]
  
expr.filename <- paste0(home.folder, "Glioma/REMBRANDT/REMBRANDT", "_normalized_filtered_nonormals_expressions.csv")


#save with genes in columns
write.csv(t(expr), file= expr.filename)

#write metadata with samples in columns
write.csv(t(metadata), file=paste0(home.folder, "Glioma/REMBRANDT/REMBRANDT", "_metadata.csv"))
