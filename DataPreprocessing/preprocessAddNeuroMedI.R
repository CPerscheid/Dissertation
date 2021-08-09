#! /usr/local/bin/Rscript
library(SummarizedExperiment)#for colData function
library(DESeq2)
library(genefilter)
library(GEOquery)
library(PCAtools)

preprocessGEOData <- function(accession.number, home.folder, fileprefix){
  data <- getGEO(accession.number, GSEMatrix=TRUE)
  
  #metadata
  metadata <- data[[1]]@phenoData@data
  #expression data
  expr <- data[[1]]@assayData[["exprs"]]
  
  expr.filename <- paste0(home.folder, fileprefix, "_normalized_expressions.csv")
  
  #remove normal/unknown/unclassified samples
  conditioned <- row.names(metadata[!metadata[["status:ch1"]] %in% c("CTL to AD", "MCI to CTL", "OTHER", "borderline MCI"),])
  intersection <- intersect(colnames(expr), conditioned)
  expr <- expr[,intersection]
  
  #filter out lowly expressed genes via genefilter
  thres <- (ncol(expr) * 70) / 100
  #filter all expression values that have a value below 6.0 in more than 30% of the samples
  expr <- expr[(rowSums(expr > 6.0)) > thres, ]
  #p <- pca(expr, metadata = metadata)
  #pca.pl <- pairsplot(p, colby="status:ch1")
  
  
  #put genes into columns
  expr <- t(expr)
  
  
  #save with genes in columns
  write.csv(expr, file= expr.filename)
  
  #write metadata with samples in columns
  write.csv(t(metadata), file=paste0(home.folder, fileprefix, "_metadata.csv"))
  
  
}

home.folder <- "/path/to/working/directory/"
preprocessGEOData('GSE63060', home.folder, "Alzheimer/BatchI/AddNeuroMed_BatchI")
