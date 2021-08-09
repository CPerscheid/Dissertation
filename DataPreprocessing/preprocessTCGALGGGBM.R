library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)

home.folder <- "/Users/Cindy/Studium/PhD/ResearchProjects/Diss/AnalysisFramework/Case_Studies/"

gbm.lgg.project <- c("TCGA-LGG", "TCGA-GBM")
gbm.lgg.subtypeColumn = "subtype_Histology"
gbm.lgg.fileprefix = paste0(home.folder, "Glioma/TCGA/GBM-LGG_TP_")
gbm.lgg.metadata.attributes = c("sample","patient","barcode","shortLetterCode","definition","age_at_diagnosis","tumor_stage","primary_diagnosis","tumor_grade","classification_of_tumor","gender","disease_type","project_id","subtype_Tissue.source.site","subtype_Study","subtype_BCR","subtype_RNAseq","subtype_SNP6","subtype_U133a","subtype_HM450","subtype_HM27","subtype_RPPA","subtype_Histology","subtype_Grade","subtype_Mutation.Count","subtype_IDH.status","subtype_IDH.codel.subtype","subtype_MGMT.promoter.status","subtype_Original.Subtype","subtype_Transcriptome.Subtype","subtype_Pan.Glioma.RNA.Expression.Cluster","subtype_IDH.specific.RNA.Expression.Cluster","subtype_Pan.Glioma.DNA.Methylation.Cluster","subtype_IDH.specific.DNA.Methylation.Cluster","subtype_Supervised.DNA.Methylation.Cluster","subtype_Random.Forest.Sturm.Cluster","subtype_RPPA.cluster")

print("query GDC")
#get only primary tumor samples 
TCGAquery <- GDCquery(project= gbm.lgg.project, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow = "HTSeq - Counts", barcode = c(""), sample.type = c("Primary solid Tumor"))

#download the actual data
GDCdownload(TCGAquery)
#prepares the data for analysis and putting it into an SummarizedExperiment object
expr1 <- GDCprepare(TCGAquery)
print("finished...")

print("remove NAs")
#get the metadata
metadata <- colData(expr1)
expr.subtype <- metadata[,c("sample", gbm.lgg.subtypeColumn)]
#remove NAs
y <- is.na(expr.subtype[[gbm.lgg.subtypeColumn]])
has.subtype <- expr.subtype[!is.na(expr.subtype[[gbm.lgg.subtypeColumn]]),]
#subset only to the patients with available subtype information
expr.subset <- subset(expr1, select = colData(expr1)$sample %in% has.subtype$sample)

#remove normal/unknown/unclassified/oligoastrocytoma samples
conditioned <- row.names(metadata[metadata[[gbm.lgg.subtypeColumn]] %in% c("astrocytoma", "oligodendroglioma", "glioblastoma"),])
intersection <- intersect(colnames(expr.subset), conditioned)
expr <- expr.subset[,intersection]

print("finished...")

print("create DGE object")
#create DGEList object
genes <-rowRanges(expr)$external_gene_name
expr.subtype <- colData(expr)[[gbm.lgg.subtypeColumn]]
expr.counts <- assay(expr) 
expr.dge <- DGEList(expr.counts, group=expr.subtype, genes=genes)
print("finished...")

print("filter dataset...")

design <- model.matrix(~ factor(expr.subtype))
voom(expr.dge, design, plot=TRUE)

keep <- filterByExpr(expr.dge, design)
#keep <- filterByExpr(expr.dge, group = expr.subtype)
print("finished...")
expr <- expr.dge[keep,]

print("normalize dataset...")
#normalize counts with TMM
expr.normalized <- calcNormFactors(expr, method = "TMM")
expr.normalized.matrix <- log2(cpm(expr.normalized) + 1)
print("finished...")

print("write dataset...")
#write metadata
metadata <- colData(expr.subset[row.names(expr.normalized.matrix)])
metadata.filtered <- metadata[,gbm.lgg.metadata.attributes]
write.table(metadata.filtered, file=paste0(gbm.lgg.fileprefix, "metadata.csv"), sep = ";")
#write counts data
write.table(expr.normalized.matrix, file=paste0(gbm.lgg.fileprefix, "expressions_normalized.csv"), sep = ";")






