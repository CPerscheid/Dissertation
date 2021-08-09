library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)


home.folder <- "/Users/Cindy/Studium/PhD/ResearchProjects/Diss/Case_Studies/"
########### BRCA ##########
brca.project <- c("TCGA-BRCA")
brca.subtypeColumn = "subtype_BRCA_Subtype_PAM50"
brca.fileprefix = paste0(home.folder, "BRCA/BRCA_TP_")
brca.metadata.attributes = c("sample","patient","barcode","shortLetterCode","definition","tumor_stage","classification_of_tumor","tumor_grade","primary_diagnosis","age_at_diagnosis","gender","disease_type","primary_site","project_id","subtype_Tumor.Type","subtype_Included_in_previous_marker_papers","subtype_Tumor_Grade","subtype_BRCA_Pathology","subtype_BRCA_Subtype_PAM50","subtype_MSI_status","subtype_CNV.Clusters","subtype_Mutation.Clusters","subtype_DNA.Methylation.Clusters","subtype_mRNA.Clusters","subtype_miRNA.Clusters","subtype_lncRNA.Clusters","subtype_Protein.Clusters")

print("query GDC")
#get only primary tumor samples 
TCGAquery <- GDCquery(project= brca.project, data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow = "HTSeq - Counts", barcode = c(""), sample.type = c("Primary solid Tumor"))
  
#download the actual data
GDCdownload(TCGAquery)
#prepares the data for analysis and putting it into an SummarizedExperiment object
expr1 <- GDCprepare(TCGAquery)
print("finished...")

print("remove NAs")
#get the metadata
metadata <- colData(expr1)
expr.subtype <- metadata[,c("sample", brca.subtypeColumn)]
#remove NAs
y <- is.na(expr.subtype[[brca.subtypeColumn]])
has.subtype <- expr.subtype[!is.na(expr.subtype[[brca.subtypeColumn]]),]
#subset only to the patients with available subtype information
expr.subset <- subset(expr1, select = colData(expr1)$sample %in% has.subtype$sample)

#remove unclassified samples
conditioned <- row.names(metadata[!is.na(metadata[[brca.subtypeColumn]]),])
intersection <- intersect(colnames(expr.subset), conditioned)
expr <- expr.subset[,intersection]
  
print("finished...")

print("create DGE object")
#create DGEList object
genes <-rowRanges(expr)$external_gene_name
expr.subtype <- colData(expr)[[brca.subtypeColumn]]
expr.counts <- assay(expr) 
expr.dge <- DGEList(expr.counts, group=expr.subtype, genes=genes)
print("finished...")

print("filter dataset...")

design <- model.matrix(~ factor(expr.subtype))
#voom(expr.dge, design, plot=TRUE)

keep <- filterByExpr(expr.dge, design, min.count = 60)
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
metadata.filtered <- metadata[,brca.metadata.attributes]
write.table(metadata.filtered, file=paste0(brca.fileprefix, "metadata.csv"), sep = ";")
#write counts data
write.table(expr.normalized.matrix, file=paste0(brca.fileprefix, "expressions_normalized.csv"), sep = ";")






