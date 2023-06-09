# script to run survival analysis using TCGA data 
# setwd()

#load packages------
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# get clinical data from TCGA------
clinical_stad <- GDCquery_clinic("TCGA-STAD")
# check out if those key informations are available
any(colnames(clinical_stad)%in%c("vital_status","days_to_last_follow_up","days_to_death"))
# which column
which(colnames(clinical_stad)%in%c("vital_status","days_to_last_follow_up","days_to_death"))
# take a quick look
clinical_stad[,c(9,38,44)]
table(clinical_stad$vital_status)

# delete not reported data------
which(clinical_stad$vital_status == "Not Reported")
clinical_stad[319,]$days_to_last_follow_up
clinical_stad2<- clinical_stad[-319,]
table(clinical_stad2$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patients death
# days_to_last_follow_up, that is the number of days passed from the initial diagnosis to the last visit

# change certain values the way they are encoded------
clinical_stad2$deceased <- ifelse(clinical_stad2$vital_status == "Alive", FALSE, TRUE)
#assign false value to all the patients that are yet alive 
#this new column will be using in the following survival analysis
# have a quick check
table(clinical_stad2$deceased)

# create an "overall survival" variable that is equal to days_to_death------
# for dead patients, and to days_to_last_follow_up for patients who are still alive
clinical_stad2$overall_survival <- ifelse(clinical_stad2$vital_status == "Alive",
                                          clinical_stad2$days_to_last_follow_up,
                                          clinical_stad2$days_to_death)

# get the gene expression data------
# build a query to get gene expression data from TCGA
query_stad_all <- GDCquery(project = "TCGA-STAD",
                           data.category = "Transcriptome Profiling",
                           experimental.strategy = "RNA-Seq",
                           workflow.type = "STAR - Counts",
                           data.type = "Gene Expression Quantification",
                           sample.type = "Primary Tumor",
                           access = "open")

output_stad <- getResults(query_stad_all)

# we can have a try using 20 samples------
# get 20 primary tissue sample barcodes
tumor <- output_stad[output_stad$sample_type == "Primary Tumor","cases"][1:20]
tumor
# get gene expression data from 20 primary tumors
query_stad <- GDCquery(project = "TCGA-STAD",
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       data.type = "Gene Expression Quantification",
                       sample.type = "Primary Tumor",
                       access = "open",
                       barcode = tumor)

# download data
GDCdownload(query_stad)

# get counts
tcga_stad_data <- GDCprepare(query_stad,summarizedExperiment = TRUE)
stad_matrix <- assay(tcga_stad_data, "unstranded")
stad_matrix[1:10,1:10]

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_stad_data))
coldata <- as.data.frame(colData(tcga_stad_data))

# vst transform counts to be used in survival analysis------
dds <- DESeqDataSetFromMatrix(countData = stad_matrix,
                              colData = coldata,
                              design = ~ 1)# there is no specific design

# removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

# variant stabilizing transformation
vsd <- vst(dds,blind = FALSE)
stad_matrix_vst <- assay(vsd)
stad_matrix_vst[1:10,1:10]

# get data for TP53 gene and add gene metadata information to it------
# we already have geneIDs matched gene symbols in gene_metadata data frame
stad_tp53 <- stad_matrix_vst %>%
  as.data.frame %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id',value = 'counts',-gene_id) %>%
  left_join(., gene_metadata, by = 'gene_id') %>%
  filter(gene_name == "TP53") 
 
# get median value
median_value <- median(stad_tp53$counts)

# denote which cases have higher or lower expression than median count
stad_tp53$strata <- ifelse(stad_tp53$counts >= median_value, "HIGH", "LOW")

# add clinical information to stad_tp53
stad_tp53$case_id <- gsub('-01.*','',stad_tp53$case_id)
stad_tp53 <- merge(stad_tp53,clinical_stad2,by.x='case_id',by.y='submitter_id')

# fitting survival curve
fit <- survfit(Surv(overall_survival,deceased) ~ strata, data = stad_tp53)
fit
ggsurvplot(fit,
           data = stad_tp53,
           pval = T,
           risk.table = T)
fit2 <- survdiff(Surv(overall_survival,deceased) ~ strata, data = stad_tp53)
fit2
