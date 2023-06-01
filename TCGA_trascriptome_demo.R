# script to run differential expression gene analysis using TCGA data

# setwd()
# load packages------
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# build a query to get gene expression data from TCGA------
query_stad <- GDCquery(project = "TCGA-STAD",
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       data.type = "Gene Expression Quantification",
                       access = "open")
output_stad <- getResults(query_stad_all)
GDCdownload(query_stad)

# get counts
tcga_stad_data <- GDCprepare(query_stad,summarizedExperiment = TRUE)
stad_all <- assay(tcga_stad_data, "unstranded")
stad_all[1:10,1:10]

saveRDS(stad_all, "stad_counts.RData")

#readRDS("stad_all_counts.RData")


# build DESeqDataSet------
coldata <- as.data.frame(colData(tcga_stad_data))
table(coldata$sample_type) # Primary Tumor(412) Solid Tissue Normal(36)

dds <- DESeqDataSetFromMatrix(countData = stad_all,
                              colData = coldata,
                              design = ~ sample_type)

# removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]

dds

# set the factor level
dds$sample_type <- relevel(dds$sample_type, ref = "Solid Tissue Normal")

# run DESeq2------
dds <- DESeq(dds)
res <- results(dds)

res

# specified the coefficient or contrast we want to build a results table for,
res <- results(dds, contrast=c("sample_type","Primary Tumor","Solid Tissue Normal"))

# Log fold change shrinkage for visualization and ranking------
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="sample_type_Primary.Tumor_vs_Solid.Tissue.Normal", type="apeglm")
resLFC

# Explore Results ----------------
summary(res)

# res0.01 <- results(dds, alpha = 0.01)
# summary(res0.01)

# save results
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
sum(res$padj < 0.1, na.rm=TRUE)  
sum(res$padj < 0.001, na.rm=TRUE)  

resSig <- subset(resOrdered, padj < 0.001)
head(resSig)

# MA plot
plotMA(resSig)

# order by log2FoldChange
resSigOrderd <- resSig[order(resSig$log2FoldChange,decreasing = TRUE),]
head(resSigOrderd)
tail(resSigOrderd)

saveRDS(resSigOrderd, "resSigOrderd_log2FC.RData")

# gene symbol conversion 
degs <- as.data.frame(resSigOrderd)
degs$ensembl_id <- rownames(degs)
rt <- as.data.frame(degs$ensembl_id)
id <- apply(rt, 1, function(x){
  strsplit(x, "[.]")[[1]][1]
})
id <- as.data.frame(id)
rt$ensembl_id <- id$id
library(org.Hs.eg.db)
rt$gene_symbol<- mapIds(org.Hs.eg.db,
                        keys = rt$ensembl_id,
                        keytype = 'ENSEMBL',
                        column = 'SYMBOL')
colnames(rt) <- c("ensembl_id","ensembl","gene_symbol")
degs <- merge(degs,rt,by="ensembl_id",all.x=T)
# delete NA
table(is.na(degs$gene_symbol))
newdegs<-degs[complete.cases(degs$gene_symbol),]
# remove duplicates
table(newdegs$gene_symbol) > 1
which(duplicated(newdegs$gene_symbol))
degs.2 <- newdegs[!duplicated(newdegs$gene_symbol),]
table(table(degs.2$gene_symbol) > 1)

saveRDS(degs.2,"degs_nodups_noNA.RData")

rownames(degs.2) <- degs.2$gene_symbol
degs.2 <- degs.2[,-c(1,8,9)]

# order by log2FC
degs.2orderd <- degs.2[order(degs.2$log2FoldChange,decreasing = TRUE),]
table(degs.2orderd$log2FoldChange > 2)#1078 uploaded genes
table(abs(degs.2orderd$log2FoldChange) >2)# 1067 downloaded genes


# DEG visulization------

#plotCounts
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
updegs <- head(degs.2orderd)
downdegs <- tail(degs.2orderd)
# up*3 down*3
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000213401.10", intgroup="sample_type", main = "MAGEA12") # MAGEA12
plotCounts(dds, gene="ENSG00000123407.4", intgroup="sample_type", main = "HOXC12") # HOXC12 
plotCounts(dds, gene="ENSG00000126890.13", intgroup="sample_type", main = "CTAG2") # CTAG2
plotCounts(dds, gene="ENSG00000229183.8", intgroup="sample_type", main = "PGA4") # PGA4
plotCounts(dds, gene="ENSG00000182156.10", intgroup="sample_type", main = "ENPP7") # ENPP7
plotCounts(dds, gene="ENSG00000172689.2", intgroup="sample_type", main = "MS4A10") # MS4A10

help("plotCounts")

#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


# PCA
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sample_type")

#heatmap
library(ph)
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(dds[1:100], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = meta, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

