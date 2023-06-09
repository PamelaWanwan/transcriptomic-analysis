# GEO chip-seq data analysis from download to differential visulization
# exmaple: GSE39582

# load packages------
library(GEOquery) # for downloading data
library(stringr) # for setting group levels
library(limma) # for normalization and differencial analysis
library(hgu133plus2.db) # for converting probes to symbols
library(tidyverse) # for manipulating data
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(ggpubr) # for DEGs visulization volcano plot
library(ggthemes) # for DEGs visulization volcano plot
library(ggplot2) # for DEGs visulization volcano plot
library(pheatmap) # for DEGs visulization  heatmap

# read in matrix data------
# options('download.file.method.GEOquery'='auto')
# options('GEOquery.inmemory.gpl'=FALSE)
# options('download.file.method.GEOquery' = 'libcurl')
# gset <- getGEO("GSE39582", destdir = "02_data input/", AnnotGPL = F, getGPL = F)

gset <- getGEO("GSE39582", destdir = "02_data input/") 


# get expression and clinical information------
# get phenoData information
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)

# set treated samples vs control samples
group_list <- ifelse(str_detect(pdata$source_name_ch1, "Adenocarcinoma"), "tumor", 
                     "normal")
# set reference level
group_list <-  factor(group_list, levels = c("normal","tumor"))
table(group_list)

# get expression matrix from gset and perform normalization------
exp <- exprs(gset[[1]])
# explore the entire expressing phenomenon
boxplot(exp, outline = FALSE, notch = T, col = group_list, las = 2)
# normalization
exp <- normalizeBetweenArrays(exp)
boxplot(exp, outline = FALSE, notch = T, col = group_list, las = 2)
range(exp)#minimum and maximum, Logarithmic conversion has been performed

# switch probes to gene symbols------
# check the platform information
index <- gset[[1]]@annotation
# GPL570 is the most common platform and it matches with hgu133plus2.db

# how much information "package:hgu133plus2.db" contains
ls("package:hgu133plus2.db")
#"hgu133plus2ENSEMBL"；"hgu133plus2ENTREZID"；"hgu133plus2ENTREZID"

# extract ids from this package
ids <- toTable(hgu133plus2SYMBOL)
head(ids) #it contains probe ids and symbols
length(ids$symbol)#43101
length(unique(ids$symbol))#20824
table(sort(table(ids$symbol)))

# match
exp <- as.data.frame(exp)

exp <- exp %>%
  mutate(probe_id = rownames(exp)) %>%
  inner_join(., ids, by = "probe_id") %>%
  select(probe_id, symbol, everything())

# remove duplicates
exp <- exp[!duplicated(exp$symbol),]

# rename rows
rownames(exp) <- exp$symbol

# remove probe ids and symbols
exp <- exp[,-(1:2)]

exp[1:3,1:3]

# PCA------
# to see if there is huge difference between normal and tumor samples
# transpose dataframe
dat <- as.data.frame(t(exp))

# do PCA
dat.pca <- PCA(dat, graph = FALSE) 

# PCA visulization
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point",
                         col.ind = group_list,
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE,
                         legend.title = "Groups")
pca_plot
# there are differences between normal and tumor to some extent

# find DEGs------
#input file 1: expression matrixz(rownames:symbol,colnames:sample)
#input file 2: group list

design <- model.matrix(~group_list)
fit <- lmFit(exp, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 2, number = Inf)

colnames(deg)

logFC <- 1.5
P.Value <- 0.05
type1 <- (deg$P.Value < P.Value)&(deg$logFC < -logFC)# down regulated
type2 <- (deg$P.Value < P.Value)&(deg$logFC > logFC)# up regulated

deg$change <- ifelse(type1, "down", ifelse(type2, "up", "stable")) # divide into up down stable three groups

table(deg$change)

# DEGs visulization volcano plot------
# volcano plot based on ggpubr
deg$logP <- -log10(deg$P.Value)
ggscatter(deg, x = "logFC", y = "logP",
          color = "change",
          palette = c("blue","black","red"),
          size = 1) + 
  theme_base() + 
  geom_hline(yintercept = -log10(P.Value), linetype = "dashed") + 
  geom_vline(xintercept = c(-logFC, logFC), linetype = "dashed")

# add details to top 5 regulated gene
deg$Label = "" #add a new column of label
deg <- deg[order(deg$logFC),]

head(deg)
tail(deg)

up.gene <- tail(rownames(deg)[which(deg$change == "up")], 5)
up.gene

down.gene <- head(rownames(deg)[which(deg$change == "down")], 5)
down.gene

# combine up and down gene and save them to label
deg.top5.genes <- c(as.character(up.gene), as.character(down.gene))
deg$Label[match(deg.top5.genes, rownames(deg))] <- deg.top5.genes

head(deg)

help("match")

# plot
ggscatter(deg, x = "logFC", y = "logP",
          color = "change",
          palette = c("blue","black","red"),
          size = 1,
          label = deg$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(P-value)") + 
  theme_base() + 
  geom_hline(yintercept = -log10(P.Value), linetype = "dashed") + 
  geom_vline(xintercept = c(-logFC, logFC), linetype = "dashed")

# volcano plot based on ggplot2
p <- ggplot(data = deg,
            aes(x = logFC,
                y = -log10(P.Value))) + 
  geom_point(alpha = 0.4, size = 3.5,
             aes(color = change)) + 
  ylab("-log10(Pvalue)") + 
  scale_color_manual(values = c("blue","black","red")) + 
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) + 
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) + 
  theme_bw()
p
# also add top gene labels
x1 <- deg %>%
  filter(change == "up") %>%
  tail(5)

x2 <- deg %>%
  filter(change == "down") %>%
  head(5)

label <- rbind(x1, x2)

volcano_plot <- p + 
  geom_point(size = 3, shape =1, data = label) + 
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color = "black")

volcano_plot

# DEGs visulization heatmap------
cg <- rownames(deg)[deg$change !="stable"]
diff <- exp[cg,]

 # extract differential expression data
annotation_col <- data.frame(group = group_list)
rownames(annotation_col) <- colnames(diff)

pheatmap(diff,
         annotation_col = annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_col = 3,
         fontsize_row = 3)

# The dowstream procedure including VENN GO KEGG GSEA PPI WGCNA etc.
