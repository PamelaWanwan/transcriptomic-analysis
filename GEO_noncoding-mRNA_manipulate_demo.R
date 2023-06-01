# this script is to manipulate Non-coding RNA profiling by array plus Expression profiling by array data
# extract expression data from it and perform deg analysis

# GSE192560 

# load packages------
library(GEOquery) # for downloading data
library(AnnoProbe) # for annotation conversion 
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
#gset2 <- getGEO("GSE192560")
#save(gset2, file = "GSE192560gset.rda")

# get expression matrix from gset and perform normalization------
load("GSE192560gset.rda")
exp <- exprs(gset2[[1]])
# explore the entire expressing phenomenon
boxplot(exp, outline = FALSE, notch = T, las = 2) 

# extract mRNA expression data
fdata <- gset2$GSE192560_series_matrix.txt.gz@featureData
des <- fdata@data
save(des, file = "GPL16956.annotation.infromation.rda")

table(des$TRANSCRIPT_TYPE)

# noncoding protein_coding 
# 25247          25008 

des <- des[,1:2];des
prot <- des[des$TRANSCRIPT_TYPE=="protein_coding",]

sub.exp <- exp[prot$ID,]
save(sub.exp, file = "GSE192560.mRNA.exp.rda")


gpl <- 'GPL16956'

probe2gene <- idmap(gpl,type = 'pipe')

head(probe2gene)
length(probe2gene$symbol)#73926
length(unique(probe2gene$symbol))#33591
table(sort(table(probe2gene$symbol)))

# match
sub.exp <- as.data.frame(sub.exp)
sub.exp <- sub.exp %>%
  mutate(probe_id = rownames(sub.exp)) %>%
  inner_join(., probe2gene, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

# remove duplicates
sub.exp <- sub.exp[!duplicated(sub.exp$symbol),]

# rename rows
rownames(sub.exp) <- sub.exp$symbol

# remove probe ids and symbols
sub.exp <- sub.exp[,-(1:2)]
sub.exp[1:3,1:3]
save(sub.exp, file = "GSE192560_sub.mrna.exp_symbols.rda")
head(sub.exp)

pdata <- pData(gset2[[1]])
table(pdata$source_name_ch1)


# set PTC_HT samples vs PTC samples
group_list <- ifelse(str_detect(pdata$source_name_ch1, "with Hashimoto's thyroiditis"), 
                     "PTC_HT", "PTC")
table(group_list)

# set reference level
group_list <-  factor(group_list, levels = c("PTC", "PTC_HT"))


# find DEGs------
#input file 1: expression matrixz(rownames:symbol,colnames:sample)
#input file 2: group list
design <- model.matrix(~group_list)
fit <- lmFit(sub.exp, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef = 2, number = Inf)

head(deg)

logFC <- 1
P.Value <- 0.05
type1 <- (deg$P.Value < P.Value)&(deg$logFC < -logFC)#down regulated
type2 <- (deg$P.Value < P.Value)&(deg$logFC > logFC)#up regulated

deg$change <- ifelse(type1, "down", ifelse(type2, "up", "stable")) # divide into up down stable three groups

table(deg$change)

save(deg, file = "GSE192560_degs.rda")

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
deg$Label = "" # add a new row of label
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

which(rownames(sub.exp)=="CXCL13")
sub.exp[3793,]
mean(as.numeric(sub.exp[3793,(1:5)]))#ptcht 10.8
mean(as.numeric(sub.exp[3793,(6:10)]))#ptc 5.5

# DEGs visulization heatmap------
cg <- rownames(deg)[deg$change !="stable"]
cg
diff <- sub.exp[cg,]
diff <- diff

# extract differential expression data
annotation_col <- data.frame(group = group_list)
group_list
rownames(annotation_col) <- c(colnames(diff))

pheatmap(diff,
         annotation_col = annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames = T,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_col = 3,
         fontsize_row = 3)
