# GSVA

# load packages
if(!"XML"%in%installed.packages())
  install.packages("XML", type = "binary")
library(GSEABase)
library(GSVA)
library(limma)

# read in data
load("01_data_input/GSE192560_exp_symbols.rda")
exp[1:6,1:6]
table(is.na(exp))
exp <- as.matrix(exp)

# using MsigDb database
msigdb_GMTs <- "msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs"
msigdb <- "msigdb.v2023.1.Hs.symbols.gmt"

geneset <- getGmt(file.path(msigdb_GMTs,msigdb))
head(geneset)

# perform gsva
ex.max <- gsva(exp,
               geneset,
               mx.diff=F,
               verbose=F,
               parallel.sz=1)
head(ex.max)
write.table(ex.max, file="03_data_out/GSVA_MsigDb_result.txt",sep="\t",
            quote=F,row.names = T)
save(ex.max, file = "03_data_out/GSVA_MsigDb_result.Rda")

# visualization
remove(list = ls())
load(file = "03_data_out/GSVA_MsigDb_result.Rda")
dim(ex.max)

# set groups
group_list <- c(rep("PTC_HT", 5), rep("PTC",5))
table(group_list)

design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(ex.max)
design

contrast.matrix<-makeContrasts("PTC_HT-PTC",#比较矩阵，PTC_HT和PTC相比
                               levels = design)#设计矩阵

contrast.matrix

# extract degs
fit <- lmFit(ex.max,design)#线性拟合
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#贝叶斯

res <- decideTests(fit2, p.value=0.05)
summary(res)

tmp <- topTable(fit2)
DEG <- na.omit(tmp)
head(DEG)

write.table(DEG, file="03_data_out/GSVA_DEG_limma_result.txt",sep="\t",
            quote=F,row.names = T)
save(DEG, group_list, file = "03_data_out/GSVA_DEG_limma.Rda")

# heatmap
#load(file = "03_data_out/GSVA_DEG_limma.Rda")
#load(file = "03_data_out/GSVA_MsigDb_result.Rda")

DEG_sig <- DEG[DEG$P.Value<0.01 & abs(DEG$logFC) > 0.3,]
dat <- ex.max[match(rownames(DEG_sig),rownames(ex.max)),]
annotation_col <- data.frame(group_list)
rownames(annotation_col) <- colnames(dat)
pheatmap::pheatmap(dat,
                   width = 20,
                   height = 11,
                   annotation_col = annotation_col,
                   show_colnames = F,
                   filename = '02_plot_out/GSVA_heatmap.png')
