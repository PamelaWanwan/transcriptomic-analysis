# GSEA

# load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(enrichplot)

# create input data------
load("03_data_out/GSE192560_degs.rda")
head(deg)

deg$symbols <- rownames(deg)
head(deg)

deg.input <- deg[,c(8,1)]
head(deg.input)

# convert gene symbol into entrzid
deg.input$ENTREZID<- deg.input$symbols%>%
  mget(org.Hs.egSYMBOL2EG, ifnotfound = NA) %>%
  as.character()

deg.input <- deg.input[,c(3,2)]
table(is.na(deg.input))
head(deg.input)

# Sort the data by logFC
geneList <- deg.input[,2]
names(geneList) <- as.character(deg.input[,1])
geneList <- sort(geneList,decreasing = T)
head(geneList)

# gseGO------
gsego <- gseGO(geneList=geneList,
               OrgDb = org.Hs.eg.db,
               ont = "ALL",
               nPerm = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(gsego)
write.table(gsego,file="03_data_out/GSE192560_degs_GSEA_GO_result.txt",sep="\t",
            quote=F,row.names = F)

# gseKEGG------
gsekk <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
head(gsekk)
write.table(gsekk,file="03_data_out/GSE192560_degs_GSEA_KEGG_result.txt",sep="\t",
            quote=F,row.names = F)

# gsea of MSigDb------
msigdb_GMTs <- "msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs"
msigdb <- "msigdb.v2023.1.Hs.entrez.gmt"

all_msigdb <- read.gmt(file.path(msigdb_GMTs,msigdb))

gsegmt <- GSEA(geneList, 
               TERM2GENE=all_msigdb, 
               verbose=FALSE)
head(gsegmt)
write.table(gsegmt,file="03_data_out/GSEA_MSigDb_result.txt",sep="\t",
            quote=F,row.names = F)

save(geneList, gsego, gsekk, gsegmt, file = "03_data_out/gsea_result.Rda")

# visualization------
remove(list = ls())
load(file = "03_data_out/gsea_result.Rda")

# dotplot
tiff(file="02_plot_out/GSEA_GO_dotplot.tiff",
     width = 120,height = 90,units ="cm",
     compression="lzw",bg="white",res=300)
dotplot(gsego, showCategory=30) + ggtitle("dotplot for GSEA")
dev.off()
dotplot(gsekk, showCategory=30)
dotplot(gsegmt, showCategory=30)

tiff(file="02_plot_out/GSEA_KEGG_cnetplot.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
cnetplot(gsekk, categorySize="pvalue", foldChange=geneList)
dev.off()
cnetplot(gsekk, foldChange=geneList)
cnetplot(gsekk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

tiff(file="02_plot_out/GSEA_GO_heatplot.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
heatplot(gsego, foldChange=geneList)
dev.off()

########## GSEA图：gseaplot、gseaplot2、gsearank ##########
## GSEA图：gseaplot
gseaplot(gsekk, geneSetID = 1, by = "runningScore", title = gsekk$Description[1])#runningScore光展示下面一个图
gseaplot(gsekk, geneSetID = 1, by = "preranked", title = gsekk$Description[1])#preranked光展示上面一个图
tiff(file="02_plot_out/GSEA_KEGG_gseaplot.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
gseaplot(gsekk, geneSetID = 1, title = gsekk$Description[1])#geneSetID的1和gsekk$Description[1]的1是对应的
dev.off()

## GSEA图：gseaplot2
gseaplot2(gsekk, geneSetID = 1, title = gsekk$Description[1])
## 在一个gsea图上绘制多个基因集
gseaplot2(gsekk, geneSetID = 1:3)#绘制三个基因集
gseaplot2(gsekk, geneSetID = 1:3, subplots = 1)#绘制第一部分
gseaplot2(gsekk, geneSetID = 1:3, subplots = 1:2)#绘制两部分
tiff(file="02_plot_out/GSEA_KEGG_gseaplot2.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
gseaplot2(gsekk, geneSetID = 1:3, subplots = 1:3,#geneSetID，gsekk中富集到的基因集编号，1：3表示前三个
          pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), #指定颜色
          ES_geom = "dot")#line 连线图；dot虚线图
dev.off()

########## 山脊图：ridgeplot ##########
tiff(file="02_plot_out/GSEA_KEGG_ridgeplot.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
ridgeplot(gsekk)
dev.off()
tiff(file="02_plot_out/GSEA_GO_ridgeplot.tiff",
     width = 35,height = 22,units ="cm",
     compression="lzw",bg="white",res=300)
ridgeplot(gsego)
dev.off()
