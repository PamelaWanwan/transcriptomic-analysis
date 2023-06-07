# this script is to perform KEGG analysis

# load packages
# example data: GSE192560

# load packages
library(pacman)
p_load(colorspace,stringi,data.table, ggplot2,tidyverse,ggvenn,DOSE,
       clusterProfiler,enrichplot,org.Hs.eg.db,KEGGgraph,Rgraphviz,pathview)
options(stringsAsFactors = F)

# read-in deg data
load(file = "03_data_out/GSE192560_degs.rda")

# extract up regulated genes
up_genes <- rownames(deg[deg$change=="up",])
length(up_genes)


# convert gene symbols to entrezIDs
entrezIDs <- up_genes %>%
  mget(org.Hs.egSYMBOL2EG, ifnotfound = NA) %>%
  as.character()
  
up_genes <- as.data.frame(up_genes)
up_genes$entrezIDs <- entrezIDs

my.data <- up_genes[up_genes$entrezIDs!="NA",]

# set filter condition
pvalueSET  <- 0.05 
qvalueSET <- 0.05 
showNum <- 30 

# perform enrichment analysis
enrich_result <- enrichKEGG(gene = my.data$entrezIDs, 
                            organism = "hsa", 
                            pvalueCutoff = pvalueSET,
                            qvalueCutoff = qvalueSET) 

KEGG <- as.data.frame(enrich_result)

# convert entrezID to gene symbols
KEGG$geneID <- as.character(sapply(KEGG$geneID,function(x){
  ids = strsplit(x,"/")[[1]]
  idx = match(ids, as.character(my.data$entrezIDs))
  symbols = my.data$up_genes[idx]
  paste(symbols, collapse = "/")
}))

write.table(KEGG,file="03_data_out/KEGG.txt",sep="\t",quote=F,row.names = F)

# Visualization of results
KEGG <- KEGG[order(KEGG$Count, decreasing = TRUE), ]
KEGG.top <- head(KEGG, n = showNum)
KEGG.top <- KEGG.top[order(KEGG.top$pvalue, decreasing = TRUE), ]
KEGG.top$Description <- factor(KEGG.top$Description,
                               levels = c(KEGG.top$Description %>% 
                                            as.data.frame() %>% 
                                            pull()))

kegg_bubble <- ggplot(KEGG.top, aes(x = -log10(pvalue), y = Description, color = -log10(pvalue))) +
  geom_point(aes(size = Count)) + # 绘制散点图，x轴为FoldEnrichment，y轴为Description，颜色表示-log10(pvalue)，点的大小表示Count
  theme_bw() + # 白色主题
  scale_y_discrete(labels = function(y) strwrap(y, width = 70)) + # 设置PATHWAY名称过长时换行
  labs(size = "Counts", x = "-log10(pvalue)", y = "KEGG PATHWAY", title = NULL) + # 添加标签和标题
  scale_color_gradient(low="blue",high ="red") + # 颜色渐变设置
  theme(axis.text = element_text(size = 10, color = "black"), # 轴标签大小和颜色
        axis.title = element_text(size = 16), # 轴标题大小
        title = element_text(size = 13)) + # 图标题大小
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE)) # 颜色渐变设置，高值颜色在下方，低值颜色在上方
#reverse size order (higher diameter on top) 
#size = guide_legend(reverse = TRUE))

kegg_bubble

ggsave(kegg_bubble,filename = "02_plot_out/kegg_bubble.pdf",width = 9,height = 9)
ggsave(kegg_bubble,filename = "02_plot_out/kegg_bubble.png",width = 9,height = 9)
ggsave(kegg_bubble, filename = "02_plot_out/kegg_bubble.tiff", width = 9, height = 9, dpi = 600)


kegg_bar <- ggplot(data = KEGG.top, aes(x = Description, y = Count, fill = -log10(pvalue))) +
  geom_col(width = 0.9) +  # 绘制柱状图
  coord_flip() +  # 坐标轴翻转
  theme_bw() +  # 白色主题
  scale_x_discrete(labels = function(x) strwrap(x, width = 70)) +  # 设置PATHWAY名称过长时换行
  labs(x = "KEGG PATHWAY", y = "Counts", title = NULL) +  # 设置横轴纵轴和标题
  theme(axis.title = element_text(size = 13),  # 设置轴标题大小
        axis.text = element_text(size = 11),  # 设置轴标签大小
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),  # 设置图标题大小和对齐方式
        legend.title = element_text(size = 13),  # 设置图例标题大小
        legend.text = element_text(size = 11, angle = 0),  # 设置图例标签大小和方向
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),  # 设置图形四周的留白
        legend.position = "right",  # 设置图例位置为右侧
        legend.box = "vertical",  # 设置图例方向为竖直
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)) +  # 设置图例四周的留白
  scale_fill_gradient(low = "blue", high = "red", name = "-log10(pvalue)") +  # 颜色渐变设置
  guides(
    #reverse fill order (higher value on top)
    fill = guide_colorbar(reverse = TRUE))  # 设置图例颜色渐变反转

kegg_bar
