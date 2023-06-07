# this script is to perform gene ontology analysis
# example data: GSE192560

# load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

options(stringsAsFactors = F)

# read-in deg data
load(file = "03_data_out/GSE192560_degs.rda")

# extract up regulated genes
up_genes <- rownames(deg[deg$change=="up",])
length(up_genes)

# set filter condition
pvalueFilter=0.05        #filter for pvalue
qvalueFilter=0.05        #filter for qvalue
showNum=10

# run GO analysis
GO_results <- enrichGO(gene = up_genes, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", ont = "all",
                       pAdjustMethod = "BH", # Benjamini and Hochberg for p-value adjusting
                       pvalueCutoff = pvalueFilter,
                       qvalueCutoff = qvalueFilter )

# save up_genes GO results
GO <- as.data.frame(GO_results)
GO <- GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="03_data_out/GSE192560_GO_ALL.txt",sep="\t",quote=F,row.names = F)

# plot GO results
# calculate FoldEnrichment 
GR_BG <- function(ratio){
  sapply(ratio,function(x) as.numeric(gsub("/.*$","",x))/as.numeric(gsub("^.*/","",x)))
}

fold_enrichment <- GR_BG(GO$GeneRatio)/GR_BG(GO$BgRatio)

GO$FoldEnrichment <- fold_enrichment

GO$pvalue <- as.numeric(GO$pvalue)

# extract "biological process", "cellular component", "molecular function" separately
go.data.bp <- GO[GO$ONTOLOGY == "BP",]
go.data.cc <- GO[GO$ONTOLOGY == "CC",]
go.data.mf <- GO[GO$ONTOLOGY == "MF",]

# Sort each table in ascending order by p-value value
go.data.bp <- go.data.bp[order(go.data.bp$pvalue),]
go.data.cc <- go.data.cc[order(go.data.cc$pvalue),]
go.data.mf <- go.data.mf[order(go.data.mf$pvalue),]

#select the first 10 rows
go.data.bp.top <- head(go.data.bp, n = showNum)
go.data.cc.top <- head(go.data.cc, n = showNum)
go.data.mf.top <- head(go.data.mf, n = showNum)

# sort each table in descending order by foldenrichment vlues
go.data.bp.top <- go.data.bp.top[order(go.data.bp.top$FoldEnrichment, decreasing = TRUE), ]
go.data.cc.top <- go.data.cc.top[order(go.data.cc.top$FoldEnrichment, decreasing = TRUE), ]
go.data.mf.top <- go.data.mf.top[order(go.data.mf.top$FoldEnrichment, decreasing = TRUE), ]

# merge them into one dataframe
go.data.top <- rbind(go.data.bp.top, go.data.cc.top, go.data.mf.top)

write.csv(go.data.top, "03_data_out/GOenrichment_top10.csv", row.names = FALSE)
write.table(go.data.top,file="03_data_out/GOenrichment_top10.txt",sep="\t",quote=F,row.names = F)

# for visualization
go.data.top$Description <- factor(go.data.top$Description,levels = rev(go.data.top$Description))
# bubble
go_bubble <- ggplot(go.data.top, aes(x = FoldEnrichment, y = Description, color = -log10(pvalue))) +
  geom_point(aes(size = Count)) + # 绘制散点图
  theme_bw() + # 白色主题
  scale_y_discrete(labels = function(y) strwrap(y, width = 60)) + # 设置Term名称过长时换行
  labs(size = "Counts", x = "Fold Enrichment", y = "GO terms", title = "GO Enrichment") +
  scale_color_gradient(low="blue",high ="red") + # 颜色渐变设置
  theme(axis.text = element_text(size = 10, color = "black"), # 轴标签大小和颜色
        axis.title = element_text(size = 16), # 轴标题大小
        title = element_text(size = 13)) + # 图标题大小
  facet_grid(ONTOLOGY~., scales = "free", space = "free") + # 按照ONTOLOGY类型分面展示
  guides(
    #reverse color order (higher value on top)
    color = guide_colorbar(reverse = TRUE))
#reverse size order (higher diameter on top) 
#size = guide_legend(reverse = TRUE))
go_bubble

ggsave(go_bubble,filename = "02_plot_out/go_bubble.pdf",width = 9,height = 9)
ggsave(go_bubble,filename = "02_plot_out/go_bubble.png",width = 9,height = 9)
ggsave(go_bubble, filename = "02_plot_out/go_bubble.tiff", width = 9, height = 9, dpi = 600)
  
# barplot
go.data.bp.top_bar <- go.data.top[go.data.top$ONTOLOGY == "BP",]
go.data.cc.top_bar <- go.data.top[go.data.top$ONTOLOGY == "CC",]
go.data.mf.top_bar <- go.data.top[go.data.top$ONTOLOGY == "MF",]

# sort each table in descending order by counts
go.data.bp.top_bar <- go.data.bp.top_bar[order(go.data.bp.top_bar$`Count`, decreasing = TRUE), ]
go.data.cc.top_bar <- go.data.cc.top_bar[order(go.data.cc.top_bar$`Count`, decreasing = TRUE), ]
go.data.mf.top_bar <- go.data.mf.top_bar[order(go.data.mf.top_bar$`Count`, decreasing = TRUE), ]

# merge them together
go.data.top_bar <- rbind(go.data.bp.top_bar, go.data.cc.top_bar, go.data.mf.top_bar)


# plot
go.data.top_bar$Description <- factor(go.data.top_bar$Description,levels = rev(go.data.top_bar$Description))
go_bar <- ggplot(data = go.data.top_bar, aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_col(width = 0.9) + # 绘制柱状图，0.9为柱子宽度
  coord_flip() + # 翻转y轴
  theme_bw() + # 黑白主题
  scale_x_discrete(labels = function(x) strwrap(x, width = 60)) + # 设置x轴标签换行，避免标签重叠
  labs(x = "GO terms", y = "Counts", title = "GO Enrichment") + # 设置x轴、y轴、标题标签
  theme(axis.title = element_text(size = 13), # 设置坐标轴标签字体大小
        axis.text = element_text(size = 11), # 设置坐标轴刻度字体大小
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 设置标题字体大小、位置和加粗
        legend.title = element_text(size = 13), # 设置图例标题字体大小
        legend.text = element_text(size = 11), # 设置图例标签字体大小
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + # 设置图形边距
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a"), # 设置图例填充颜色
                    breaks = c("BP", "CC", "MF"), # 设置图例标签
                    labels = c("Biological Process", "Cellular Component", "Molecular Function")) # 设置图例标签文字
go_bar
# 保存柱状图
ggsave(go_bar,filename = "02_plot_out/go_barplot.pdf",width = 9,height = 9)
ggsave(go_bar,filename = "02_plot_out/go_barplot.png",width = 9,height = 9)
ggsave(go_bar, filename = "02_plot_out/go_bubble.tiff", width = 9, height = 9, dpi = 600)
