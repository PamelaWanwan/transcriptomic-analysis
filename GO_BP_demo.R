# this script is to perform gene ontology analysis
# example data: GSE192560

# load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# read-in deg data
load(file = "03_data_out/GSE192560_degs.rda")

# extract up regulated genes
up_genes <- rownames(deg[deg$change=="up",])
length(up_genes)

# run GO analysis
GO_results_BP <- enrichGO(gene = up_genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

# set filter condition
pvalueFilter=0.05        #filter for pvalue
qvalueFilter=0.05        #filter for qvalue

# set colors
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

# save up_genes GO BP results
GO_BP <- as.data.frame(GO_results_BP)
GO_BP <- GO_BP[(GO_BP$pvalue<pvalueFilter & GO_BP$qvalue<qvalueFilter),]
write.table(GO_BP,file="03_data_out/GSE192560.txt",sep="\t",quote=F,row.names = F)

# plot GO BP results
# barplot
barplot(GO_results_BP, showCategory = 15, color = colorSel)
help("barplot")
