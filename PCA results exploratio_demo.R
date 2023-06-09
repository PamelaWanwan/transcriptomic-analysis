# this script is for PCA and PCA results exploration
library(stringr)
library(ggplot2)

# read in data
load(file = "03_data_out/GSE192560_sub.mrna.exp_symbols.rda")
head(sub.exp)
data.matrix <- as.matrix(sub.exp)
rm(sub.exp)

# set groups
load("01_data_input/GSE192560gset.rda")

pdata <- pData(gset2[[1]])
table(pdata$source_name_ch1)



# run pca
pca <- prcomp(t(data.matrix), scale. = TRUE)
# by default, the prcomp() expects our sample to be rows and the genes(variables)
# to be columns, so we have to transpose the matrix using the t() function

# plot
plot(pca$x[,1], pca$x[,2])#using the first two PCs to draw a 2-d plot

# we use the square of sdev to see how much variation in the original data PC1 accounts for
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
pca.data
pca.data$Sample <- c(paste0(rep("PTC_HT", 5), 1:5), 
                     paste0(rep("PTC",5),1:5))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label= Sample)) + 
  geom_text() + 
  xlab(paste0("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste0("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() + 
  ggtitle("PCA Graph")

# use loading score(rotation) to determine which genes have the largest effect
# on where samples are plotted in the PCA plot

loading_scores <- pca$rotation[,1] # for PC1
# genes that push samples to the left of the graph  will have negative values
# genes that push samples to the right will have large positive values
gene_rank <- sort(loading_scores, decreasing = TRUE)

pos_10_genes <- names(gene_rank[1:10])
pca$rotation[pos_10_genes,1]

neg_10_genes <- names(gene_rank[20375:20384])
pca$rotation[neg_10_genes,1]
