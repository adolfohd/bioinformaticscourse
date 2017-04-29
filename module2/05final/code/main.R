setwd("/media/fito/Windows/Users/fitoh/Documents/code/doc/bioinformaticscourse/module2/05final/code")
rm(list = ls())
source("/home/fito/code/fito_lib/r/pkgtest.R")

my.install.packages <- function(pkg) if (!require(pkg)) install.packages(pkg);

my.hclust <- function(x, method = "complete", dmeth = "euclidean"){
  hclust(dist(x, method = dmeth), method = method)
}

#### dependencies

pkgTest('dendextend')
pkgTest('colorspace')
pkgTest('knitr')
pkgTest('ggplot2')
library(dendextend)
library(colorspace)
library(knitr)
#### main ##############
k <- 4

pathways.scores <- as.data.frame(read.table("../data/pathways.scores"))
summary(pathways.scores)
colnames(pathways.scores)

# --- Hierarchical clustering
hierarch.clusters <- my.hclust(pathways.scores)

hierarch.cluster.dendrogram <- as.dendrogram(hierarch.clusters)
hierarch.cluster.dendrogram <- color_branches(hierarch.cluster.dendrogram, k = k)
# plot(cut(cluster.dendrogram, h=4)$upper, main = "cut dendrogram")

hierarch.clusters.cut <- cutree(hierarch.cluster.dendrogram, k)
hierarch.clusters.table <- table(hierarch.clusters.cut)
# kable(hierarch.clusters.table)

#------ k-means

kmeans.clusters <- kmeans(pathways.scores, centers = 4, nstart = 30)
plot(kmeans.clusters)

kmeans.clusters.table <- table(kmeans.clusters$cluster)

kmeans.clusters$cluster <- as.factor(kmeans.clusters$cluster)


# ----- Jackard Coefficient




# ---- PLOTS -----

plot.colors <- rainbow_hcl(k)

jpeg('../images/hierarch.dendrogram.jpg')
plot(hierarch.cluster.dendrogram, main = "original dendrogram")
dev.off()

jpeg('../images/hierarch.histogram.jpg')
bp <- barplot(hierarch.clusters.table, col = plot.colors, xlab = "clusters", ylab = "number of elements", main = "histogram of resulting clustering")
text(bp, 0, hierarch.clusters.table, cex = 1, pos = 3)
dev.off()

jpeg('../images/hierarchical.clustering.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = as.factor(hierarch.clusters.cut))) + geom_point() + ggtitle("Hierarchical clustering results")
dev.off()

jpeg('../images/kmeans.histogram.jpg')
bp <- barplot(kmeans.clusters.table, col = plot.colors, xlab = "clusters", ylab = "number of elements", main = "histogram - kmeans")
text(bp, 0, kmeans.clusters.table, cex = 1, pos = 3)
dev.off()

jpeg('../images/kmeans.clustering.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = kmeans.clusters$cluster)) + geom_point() + ggtitle("K-means clustering results")
dev.off()
