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

pathways.values <- as.data.frame(read.table("../data/pathways.values"))
pathways.scores <- as.data.frame(read.table("../data/pathways.scores"))

summary(pathways.scores)
colnames(pathways.scores)

# Hierarchical clustering

hierarch.clusters.reduced <- my.hclust(pathways.scores)
hierarch.clusters.structural <- my.hclust(pathways.values)

hierarch.cluster.dendrogram.reduced <- as.dendrogram(hierarch.clusters.reduced)
hierarch.cluster.dendrogram.reduced <- color_branches(hierarch.cluster.dendrogram.reduced, k = k)
hierarch.cluster.dendrogram.structural <- as.dendrogram(hierarch.clusters.structural)
hierarch.cluster.dendrogram.structural <- color_branches(hierarch.cluster.dendrogram.structural, k = k)

hierarch.clusters.cut.reduced <- cutree(hierarch.cluster.dendrogram.reduced, k)
hierarch.clusters.cut.structural <- cutree(hierarch.cluster.dendrogram.structural, k)

hierarch.clusters.table.reduced <- table(hierarch.clusters.cut.reduced)
hierarch.clusters.table.structural <- table(hierarch.clusters.cut.structural)

## Jackard coefficient

my.jaccard.coefficient <- function(labels1, cluster1, labels2, cluster2){
  labels1.comparisson <- labels1 == cluster1
  labels2.comparisson <- labels2 == cluster2
  return (sum(labels1.comparisson & labels2.comparisson ) / sum(labels1.comparisson | labels2.comparisson ))
}

my.jaccard.matrix <- function(labels1, labels2){
  labels1 <- as.array(labels1)
  labels2 <- as.array(labels2)
  jaccard.hierarchical <- matrix(nrow = k, ncol = k)
  max.jaccards <- array()
  for (i in 1:k){
    max.jaccards[i] <- 0
    for(j in 1:k){
      jaccard.hierarchical[i,j] <- my.jaccard.coefficient(labels1, i, labels2, j)
      print(max.jaccards[i])
      if (jaccard.hierarchical[i,j] > max.jaccards[i])
        max.jaccards[i] <- j
    }
  }
  return(list(jaccard.hierarchical, max.jaccards))
}

jaccard.matrix <- my.jaccard.matrix(hierarch.clusters.cut.structural, hierarch.clusters.cut.reduced)
print(jaccard.matrix )

# ---- PLOTS -----

plot.colors <- rainbow_hcl(k)

jpeg('../images/hierarch.dendrogram.structural.jpg')
plot(hierarch.cluster.dendrogram.structural, main = "Hierarchical clustering, structural data")
dev.off()

jpeg('../images/hierarch.dendrogram.reduced.jpg')
plot(hierarch.cluster.dendrogram.reduced, main = "Hierarchical clustering, reduced data")
dev.off()

jpeg('../images/hierarch.histogram.structural.jpg')
bp <- barplot(hierarch.clusters.table.structural, col = plot.colors,
              xlab = "clusters", ylab = "number of elements", 
              main = "structural data: hierarchical clustering histogram")
text(bp, 0, hierarch.clusters.table.structural, cex = 1, pos = 3)
dev.off()

jpeg('../images/hierarch.histogram.reduced.jpg')
bp <- barplot(hierarch.clusters.table.reduced, col = plot.colors, 
              xlab = "clusters", ylab = "number of elements", 
              main = "reduced data: hierarchical clustering histogram")
text(bp, 0, hierarch.clusters.table.reduced, cex = 1, pos = 3)
dev.off()

jpeg('../images/hierarchical.clustering.structural.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = as.factor(hierarch.clusters.cut.structural))
       ) + geom_point() + ggtitle("Structural data: Hierarchical clustering")
dev.off()

jpeg('../images/hierarchical.clustering.reduced.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = as.factor(hierarch.clusters.cut.reduced))
       ) + geom_point() + ggtitle("Reduced data: Hierarchical clustering")
dev.off()

# K-means

kmeans.clusters.reduced <- kmeans(pathways.scores, centers = 4, nstart = 30)
plot(kmeans.clusters.reduced)

kmeans.clusters.table.reduced <- table(kmeans.clusters.reduced$cluster)

kmeans.clusters.reduced$cluster <- as.factor(kmeans.clusters.reduced$cluster)



## Plots

jpeg('../images/kmeans.histogram.jpg')
bp <- barplot(kmeans.clusters.table, col = plot.colors, xlab = "clusters", ylab = "number of elements", main = "histogram - kmeans")
text(bp, 0, kmeans.clusters.table, cex = 1, pos = 3)
dev.off()

jpeg('../images/kmeans.clustering.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = kmeans.clusters$cluster)) + geom_point() + ggtitle("K-means clustering results")
dev.off()
