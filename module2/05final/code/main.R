
setwd("~/code/doc/bioinformaticscourse/module2/05final/code")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")

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
hierarch.clusters.complete <- my.hclust(pathways.values)

hierarch.cluster.dendrogram.reduced <- as.dendrogram(hierarch.clusters.reduced)
hierarch.cluster.dendrogram.reduced <- color_branches(hierarch.cluster.dendrogram.reduced, k = k)
hierarch.cluster.dendrogram.complete <- as.dendrogram(hierarch.clusters.complete)
hierarch.cluster.dendrogram.complete <- color_branches(hierarch.cluster.dendrogram.complete, k = k)

hierarch.clusters.cut.reduced <- cutree(hierarch.cluster.dendrogram.reduced, k)
hierarch.clusters.cut.complete <- cutree(hierarch.cluster.dendrogram.complete, k)

hierarch.clusters.table.reduced <- table(hierarch.clusters.cut.reduced)
hierarch.clusters.table.complete <- table(hierarch.clusters.cut.complete)

# K-means
kmeans.clusters.reduced <- kmeans(pathways.scores, centers = 4, nstart = 30)
kmeans.clusters.reduced$cluster <- as.factor(kmeans.clusters.reduced$cluster)
kmeans.clusters.complete <- kmeans(pathways.values, centers = 4, nstart = 30)
kmeans.clusters.complete$cluster <- as.factor(kmeans.clusters.complete$cluster)

kmeans.clusters.table.reduced <- table(kmeans.clusters.reduced$cluster)
kmeans.clusters.table.complete <- table(kmeans.clusters.complete$cluster)

# Jaccard coefficient

my.jaccard.coefficient <- function(labels1, cluster1, labels2, cluster2){
  labels1.comparisson <- labels1 == cluster1
  labels2.comparisson <- labels2 == cluster2
  return (sum(labels1.comparisson & labels2.comparisson ) / sum(labels1.comparisson | labels2.comparisson ))
}

my.jaccard.matrix <- function(labels1, labels2, k){
  labels1 <- as.array(labels1)
  labels2 <- as.array(labels2)
  jaccard.hierarchical  <- matrix(nrow = k, ncol = k)
  max.coeffs.positions  <- array()
  max.coeffs            <- array()
  for (i in 1:k){
    max.coeffs.positions[i] <- 0
    max.coeffs[i] <- 0
    for(j in 1:k){
      jaccard.hierarchical[i,j] <- my.jaccard.coefficient(labels1, i, labels2, j)
      if (jaccard.hierarchical[i,j] > max.coeffs[i]){
        max.coeffs[i]           <- jaccard.hierarchical[i,j]
        max.coeffs.positions[i] <- j
      }
    }
  }
  return(list(jaccard.hierarchical, max.coeffs, max.coeffs.positions))
}

jaccard.matrix.hierarchical <- my.jaccard.matrix(hierarch.clusters.cut.complete, kmeans.clusters.complete$cluster, k)
print(jaccard.matrix.hierarchical)

jaccard.matrix.kmeans <- my.jaccard.matrix(hierarch.clusters.cut.reduced, kmeans.clusters.reduced$cluster, k)
print(jaccard.matrix.kmeans)
# ---- PLOTS -----
source("plots.R")
