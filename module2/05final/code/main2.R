setwd("~/code/doc/bioinformaticscourse/module2/05final/code")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")
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

avillin.join <- as.data.frame(read.table("../data/avillin-join-lnx.values"))

# 1

N = 1000

datasets <- split(
              avillin.join, 
              cut(
                1:nrow(avillin.join), 
                seq(
                  0, 
                  nrow(avillin.join), 
                  length.out = N
                )
              )
            )

# 2

K = length(datasets)

medoid <- function(d){
  m <- as.array(colMeans(d))
  d.min <- 20
  med<- 0
  n <- nrow(d)
  for (i in 1:n){
    ma <- rbind(as.array(m),d[i,])
    di <- dist(ma)
    if (di < d.min){
      d.min <- di
      med <- i
    }
  }
  return(med)
}

get.medoids <- function(datasets){
  indexes <- matrix()
  values <- matrix(0,0,3)
  colnames(values) <- colnames(datasets[[1]])
  print(values)
  K = length(datasets)
  for (i in 1:K){
    indexes[i] <- medoid(datasets[[i]])  
    values <- rbind(values, datasets[[i]][indexes[i],])
  }
  return(list(indexes, values))
}

medoids <- get.medoids(datasets)

# 3

medoids.dataset <- as.data.frame(medoids[[2]])

summary(avillin.join)
summary(medoids.dataset)

# 4

k <- 4

## Hierarchical clustering
hierarch.clusters             <- my.hclust(medoids.dataset)

hierarch.cluster.dendrogram   <- as.dendrogram(hierarch.clusters)
hierarch.cluster.dendrogram   <- color_branches(hierarch.cluster.dendrogram, k = k)

hierarch.clusters.cut         <- cutree(hierarch.cluster.dendrogram, k)

hierarch.clusters.table       <- table(hierarch.clusters.cut)

## K-means
kmeans.clusters               <- kmeans(medoids.dataset, centers = 4, nstart = 30)
kmeans.clusters$cluster       <- as.factor(kmeans.clusters$cluster)

kmeans.clusters.table         <- table(kmeans.clusters$cluster)


## Plots

### Hierarchichal

plot.colors <- rainbow_hcl(k)
jpeg('../images/02/hierarch.dendrogram.jpg')
plot(hierarch.cluster.dendrogram, main = "Hierarchical clustering")
dev.off()

jpeg('../images/02/hierarch.histogram.jpg')
bp <- barplot(hierarch.clusters.table, col = plot.colors,
              xlab = "clusters", ylab = "number of elements", 
              main = "Hierarchical clustering histogram")
text(bp, 0, hierarch.clusters.table, cex = 1, pos = 3)
dev.off()

jpeg('../images/02/hierarchical.clustering.jpg')
ggplot(medoids.dataset, aes(PC1, PC2, color = as.factor(hierarch.clusters.cut))
) + geom_point() + ggtitle("Hierarchical clustering")
dev.off()

### k-means

jpeg('../images/02/kmeans.histogram.jpg')
bp <- barplot(kmeans.clusters.table, col = plot.colors, xlab = "clusters", ylab = "number of elements", 
              main = "k-means histrogram")
text(bp, 0, kmeans.clusters.table, cex = 1, pos = 3)
dev.off()

jpeg('../images/02/kmeans.clustering.jpg')
ggplot(medoids.dataset, aes(PC1, PC2, color = kmeans.clusters$cluster)
) + geom_point() + ggtitle("k-means clustering")
dev.off()



