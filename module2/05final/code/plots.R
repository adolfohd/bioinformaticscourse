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

#####################
jpeg('../images/kmeans.histogram.structural.jpg')
bp <- barplot(kmeans.clusters.table.structural, col = plot.colors, xlab = "clusters", ylab = "number of elements", 
              main = "Structural data: k-means histrogram")
text(bp, 0, kmeans.clusters.table.structural, cex = 1, pos = 3)
dev.off()

jpeg('../images/kmeans.histogram.reduced.jpg')
bp <- barplot(kmeans.clusters.table.reduced, col = plot.colors, xlab = "clusters", ylab = "number of elements", 
              main = "Reduced data: k-means histrogram")
text(bp, 0, kmeans.clusters.table.reduced, cex = 1, pos = 3)
dev.off()

jpeg('../images/kmeans.clustering.structural.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = kmeans.clusters.structural$cluster)
       ) + geom_point() + ggtitle("Structural data: K-means clustering")
dev.off()

jpeg('../images/kmeans.clustering.reduced.jpg')
ggplot(pathways.scores, aes(RC1, RC2, color = kmeans.clusters.reduced$cluster)
) + geom_point() + ggtitle("Reduced data: K-means clustering")
dev.off()
