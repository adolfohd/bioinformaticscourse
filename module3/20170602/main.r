rm(list = ls())
setwd("/media/fito/Windows/Users/fitoh/Documents/code/doc/bioinformaticscourse/module3/20170602")

blast.result <- read.csv("data/resultblast_withcolumnames", sep = "\t")
source("get.interactions.dataset.r")
n.interactions.dataset <- get.interactions.dataset(blast.result)

# Order
n.interactions.dataset <- n.interactions.dataset[order(n.interactions.dataset$interactions, decreasing = T),]

# Log
n.interactions.dataset$log.interactions     <- as.numeric(log(as.numeric(n.interactions.dataset$interactions)))
n.interactions.dataset$log.length           <- as.numeric(log(as.numeric(n.interactions.dataset$mean.query.length)))
n.interactions.dataset$log.bitscore         <- as.numeric(log(as.numeric(n.interactions.dataset$mean.bitscore)))

# All
source("plot.interactions.histogram.r")
#plot.interactions.histogram(n.interactions.dataset, 25)

source("plot.length.histogram.r")
#plot.length.histogram(n.interactions.dataset, 25)

#plot(x=n.interactions.dataset$interactions, n.interactions.dataset$mean.query.length, log = "xy")

k <- 5
d<-n.interactions.dataset[,c(8,9)]
d <- as.data.frame(d)
kmeans.clustering <- kmeans(d, centers = k)

# n.interactions.dataset <- n.interactions.dataset[order(kmeans.clustering$cluster, decreasing = T),]


# plot(kmeans.clustering)
# Extract nodes with only one connection
# plot(x = n.interactions.dataset$)
library("colorspace")
jpeg(paste("images/clusters/clustering.jpg"))
plot(n.interactions.dataset[,c(1,4)], col= rainbow_hcl(k)[kmeans.clustering$cluster], log = "xy",
     main = "k-means clustering",
     xlab = ("Number of interactions"),
     ylab = ("Average Length")
)
legend("topright", pch=16, col=rainbow_hcl(k),
       legend=unique(kmeans.clustering$cluster))
dev.off()

# points(kmeans.clustering$centers, col=1:5,pch=8,cex=1)

n.items <- array(NA, k)
for (i in 1:k){
  n.items[i] <- sum(kmeans.clustering$cluster==i)
  jpeg(paste("images/clusters/cluster", i, "_", n.items[i] , "items.jpg"))
  plot(n.interactions.dataset[, c(1,4)], 
       col = rainbow_hcl(2)[((kmeans.clustering$cluster == i) + 1)], 
       log = "xy",
       main = (paste("Cluster", i, ": ", n.items[i] , " items")),
       xlab = ("Number of interactions"),
       ylab = ("Average Length")
  )
  dev.off()
}

clusters.ordered.by.nitems <- order(n.items, decreasing = T)
n.items <- n.items[order(n.items, decreasing = T)]

for (i in 1:k){
  print(paste("The cluster ", clusters.ordered.by.nitems[i], "has ", n.items[i], " items"))
}

best.cluster <- readline("Please analize the plots in images/clusters/, and input the label of one chosen cluster : \n") 

filtered.dataset <- n.interactions.dataset[kmeans.clustering$cluster==best.cluster,]

write.csv(filtered.dataset, "filtered.dataset.csv")


