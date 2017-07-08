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
jpeg("images/histograms/interaction.histogram.jpg")
plot.interactions.histogram(n.interactions.dataset, 25)
dev.off()
source("plot.length.histogram.r")
jpeg("images/histograms/length.histogram.jpg")
plot.length.histogram(n.interactions.dataset, 25)
dev.off()
#plot(x=n.interactions.dataset$interactions, n.interactions.dataset$mean.query.length, log = "xy")

## <!-- Code provided by J.Murillo
# wss<- (nrow(scores[1:199, ]-1)*sum(apply(scores[1:199, ], 2, var)))
# for (i in 2:15)
#   wss[i] <- sum(kmeans(scores[200:397, ], centers = i)$withinss)
# plot(1:15, wss, type ="b", xlab = "Number of cluster", ylab = "withinss groups sum of error squares")
# 
# library(vegan)
# fit <- cascadeKM(scale(scores[1:199,], center =  TRUE, scale = TRUE), 1, 10, iter = 1000)
# plot(fit, sortg = TRUE, grpmts.plot = TRUE)
# calinski.best <- as.numeric(which.max(fit$results[2,]))
# calinski.best -->
d<-n.interactions.dataset[,c(8,9)]
wss <- array(NA, 15)
for (i in 2:15)
    wss[i] <- sum(kmeans(n.interactions.dataset, centers = i)$withinss)
jpeg("images/sumofsquaresVSnumberofgroups.jpg")
  plot(1:15, wss, type ="b", xlab = "Number of clusters", ylab = "withinss groups sum of error squares")
dev.off()

library(vegan)
fit <- cascadeKM(scale(d, center =  TRUE, scale = TRUE), 1, 10, iter = 1000)
jpeg("images/kalinskicriterion.jpg")
  plot(fit, sortg = TRUE, grpmts.plot = TRUE)
dev.off()


k <-6 # Best K accorrding to both methods
d <- as.data.frame(d)
# d <- n.interactions.dataset[,2:8]
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
       legend=1:k)
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

'
# Now we need to look for these IDs in the next source file ???
source.cytos.dataset <- read.csv("../20170525/reference/blastresumen.cytos.txt", sep = "\t")
grepl(source.cytos.dataset$Nodo1[1], rownames(filtered.dataset))
'
###

filtered.dataset <- read.csv("filtered.dataset.csv")

###
library(vegan)
fit <- cascadeKM(scale(filtered.dataset[,-1], center =  TRUE, scale = TRUE), 1, 10, iter = 1000)
jpeg("images/kalinskicriterion_filtered.jpg")
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
dev.off()

k.f <- 20 # Best K accorrding to both methods
d.f <- as.data.frame(filtered.dataset)
kmeans.clustering.filtered <- kmeans(d.f[,2:8], centers = k.f)

n.items <- array(NA, k.f)
for (i in 1:k.f){
  library(colorspace)
  n.items[i] <- sum(kmeans.clustering.filtered$cluster==i)
  jpeg(paste("images/clusters/filtered/cluster", i, "_", n.items[i] , "items.jpg"))
  plot(filtered.dataset[, c(2,5)], 
       col = rainbow_hcl(2)[((kmeans.clustering.filtered$cluster == i) + 1)], 
       log = "xy",
       main = (paste("Cluster", i, ": ", n.items[i] , " items")),
       xlab = ("Number of interactions"),
       ylab = ("Average Length")
  )
  dev.off()
}



## knn

library(class)
k.f <- 20

knn.filtered <- knn()



