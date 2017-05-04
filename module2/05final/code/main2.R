setwd("~/code/doc/bioinformaticscourse/module2/05final/code")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")
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


