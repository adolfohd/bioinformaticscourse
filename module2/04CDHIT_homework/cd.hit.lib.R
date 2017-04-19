###################################
generate.random.dna.sequences <- function(max.length, min.length, n){
  sequences <- array()
  
  if (max.length%%1!=0 || min.length%%1!=0 || n%%1!=0){
    stop("all parameters must be integer")
  }
  if(max.length<min.length)
    stop("maximum length must not be lower than minimum length")
  
  for (i in 1:n){
    sequence.size <- sample(min.length:max.length,1)
    sequences[i] <- paste(sample(x = c("A", "C", "G", "T"), size = sequence.size, replace = T), sep = "", collapse = "")
  }
  return(sequences)
}
order.sequences.by.length <- function(seq){
  return(seq[order(nchar(seq), seq, decreasing = T)])
}
###################################

source("/home/fito/code/fito_lib/r/pkgtest.R")
generate.kmers.dictionary <- function(alphabet, k){
  pkgTest("data.table")
  library(data.table)
  lookup.table <- do.call(CJ, replicate(k, alphabet, FALSE))
  return(apply(lookup.table,1,paste,collapse=""))
}
###################################

get.minimum.similarities <- function(L, K, p){
  return(ceiling(L - K + 1 - (1 - p) * K * L))
}

kmers.comparisson <- function(seq.a, seq.b, k, threshold){
  lookup.table <- generate.kmers.dictionary(c("A", "C", "G", "T"),k)
  # truncate sequences to the shortest length
  # print(k)
  # print(c("seq.a=", seq.a))
  # print(c("seq.b=", seq.b))
  n <- min(c(nchar(seq.a), nchar(seq.b)))
  seq.a <- substr(seq.a, 1, n)
  seq.b <- substr(seq.b, 1, n)
  seq.a.strings <- array()
  seq.b.strings <- array()
  for(i in 1:n-k+1){
    seq.a.strings[i] <- substr(seq.a,i,i+k-1)
    seq.b.strings[i] <- substr(seq.b,i,i+k-1)
  }
  seq.a.intable <- as.numeric(factor(seq.a.strings, levels = lookup.table))
  seq.b.intable <- as.numeric(factor(seq.b.strings, levels = lookup.table))
  
  # print(c("seq.a.intable:",seq.a.intable))
  # print(c("seq.b.intable", seq.b.intable))
  
  similar.kmers <- sum(seq.a.intable==seq.b.intable)
  
  minimum.similarities <- get.minimum.similarities(n,k,threshold)
  return(similar.kmers >= minimum.similarities)
}

compare.sequences <- function(seq.a, seq.b, kmers, threshold){ # kmers is an array with all the "k" values to be evaluated (eg: 2-mers and 5-mers)
  kmers.similarities <- array()
  for (i in 1:length(kmers)){
    kmers.similarities[i] <- kmers.comparisson(seq.a, seq.b, kmers[i], threshold) # store all evaluated similarities
  }
  # print(c("similarities=", kmers.similarities))
  return(all(kmers.similarities)) # return TRUE only if all similarities are true 
}

cd.hit <- function(sequences, threshold, k.mers){
  N = length(sequences)
  print (c("N=", N))
  clusters    <- array()
  cluster.sets <- list()
  non.evaluated.sequences <- 1:N # Here we store the indexes of the sequences that haven't been evaluated
  
  i <- 1
  sequence.labels <- array()
  while(length(non.evaluated.sequences)>0){ # As long as we still have sequences to evaluate
    print(c("i=", i))
    clusters[i]<- non.evaluated.sequences[1] # store the first non-evaluated sequence as a new cluster
    sequence.labels[clusters[i]] <- clusters[i] # this cluster gives label to the corresponding sequence
    cluster.sets[[i]]<-array()
    cluster.sets[[i]][1] <- names(sequences)[clusters[i]]
    non.evaluated.sequences <- non.evaluated.sequences[-1] # and remove it from the sequences to be evaluated 
    sequences.to.evaluate <- non.evaluated.sequences # this is to avoid the in-loop elimination to affect the loop's indexes
    for (j in sequences.to.evaluate){ # every non-evaluated sequence
      seq.a <- sequences[clusters[i]]
      seq.b <- sequences[j]
      are.similar.sequences = compare.sequences(seq.a, seq.b, k.mers, threshold)
      if(are.similar.sequences){
        print(c("similar!!, j=", j))
        sequence.labels[j] <- clusters[i] # assign the current cluster as label of the currently evaluated sequence
        cluster.sets[[i]] <- c(cluster.sets[[i]], names(sequences)[j]) # add this sequence to the current cluster
        non.evaluated.sequences <- non.evaluated.sequences[ ! non.evaluated.sequences %in% j] # remove this sequence so we don't evaluate it again
      } # don't do anything if evaluated sequences do not belong to the evaluated clusters
    }
    i <- i + 1
  }
  return(list(clusters,sequence.labels, cluster.sets))
}
