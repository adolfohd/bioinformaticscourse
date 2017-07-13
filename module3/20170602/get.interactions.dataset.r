
get.interactions.dataset <- function(blast.result){
  n.interactions.dataset <- data.frame(row.names = as.factor(levels(blast.result$query)))
  n <- nrow(n.interactions.dataset)
  interactions        = as.numeric(array(0, n))
  mean.similarity     = as.numeric(array(0, n))
  mean.expected.value = as.numeric(array(0, n))
  mean.query.length   = as.numeric(array(0, n))
  mean.mismatches     = as.numeric(array(0, n))
  mean.gaps           = as.numeric(array(0, n))
  mean.bitscore       = as.numeric(array(0, n))
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  one.percent <- ceiling(n/100)
  for ( i in 1:n){
    same.level <- (blast.result$query == levels(blast.result$query)[i])
    interactions[i]        <- sum(same.level)
    mean.similarity[i]     <- sum(as.numeric(blast.result$`%id`           [same.level])) / interactions[i]
    mean.expected.value[i] <- sum(as.numeric(blast.result$E.value          [same.level])) / interactions[i]
    mean.query.length[i]   <- sum(as.numeric(blast.result$alignment.length [same.level])) / interactions[i]
    mean.mismatches[i]     <- sum(as.numeric(blast.result$mismatches       [same.level])) / interactions[i]
    mean.gaps[i]           <- sum(as.numeric(blast.result$gap.openings     [same.level])) / interactions[i]
    mean.bitscore[i]       <- sum(as.numeric(blast.result$bit.score        [same.level])) / interactions[i]
    if (! (i %% one.percent)) # limit updates of progress bar every time a percent of the job is complete
      setTxtProgressBar(pb, i)
  }
  close(pb)
  
  n.interactions.dataset$interactions        <- interactions
  n.interactions.dataset$mean.similarity     <- mean.similarity
  n.interactions.dataset$mean.expected.value <- mean.expected.value
  
  n.interactions.dataset$mean.query.length   <- mean.query.length
  n.interactions.dataset$mean.mismatches     <- mean.mismatches
  n.interactions.dataset$mean.gaps           <- mean.gaps
  n.interactions.dataset$mean.bitscore       <- mean.bitscore

  return(n.interactions.dataset)
}
