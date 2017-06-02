
plot.interactions.histogram <- function(n.interactions.dataset, n.intervals){
  maximum.xaxis <- max(n.interactions.dataset$interactions)
  minimum.xaxis <- min(n.interactions.dataset$interactions)
  xaxis.range <- maximum.xaxis - minimum.xaxis
  step.size <- xaxis.range/n.intervals
  
  breaks.plot = seq(minimum.xaxis, maximum.xaxis,by=step.size)
  
  library(ggplot2) #load the ggplot2 graph package
  ggplot(n.interactions.dataset, aes(x=interactions)) + 
    geom_histogram(
      breaks=breaks.plot,
      colour='black') +
    stat_bin(breaks = breaks.plot, binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) +
    scale_x_continuous(breaks=round(breaks.plot,0)) 
}

