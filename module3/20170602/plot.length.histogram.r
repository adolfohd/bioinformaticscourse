
plot.length.histogram <- function(n.interactions.dataset, n.intervals){
  maximum.xaxis <- max(n.interactions.dataset$mean.query.length)
  minimum.xaxis <- min(n.interactions.dataset$mean.query.length)
  xaxis.range <- maximum.xaxis - minimum.xaxis
  step.size <- xaxis.range/n.intervals
  
  breaks.plot = seq(minimum.xaxis, maximum.xaxis,by=step.size)
  
  library(ggplot2) #load the ggplot2 graph package
  ggplot(n.interactions.dataset, aes(x=mean.query.length)) + 
    labs(x = "Mean Query Length", y = "Node Count", title = "Histogram of Mean Query Length of Node Interactions") +
    theme(axis.text.x = element_text(angle=90)) +
    geom_histogram(
      breaks=breaks.plot,
      colour= "orange") +
    stat_bin(breaks = breaks.plot, binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) +
    scale_x_continuous(breaks=round(breaks.plot,0)) 
  
}

