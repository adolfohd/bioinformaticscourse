to.nwb <- function(filename){
  ## region for all node properties 
  
  lines <- c("#Converted to NWB format by Adolfo Hoyos",
             "*Nodes",
             "id*int	label*string	height*float	width*float")
  dataset <- read.csv(filename, sep = "\t")

  # colnames(dataset)[1] <- "node1"
  # colnames(dataset)[2] <- "node2"
  # colnames(dataset)[3] <- "weight"
  # colnames(dataset)[4] <- "repeat"

  # erase comment and empty lines
  unique.levels <- union(levels(dataset[,1]), levels(dataset[,2]))
  n = nrow(dataset)
  lines[2] <- paste(lines[2], length(unique.levels))
  i <- 1
  for (unique.node in unique.levels){
    w <- sum((dataset[,1]==unique.node)| (dataset[,2] == unique.node))/n*100
    lines <- c(lines, paste(i, '\t"',unique.node, '"\t', w, "\t", w, sep = ""))
    i <- i+1
  }
  
  ## Region for interactions between pairs of nodes
  
  
  
  lines <- c(lines, paste("*UndirectedEdges", nrow(dataset)))
  lines <- c(lines, paste("source*int","target*int", paste(colnames(dataset)[3:ncol(dataset)], "*float", sep = "", collapse = "\t"), sep = "\t"))
  
  print(paste(nrow(dataset), "rows read."))
  node.relationship.columns <- c(
    as.character(as.numeric(dataset[,1])),
    as.character(as.numeric(dataset[,2]))#,
    #as.character(dataset$weight)
  )
  
  for (i in 3:ncol(dataset)){
    node.relationship.columns <- c(node.relationship.columns, as.character(dataset[,i]))
  }
  
  pb <- txtProgressBar(min = 0, max = nrow(dataset), style = 3)
  morelines <- array(NA, nrow(dataset))
  one.percent <- ceiling(nrow(dataset)/100)
  for (i in 1: nrow(dataset)){
    morelines[i] <- node.relationship.columns[i]
    for (j in 1:((ncol(dataset))-1)){
    morelines[i] <- paste(morelines[i], node.relationship.columns[i + j*nrow(dataset)], sep = "\t")
     # node.relationship.columns[i],                  # Source node indexes
    #  node.relationship.columns[i + nrow(dataset)],  # Target node indexes
    #  node.relationship.columns[i + 2*nrow(dataset)] # Weights
    #  , sep = "\t")
    }
    if (! (i %% one.percent)) # limit updates of progress bar every time a percent of the job is complete
      setTxtProgressBar(pb, i)
  }
  Sys.sleep(1)
  close(pb)
  print("Processing complete")
  lines <- c(lines, morelines)
  return(lines)
}

save.lines.tofile <- function(lines, filename){
  library(stringr)
  filename.output <- paste(str_extract(filename, "(?<=/)(.*?)(?=\\.)"), ".nwb", sep = "")
  file.connection <-file(filename.output)
  writeLines(lines, filename.output)
  close(file.connection)
  print(paste("File", filename.output, "saved in work directory."))
}
