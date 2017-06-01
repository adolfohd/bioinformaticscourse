to.nwb <- function(filename){
  lines <- c("#Converted to NWB format by Adolfo Hoyos",
             "*Nodes",
             "id*int	label*string")
  dataset <- read.csv(filename, sep = "\t")
  colnames(dataset)[1] <- "node1"
  colnames(dataset)[2] <- "node2"
  colnames(dataset)[3] <- "weight"
  colnames(dataset)[4] <- "repeat"
  # erase comment and empty lines
  
  if(sum(levels(dataset$node1) != levels(dataset$node2)))
    return (NA)
  else
    unique.levels <- levels(dataset$node1)
  lines[2] <- paste(lines[2], length(unique.levels))
  i <- 1
  for (unique.node in unique.levels){
    lines <- c(lines, paste(i, ' "',unique.node, '"', sep = "\t"))
    i <- i+1
  }
  lines <- c(lines, paste("*UndirectedEdges", nrow(dataset)))
  lines <- c(lines, paste("source*int target*int ", colnames(dataset)[3], "*float", sep = ""))
  
  print(paste(nrow(dataset), "rows read."))
  
  node.relationship.columns <- c(
    as.character(as.numeric(dataset$node1)),
    as.character(as.numeric(dataset$node2)),
    as.character(dataset$weight)
  )
  pb <- txtProgressBar(min = 0, max = nrow(dataset), style = 3)
  morelines <- array(NA, nrow(dataset))
  one.percent <- ceiling(nrow(dataset)/100)
  for (i in 1: nrow(dataset)){
    morelines[i] <- paste(
      node.relationship.columns[i],                  # Source node indexes
      node.relationship.columns[i + nrow(dataset)],  # Target node indexes
      node.relationship.columns[i + 2*nrow(dataset)] # Weights
    , sep = "\t")
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
