# load("/media/fito/Windows/Users/fitoh/Documents/code/doc/bioinformaticscourse/module3/20170505/.RData")

rm(list = ls())
# source("~/code/fito_lib/r/pkgtest.R")

# rna.q           <- readLines("../data/humangenome/rna.q")
rna.q           <- readLines("rna.q")
# seq.gene.md           <- read.csv("../data/humangenome/seq_gene.md", sep = "\t")
seq.gene.md           <- read.csv("seq_gene.md", sep = "\t")


non.coding.rna.brute  <- rna.q[grepl("non-coding RNA", rna.q, fixed = T)]

library(stringr)
gene.ids <- str_extract(non.coding.rna.brute, "(?<=GeneID:)\\d+")
gene.info <- str_extract(non.coding.rna.brute, "([^\t]*){1}$")

non.coding.rna <- as.data.frame(matrix(c(gene.ids, gene.info), ncol = 2, byrow = F))
colnames(non.coding.rna) <- c("gene.id", "gene.info")
head(non.coding.rna)

noncodingrna.geneinfoindexes <- array()
pb <- txtProgressBar(min = 0, max = nrow(non.coding.rna), style = 3)

seq.gene.md.relevant <- seq.gene.md[
    grepl("^GENE", seq.gene.md$feature_type) &
    grepl("^GRCh38.", seq.gene.md$group_label),] 


test <- array()
# for (i in 1:nrow(non.coding.rna)){
for (i in 1:2){
    query.string <- paste("GeneID:", non.coding.rna$gene.id[i], "$", sep="")
    print(query.string)
    found.gene <-    which(grepl(    query.string,      seq.gene.md.relevant$feature_id))
    print(found.gene)
  if (length(found.gene)!=0){
    # noncodingrna.geneinfoindexes[i] <- found.gene
    test[i] <- found.gene
    # print(noncodingrna.geneinfoindexes[i])
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
Sys.sleep(1)
close(pb)

# Let's generate the dataset we are interested in: "old" genome sequence information of noncoding genes. 
noncodingrna.geneinfo <- seq.gene.md.relevant[noncodingrna.geneinfoindexes,]


# Now we need to take only the distances

positioninfo.columns = c("feature_id", "chromosome", "chr_start", 
                         "chr_stop", "chr_orient", "contig", "feature_name")
noncodingrna.positioninfo <- noncodingrna.geneinfo[,positioninfo.columns]

# Then we add a column containing the length of those genes

noncodinggenes.length = as.numeric(noncodingrna.positioninfo[,"chr_stop"]) - 
                        as.numeric(noncodingrna.positioninfo[,"chr_start"])
colnames(noncodinggenes.length) <- "chr_length"
noncodingrna.positioninfo$chr_length <- noncodinggenes.length
head(noncodingrna.positioninfo)

# lets divide the range of lengths in N segments so we can plot a histogram
n.intervals <- 100
max_length  <- max(noncodingrna.positioninfo$chr_length, na.rm = T)
min_length  <- min(noncodingrna.positioninfo$chr_length, na.rm = T)
length.range <- max_length - min_length
interval    <- length.range/n.intervals

# count values in every interval

histogram <- array()
x.axis <- array()
for (i in 1:nrow(noncodingrna.positioninfo)){
  current.position <- noncodingrna.positioninfo$chr_length[i]
  target.histogram.value <- 1 + 
                            round(
                              (current.position - min_length)/length.range * n.intervals
                              , 0)
  if (is.na(histogram[target.histogram.value]))
    histogram[target.histogram.value] <- 1
  else
    histogram[target.histogram.value] <- histogram[target.histogram.value] + 1
}

for (i in 1:n.intervals+1)
  x.axis[i] <- min_length + interval*(i-1)

plot(x.axis[1:20], histogram[1:20], type = "l")


# This showed us that the majority of non-coding genes are relatively short, so we are only 
# going to compute the histogram in lengths in steps of only 100, up to 10,000.

histogram <- array()
x.axis <- array()
n.intervals <- 100
max_length  <- 10000
min_length  <- 0
length.range <- max_length - min_length
interval    <- length.range/n.intervals
for (i in 1:nrow(noncodingrna.positioninfo)){
  current.length <- noncodingrna.positioninfo$chr_length[i]
  if (current.length < max_length && !is.na(current.length)){
    target.histogram.value <- 1 + 
      floor(
        (current.length - min_length)/length.range * n.intervals)
    if (is.na(histogram[target.histogram.value]))
      histogram[target.histogram.value] <- 1
    else
      histogram[target.histogram.value] <- histogram[target.histogram.value] + 1

  }
}

for (i in 1:n.intervals+1)
  x.axis[i] <- min_length + interval*(i-1)

plot(x.axis[1:20], histogram[1:20], type = "l")




###


# plot(hist(noncodingrna.positioninfo$chr_length, xlim = c(0,1000)))
# plot(density(noncodingrna.positioninfo$chr_length, na.rm = T, from = 0, to = 100000, n = 100))

data <- as.data.frame(noncodingrna.positioninfo$chr_length)

library(ggplot2) #load the ggplot2 graph package

ggplot(noncodingrna.positioninfo, aes(x=chr_length)) + 
  geom_histogram(breaks=seq(0,100000,by=100),colour='black') + 
  xlab("Average Frequency") + ylab("Count of USER.ID") +
  scale_fill_manual("Group", breaks = c("1","2","3"), values = c("grey30","grey50", "grey70")) +
  theme_bw()

ggplot(noncodingrna.positioninfo, aes(x=chr_length, fill=factor(chromosome))) + 
  geom_histogram(bresaks=seq(0,100,by=100),colour='black') + 
  xlab("Average Frequency") + ylab("Count of gene lengths") +
  scale_fill_manual("Chromosome", breaks = 1:39, values = colors(1:39)) +
  theme_bw()


save.image()

