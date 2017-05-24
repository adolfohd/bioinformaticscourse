# load("/media/fito/Windows/Users/fitoh/Documents/code/doc/bioinformaticscourse/module3/20170505/.RData")
setwd("~/code/doc/bioinformaticscourse/module3/20170505")
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

# Since there are columns filled with NA values:

noncodingrna.positioninfo[is.na(noncodingrna.positioninfo$chromosome),]

# Let's take them away
noncodingrna.positioninfo <- noncodingrna.positioninfo[!is.na(noncodingrna.positioninfo),]

# In order to be able to show highlight which genes belong to each chromosomo, we're going to
# clean the chromosome labels, since there are some duplicates in its unique values:
uni <- unique(noncodingrna.positioninfo$chromosome)
print(uni)

only.numeric.chromosome.indexes <- uni[grepl("^\\d+$", uni)]
for (i in only.numeric.chromosome.indexes){
  noncodingrna.positioninfo$chromosome[grepl(paste("^",i,"\\|", sep = ""), noncodingrna.positioninfo$chromosome)] <- i
}
# Assign the label "Un" to all choromosomes with names starting with "Un"
noncodingrna.positioninfo$chromosome[grepl("Un", noncodingrna.positioninfo$chromosome)]

noncodingrna.positioninfo$chromosome <- as.character(noncodingrna.positioninfo$chromosome)
noncodingrna.positioninfo$chromosome[grepl("Un", noncodingrna.positioninfo$chromosome)] <- "Un"
unique(noncodingrna.positioninfo$chromosome)

# lets divide the range of lengths in N segments so we can plot a histogram
n.intervals <- 100
maximum.length  <- max(noncodingrna.positioninfo$chr_length, na.rm = T)
minimum.length  <- min(noncodingrna.positioninfo$chr_length, na.rm = T)
length.range <- max_length - min_length
step.size    <- length.range/n.intervals

library(ggplot2) #load the ggplot2 graph package
ggplot(noncodingrna.positioninfo, aes(x=chr_length, fill=factor(chromosome))) + 
  geom_histogram(
    breaks=seq(minimum.length, maximum.length,by=step.size),
    colour='black') + 
  xlab("Average Frequency") + ylab("Count of gene lengths") +
  # scale_fill_manual("Chromosome", breaks = 1:39, values = colors(1:39)) +
  theme_bw()

# This showed us that the majority of non-coding genes are relatively short, so we are only 
# going to compute the histogram in lengths in steps of only 100, up to 10,000.

n.intervals <- 30
maximum.length <- 20000
minimum.length <- 0

length.range <- maximum.length - minimum.length
step.size <- length.range/n.intervals

# Assure that the chromosome column is taken as factor
noncodingrna.positioninfo$chromosome <- as.factor(noncodingrna.positioninfo$chromosome)

ordered.chromosome.indexes =c(1,12,
  16:22,
  2:11,
  13:15,24, 25, 23)
ordered.chr.levels <- levels(noncodingrna.positioninfo$chromosome)[ordered.chromosome.indexes]
noncodingrna.positioninfo$chromosome <- factor(noncodingrna.positioninfo$chromosome, ordered.chr.levels)

jpeg('images/histogram.jpg')
library(ggplot2) #load the ggplot2 graph package
breaks.plot = seq(minimum.length, maximum.length,by=step.size)
ggplot(noncodingrna.positioninfo, aes(x=chr_length, fill=factor(chromosome))) + 
  geom_histogram(
    breaks=breaks.plot,
    colour='black') +
  stat_bin(breaks = breaks.plot, geom="text", colour="white", size=2, 
           aes(label=factor(chromosome), fill=factor(chromosome), y=..count.. - .2), drop = T) +
  scale_x_continuous(breaks=round(breaks.plot,0)) +
  theme(axis.text.x = element_text(angle=90))
          #axis.text.y = element_text(face="bold", color="#993333", 
           #                          size=14, angle=45))
dev.off()
