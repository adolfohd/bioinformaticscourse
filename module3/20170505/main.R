# load("/media/fito/Windows/Users/fitoh/Documents/code/doc/bioinformaticscourse/module3/20170505/.RData")
setwd("~/code/doc/bioinformaticscourse/module3/20170505")
rm(list = ls())
# source("~/code/fito_lib/r/pkgtest.R")

# rna.q           <- readLines("../data/humangenome/rna.q", sep = "\t")
rna.q           <- read.csv("../data/humangenome/rna.q", sep = "\t")
# rna.q           <- readLines("rna.q")
seq.gene.md           <- read.csv("../data/humangenome/seq_gene.md", sep = "\t")
# seq.gene.md           <- read.csv("seq_gene.md", sep = "\t")


non.coding.rna.rows  <- which(grepl("non-coding RNA", rna.q$defline, fixed = T))

'
library(stringr)
gene.ids <- str_extract(rna.q[non.coding.rna.rows,]$gene_id, "(?<=GeneID:)\\d+")
gene.ids <- rna.q[non.coding.rna.rows,]$gene_id
gene.info <- rna.q[non.coding.rna.rows,]$defline
'

non.coding.rna <- rna.q[non.coding.rna.rows, c(1,5,7,10)]

ncNRA <- non.coding.rna[ order(non.coding.rna[,2], -non.coding.rna[,3]), ]

unique.ncNRA.rows <- array(NA, length(unique(ncNRA$gene_id)))

unique.ncNRA.rows[1] <- 1
unique.index <- 1
last.unique <- ncNRA$gene_id[unique.index]

gene.ids <- as.character(ncNRA$gene_id)
last.unique <- gene.ids[unique.index]
for (i in 2:nrow(ncNRA)){
  if (gene.ids[i] != last.unique){
    unique.index <- unique.index + 1
    unique.ncNRA.rows[unique.index] <- i
    last.unique <- ncNRA$gene_id[i]
  }
}
ncNRA.unique <- ncNRA[unique.ncNRA.rows,]

seq.gene.md.relevant <- seq.gene.md[
  grepl("^GENE", seq.gene.md$feature_type) &
    grepl("^GRCh38.", seq.gene.md$group_label),] 

noncodingrna.geneinfoindexes <- array()
pb <- txtProgressBar(min = 0, max = nrow(non.coding.rna), style = 3)

# test <- array()
for (i in 1:nrow(non.coding.rna)){
# for (i in 1:2){
    query.string <- paste(ncNRA.unique$gene_id[i], "$", sep="")
    # query.string <- as.character(ncNRA.unique$gene_id[i])
    # print(query.string)
    found.gene <-    which(grepl(    query.string,      seq.gene.md.relevant$feature_id))
    # print(found.gene)
  if (length(found.gene)!=0){
    noncodingrna.geneinfoindexes[i] <- found.gene
    # test[i] <- found.gene
    # print(noncodingrna.geneinfoindexes[i])
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
Sys.sleep(1)
close(pb)

# Let's generate the dataset we are interested in: "old" genome sequence information of noncoding genes. 
noncodingrna.geneinfo <- seq.gene.md.relevant[noncodingrna.geneinfoindexes,]

write.csv(noncodingrna.geneinfo, "ncNRA.seqgene.csv")


#####################################################

noncodingrna.geneinfo <- read.csv("ncNRA.seqgene.csv")

# Now we need to take only the distances

positioninfo.columns = c("feature_id", "chromosome", "chr_start", 
                         "chr_stop", "chr_orient", "contig", "feature_name")
noncodingrna.positioninfo <- noncodingrna.geneinfo[,positioninfo.columns]

# Then we add a column containing the length of those genes

noncodinggenes.length = as.numeric(noncodingrna.positioninfo[,"chr_stop"]) - 
                        as.numeric(noncodingrna.positioninfo[,"chr_start"])
# colnames(noncodinggenes.length) <- "chr_length"
noncodingrna.positioninfo$chr_length <- noncodinggenes.length
head(noncodingrna.positioninfo)

# Since there are columns filled with NA values:

# tail(noncodingrna.positioninfo[-is.na(noncodingrna.positioninfo$chr_length),])

# Let's take them away
# noncodingrna.positioninfo <- noncodingrna.positioninfo[!is.na(noncodingrna.positioninfo),]

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

# Assure that the chromosome column is taken as factor
noncodingrna.positioninfo$chromosome <- as.factor(noncodingrna.positioninfo$chromosome)
order.levels <- function(factors){
  indexes  <- array()
  for (i in 1:22){ # First we look for numbers from 1 to 22
    found <- which(factors==i)
    if (i > 1)
      indexes <- c(indexes, found)
    else
      indexes <- found
  }
  text.factors <- c("X", "Y", "Un") # then, for the other chromosomes
  for (i in text.factors){
    indexes <- c(indexes, which(factors==i))
  }
  return (indexes)
}
ordered.chromosome.indexes <- order.levels(levels(noncodingrna.positioninfo$chromosome))
print(ordered.chromosome.indexes)
ordered.chr.levels <- levels(noncodingrna.positioninfo$chromosome)[ordered.chromosome.indexes]
noncodingrna.positioninfo$chromosome <- factor(noncodingrna.positioninfo$chromosome, ordered.chr.levels)


# lets divide the range of lengths in N segments so we can plot a histogram
n.intervals <- 100
maximum.length  <- max(noncodingrna.positioninfo$chr_length, na.rm = T)
minimum.length  <- min(noncodingrna.positioninfo$chr_length, na.rm = T)
length.range <- maximum.length - minimum.length
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
maximum.length <- 30000
minimum.length <- 0

length.range <- maximum.length - minimum.length
step.size <- length.range/n.intervals



library(ggplot2) #load the ggplot2 graph package
breaks.plot = seq(minimum.length, maximum.length,by=step.size)
outside.members   <- sum(noncodingrna.positioninfo$chr_length >  maximum.length, na.rm = T)
inside.members    <- sum(noncodingrna.positioninfo$chr_length <= maximum.length, na.rm = T)
pkgTest("grob")

pkgTest("scales")
require(scales)
n <- 25
#plot.colors <- hue_pal(h = c(0, 360) + 15, 
#                               c = 100, l = 65, 
#                               h.start = 0, direction = 1)(n)[order(sample(1:n, n))]
plot.colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
"#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
"#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
"#8A7C64")

# jpeg('images/histogram.jpg', width = 8, height = 8, units = 'in', res = 300)
svg('images/histogram.svg')
ggplot(noncodingrna.positioninfo, aes(x=chr_length, fill=factor(chromosome))) + 
  geom_histogram(
    breaks=breaks.plot,
    colour='black') +
  guides(fill=guide_legend(title="Chromosome")) +
  # scale_fill_hue(c=70, l=70) +
  scale_fill_manual(values = plot.colors, drop = F) +
  stat_bin(breaks = breaks.plot, geom="text", colour="white", size=2, 
           aes(label=factor(chromosome), fill=factor(chromosome), y=..count.. - .2), drop = T) +
  scale_x_continuous(breaks=round(breaks.plot,0)) +
  theme(axis.text.x = element_text(angle=90)) + xlab("Non-coding RNA length") + ylab("Count of gene lengths") +
  annotate("text", x = 0.8*maximum.length,  y=300, label = paste(
    "Members in range = \n", inside.members, " ncRNA \nMembers outside range = \n", outside.members, " ncRNA", sep = ""), size=3)
dev.off()
