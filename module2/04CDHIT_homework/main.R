setwd("~/code/doc/bioinformaticscourse/module2/04CDHIT_homework")
rm(list = ls())


source("/home/fito/code/fito_lib/r/pkgtest.R")
source("cd.hit.lib.R") # all methods used here for CDHIT
## ------------------- MAIN -------------------------- ##

pkgTest("data.table")
pkgTest ("bio3d")

# sequences<-generate.random.dna.sequences(10,7,4) # in case we don't have any sequences to work with
# sequences <- c(sequences, sequences, sequences, sequences) # duplicate just for fun

library(bio3d)
fastaFile = read.fasta ("DiceMini01.fasta")
fasta.sequences <- as.matrix(fastaFile[[2]])
sequences <- array()

for (i in 1:nrow(fasta.sequences)){
  sequences[i] <- paste(fasta.sequences[i,], sep = "", collapse = "")
  sequences[i] <- unlist(strsplit(sequences[i], split='-', fixed=TRUE))[1]
}
names(sequences) <- rownames(fasta.sequences)


#---- the algorithm
ordered.sequences <- order.sequences.by.length(sequences) # this should be inside "cd.hit" but let's keep it here

print("ordered sequences:")
print(ordered.sequences)

threshold <- 0.8
k.mers <- c(2,5) # the "k" values to be evaluated, 2-mers and 5-mers are the default in the CD-HIT algorithm
cd.hit.results <- cd.hit(ordered.sequences, threshold, k.mers)

cluster.sets <- cd.hit.results["cluster.sets"]

print((cluster.sets))
