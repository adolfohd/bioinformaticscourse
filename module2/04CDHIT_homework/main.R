setwd("~/code/doc/bioinformaticscourse/module2/04CDHIT_homework")
rm(list = ls())


source("/home/fito/code/fito_lib/r/pkgtest.R")
source("cd.hit.lib.R") # all methods used here for CDHIT
## ------------------- MAIN -------------------------- ##

pkgTest("data.table")

sequences<-generate.random.dna.sequences(10,7,4)
sequences <- c(sequences, sequences, sequences, sequences) # duplicate just for fun

print("original sequences are:")
print(sequences)


#---- the algorithm
ordered.sequences <- order.sequences.by.length(sequences) # this should be inside "cd.hit" but let's keep it here

print("ordered sequences:")
print(ordered.sequences)

threshold <- 0.9
k.mers <- c(2,5) # the "k" values to be evaluated, 2-mers and 5-mers are the default in the CD-HIT algorithm
cd.hit.results <- cd.hit(ordered.sequences, threshold, k.mers)

print(cd.hit.results)






