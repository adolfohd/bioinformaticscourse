setwd("~/code/doc/bioinformaticscourse/module3/20170505/")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")

rna.q           <- readLines("../data/humangenome/rna.q")
# seq.gene.md     <- readLines("../data/humangenome/seq_gene.md")
seq.gene.md           <- read.csv("../data/humangenome/seq_gene.md", sep = "\t")


non.coding.rna.brute  <- rna.q[grepl("non-coding RNA", rna.q, fixed = T)]

library(stringr)
gene.ids <- str_extract(non.coding.rna.brute, "(?<=GeneID:)\\d+")
gene.info <- str_extract(non.coding.rna.brute, "([^\t]*){1}$")

non.coding.rna <- as.data.frame(matrix(c(gene.ids, gene.info), ncol = 2, byrow = F))
colnames(non.coding.rna) <- c("gene.id", "gene.info")
head(non.coding.rna)

rna.gene.info <- list()
pb <- txtProgressBar(min = 0, max = nrow(non.coding.rna), style = 3)
# for (i in 1:nrow(non.coding.rna)){
for (i in 1:2){
  rna.gene.info[[i]] <- seq.gene.md[grepl(paste("GeneID:", gene.ids[i], "$", sep=""), seq.gene.md$feature_id),] 
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
Sys.sleep(1)
close(pb)
save.image(file = "data/.Rdata")

