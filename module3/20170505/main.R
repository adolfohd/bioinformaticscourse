setwd("~/code/doc/bioinformaticscourse/module3/20170505/")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")

rna.q           <- readLines("../data/humangenome/rna.q")
seq.gene.md     <- readLines("../data/humangenome/seq_gene.md")

non.coding.rna.brute <- rna.q[grepl("non-coding RNA", rna.q, fixed = T)]

library(stringr)
gene.ids <- str_extract(non.coding.rna.brute, "(?<=GeneID:)\\d+")
gene.info <- str_extract(non.coding.rna.brute, "([^\t]*){1}$")

non.coding.rna <- as.data.frame(matrix(c(gene.ids, gene.info), ncol = 2, byrow = F))
colnames(non.coding.rna) <- c("gene.id", "gene.info")
head(non.coding.rna)

gene.info <- list()
for (gene in gene.ids){
  gene.info <- seq.gene.md[grepl(gene, seq.gene.md, fixed = T)]  
}



