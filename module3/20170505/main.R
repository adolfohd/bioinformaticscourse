setwd("~/code/doc/bioinformaticscourse/module3/20170505/")
rm(list = ls())
source("~/code/fito_lib/r/pkgtest.R")

seq.gene.md <- as.data.frame(read.table("../data/humangenome/seq_gene.md"))
