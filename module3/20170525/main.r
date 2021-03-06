#!/usr/bin/Rscript
rm(list = ls())

setwd("~/code/doc/bioinformaticscourse/module3/20170525")
# reference.lines <- readLines("reference/humanprot.nwb")

args <- commandArgs(TRUE)
filename <- args[1]
# filename <- "reference/blastresumen.cytos.txt"
filename <- "reference/g.genetype.txt"
# filename <- "reference/blastresumen.cytos.sample.txt"

source("to.nwb.r")
lines.tosave <- to.nwb(filename)
save.lines.tofile(lines = lines.tosave, filename = filename)
