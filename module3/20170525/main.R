rm(list = ls())

setwd("~/code/doc/bioinformaticscourse/module3/20170525")
reference.lines <- readLines("reference/humanprot.nwb")

filename <- "reference/blastresumen.cytos.txt"
# filename <- "reference/blastresumen.cytos.sample.txt"

source("csv2nwb.r")
lines.tosave <- to.nwb(filename)
save.lines.tofile(lines = lines.tosave, filename = filename)
