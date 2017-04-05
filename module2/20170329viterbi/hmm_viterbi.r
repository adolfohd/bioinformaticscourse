rm(list = ls())
setwd('~/code/doc/bioinformaticscourse/module2/20170329viterbi')
nucleotides <- c("A","C", "G", "T")
states <- c("H", "L")

# HMM parameters

a.high <- c(0.9, 0.1)
a.low <- c(0.3, 0.7)
a   <- matrix(c(a.high, a.low), 2,2 , byrow = T)
rownames(a) <- states
colnames(a) <- states

pi  <- as.matrix( c(0.5, 0.5) )
rownames(pi) <- states

e.high  <- c(0.1, 0.4, 0.4, 0.1)
e.low   <- c(0.4, 0.1, 0.1, 0.4)
e       <- matrix(c(e.high , e.low), 2, 4, byrow= T)
rownames(e) <- states
colnames(e) <- nucleotides

# sequence generation

N=25

source("cg.rich.sequence.generator.r")
generated.sequences <- cg.rich.sequence.generator(a, e, pi, N)
measures.sequence <- generated.sequences[[1]]
states.sequence <- generated.sequences[[2]]
cat("\n", measures.sequence)
cat("\n", states.sequence)

# Viterbi
path <- character(length(measures.sequence))
source("Vk.r")
Vk.H <- Vk("H", N)
Vk.L <- Vk("L", N)
print(path)

