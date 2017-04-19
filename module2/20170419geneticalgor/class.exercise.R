setwd("~/code/doc/bioinformaticscourse/module2/20170419geneticalgor")

rm(list = ls())

evaluate.chromosome.fitness <- function(chromosome, phenotype){
  profit <- sum(as.numeric(phenotype["value",]) * as.numeric(chromosome))
  total.weight <- sum(as.numeric(phenotype["value",]) * as.numeric(chromosome))
  
  return(profit * (total.weight <=15 ))
}


element <- c("A","P", "L", "F", "C")
weight <- c(2,1,4,1,3)
value <- c(10,5,1,7,3)

phenotype <- matrix(c(element, weight, value), 3,  byrow = T)
rownames(phenotype) <- c("element", "weight", "value")
print(phenotype)

chromosomes <- list()

chromosomes[[1]] <- c(0,1,0,0,1)
chromosomes[[2]] <- c(1,0,1,0,1)
chromosomes[[3]] <- c(0,1,0,1,1)
chromosomes[[4]] <- c(1,1,1,0,0)
chromosomes[[5]] <- c(0,1,1,0,0)

names(chromosomes) <- c(toupper(letters[1:length(chromosomes)]))
##### Evaluate all fitnesses

fitnesses <- array()
for (i in 1:length(chromosomes)){
  fitnesses[i] <- evaluate.chromosome.fitness(chromosomes[[i]], phenotype)
}
names(fitnesses) <- names(chromosomes)
print(fitnesses)

##### Selection of parents

#proportion.of.parents <- 0.5

#sample(fitnesses, , replace = F)
#next.generation[i]

parents <- names(fitnesses)[order(fitnesses, decreasing=TRUE)[1:2]] # take the best 2 parents




##### Crossover
children <- list()
children[[1]] <- c(chromosomes[[parents[1]]][1:2], chromosomes[[parents[2]]][3:5])
children[[2]] <- c(chromosomes[[parents[2]]][1:2], chromosomes[[parents[1]]][3:5])

next.generation <- chromosomes
next.generation[[parents[1]]] <- children[[1]]
next.generation[[parents[2]]] <- children[[2]]

print(chromosomes)
print(next.generation)

##### Mutation