
cg.rich.sequence.generator <- function(a, e, pi, N){
  measures.sequence <- character()
  states.sequence   <- character()
  
  state.first       <- sample(states, 1, rep= T, prob = pi) # first state probability given by initial state
  nucleotide.first  <- sample(nucleotides, 1,prob = e[state.first,])
  
  measures.sequence[1]  <- nucleotide.first
  states.sequence[1]    <- state.first
  
  for (i in 2:N){
    state.previous        <- states.sequence[i-1] 
    print(state.previous)
    state.current         <-  sample(states, 1, rep= T, prob = a[state.previous,]) # assign the probability of the previous state
    nucleotide.current    <- sample(nucleotides, 1, prob = e[state.current,])
    
    measures.sequence[i]  <- nucleotide.current
    
    states.sequence[i]    <- state.current
  }
  return( list(measures.sequence, states.sequence))
}
