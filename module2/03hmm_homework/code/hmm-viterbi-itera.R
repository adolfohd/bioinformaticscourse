#--------------------------------------------------
# Implements an iterative Viterbi algorithm 
#--------------------------------------------------


#--------------------------------------------------
# A Hidden Markov Model of DNA sequence evolution 
#--------------------------------------------------
# The transition matrix and emission matrix for a HMM
nucleotides         <- c("A", "C", "G", "T")   # Define the alphabet of nucleotides
states              <- c("H", "L")             # Define the names of the states Hight (CG-rich) and Low (AT-rich)
# GCrichprobs         <- c(0.9, 0.1)             # Set the probabilities of switching states, where the previous state was "AT-rich"
# ATrichprobs         <- c(0.3, 0.7)             # Set the probabilities of switching states, where the previous state was "GC-rich"
GCrichprobs         <- c(0.9, 0.1)
ATrichprobs         <- c(0.3, 0.7)
#--------------------------------------------------------------
# Transitions
#--------------------------------------------------------------
transitions <- matrix(c(GCrichprobs, ATrichprobs), 2, 2, byrow = TRUE) # Create a 2 x 2 matrix
rownames(transitions) <- states
colnames(transitions) <- states

cat (">>> The transition matriz for the High GC and Low GC:\n")
print (transitions)

#--------------------------------------------------------------
# Emissions
#--------------------------------------------------------------
# HighGCProbs    <- c(0.1, 0.41, 0.39, 0.1) # Set the values of the probabilities, for the GC-rich state
# LowGCProbs    <- c(0.39, 0.1, 0.1, 0.41) # Set the values of the probabilities, for the AT-rich state
HighGCProbs    <- c(0.1, 0.41, 0.39, 0.1)
LowGCProbs    <- c(0.39, 0.1, 0.1, 0.41) 
emissions <- matrix(c(HighGCProbs, LowGCProbs), 2, 4, byrow = TRUE) # Create a 2 x 4 matrix
rownames(emissions) <- states
colnames(emissions) <- nucleotides
 
cat ("\n>>> The HMM emission matrix:\n")
print (emissions)

#-----------------------------------------------------------------------
# This fill the viterbi graph using a matrix v 
# This function is called by the viterbi function (below)
#-----------------------------------------------------------------------
fillViterbiGraph <- function(sequence, transitions, emissions) {
	# Find out how many states are in the HMM
	numstates <- dim(transitions)[1]
	# Make a matrix with as many rows as positions in the sequence, and as many
	# columns as states in the HMM
	v <- matrix(NA, nrow =  dim(transitions)[1], ncol = length(sequence) )
	# Set the values in the first column of matrix v (representing the first position of the sequence) to 0
	v[,1] <- 0
	# Set the value in the first row of matrix v, first column to 1
	v[1,1] <- 1
	# Fill in the matrix v:

	for (i in 2:length(sequence)) { # For each position in the DNA sequence: 
		for (l in 1:numstates) {  # For each of the states of in the HMM:
		   # Find the probabilility, if we are in state l, of choosing the nucleotide at position in the sequence
		   statelprobnucleotidei <- emissions[l,sequence[i]]

		   v[l,i] <-  statelprobnucleotidei * max(v[,(i-1)] * transitions[,l])
		}
	}
	return(v)
}

	
#-----------------------------------------------------------------------
# This carries out the Viterbi algorithm by calling the fillViterbiGraph function
#-----------------------------------------------------------------------
viterbi <- function(sequence, transitions, emissions) {
     theMostProbPath = character()
     # Get the names of the states in the HMM:
     states <- rownames(emissions)

     # Make the Viterbi matrix v:
     v <- fillViterbiGraph(sequence, transitions, emissions)

     # Go through each of the rows of the matrix v (where each row represents
     # a position in the DNA sequence), and find out which column has the
     # maximum value for that row (where each column represents one state of
     # the HMM):
     mostprobablestatepath <- apply(v, 2, function(x) which.max(x))

     # Print out the most probable state path:
     prevnucleotide <- sequence[1]
     prevmostprobablestate <- mostprobablestatepath[1]
     prevmostprobablestatename <- states[prevmostprobablestate]
	 theMostProbPath = c(theMostProbPath, prevmostprobablestatename)
     startpos <- 1
     for (i in 2:length(sequence)) {
        nucleotide <- sequence[i]
        mostprobablestate <- mostprobablestatepath[i]
        mostprobablestatename <- states[mostprobablestate]
        if (mostprobablestatename != prevmostprobablestatename) {
           print(paste("Positions",startpos,"-",(i-1), "Most probable state = ", prevmostprobablestatename))
           startpos <- i
        }
        prevnucleotide <- nucleotide
        prevmostprobablestatename <- mostprobablestatename
		theMostProbPath = c(theMostProbPath, mostprobablestatename)
     }
     print(paste("Positions",startpos,"-",i, "Most probable state = ", prevmostprobablestatename))
	 return (theMostProbPath)
}


#--------------------------------------------------------------
# Main
#--------------------------------------------------------------
stringSequence = "CCCGGGGTTTCCC"
mySequence     = strsplit (stringSequence, split="")[[1]]

# Call to Viterbi to compute the most probable path
theMostProbPath = viterbi (mySequence, transitions, emissions)

cat ("\n>>> The most probable path is:\n")
cat (mySequence)
cat ("\n")
cat (theMostProbPath)
cat ("\n")

