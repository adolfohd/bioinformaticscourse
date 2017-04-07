# Code that shows how to create a HMM and define a function to
# generate a sequence of observationa and states based on the
# model parameters (transitions, emissions, and initials)

nucleotides <- c ("A","C","G","T")
states = c("H","L")
initials = c (0.5, 0.5)

# Create Transitions Matrix
GCHight = c(0.9,0.1)
GCLow = c(0.5, 0.5)
transitions = matrix (c (GCHight, GCLow), 2,2, byrow =TRUE)
rownames (transitions) = states
colnames (transitions) = states

# Create Emissions Matrix
highEmissions = c (0.1, 0.4, 0.4, 0.1)
lowEmissions = c (0.4,0.1,0.1,0.4)
emissions = matrix (c (highEmissions, lowEmissions), 2, 4, byrow=T)
rownames (emissions) = states
colnames (emissions) = nucleotides

#----------------------------------------------------------------------
# Function to generate a sequence based on a HMM parameters
#----------------------------------------------------------------------
generaSeqsHMM <- function (transitions, emissions, initials, N) {
	vectorObservaciones <- character ()
	vectorStates     <- character ()

	# Get the value for the first element (nucleotide and state)
	firstState <- sample (states, 1, rep=TRUE, prob=initials)
	probabilities = emissions [firstState,]
	firstNucleotide = sample (nucleotides, 1, rep=T, prob=probabilities)

	vectorObservaciones [1] = firstNucleotide
	vectorStates [1] = firstState

	# Get the values for the rest of elements (nucleotides and states)
	for (i in 2:N) {
		# Get the current state using the previous state probabilities
		prevState = vectorStates [i-1]             # Previous state (PS)
		prevStateProbs = transitions [prevState,]  # Probs for the PS
		currentState = sample (states, 1, rep=T, prob=prevStateProbs)


		# Get the current nucleotide using the emissions of the current state
		emissionProbs = emissions [currentState,]
		currentNucleotide = sample (nucleotides, 1, rep=T, prob=emissionProbs)

		# Add values to the both vectors
		vectorStates [i] = currentState
		vectorObservaciones [i] = currentNucleotide
	}
	# Packs the results in a list
	return (list (vectorObservaciones, vectorStates))
}

#----------------------------------------------------------------------
# Main
#----------------------------------------------------------------------
observationStates = generaSeqsHMM (transitions, emissions, initials, 10)
observaciones = observationStates [[1]]
states = observationStates [[2]]

cat ("\n",observaciones)
cat ("\n",states,"\n")

