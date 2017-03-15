# Module 2

Welcome to my documentation of module 2, lectured by Dr. Luis Garreta.

## Folder structure

- In the folder [material_curso_garreta][materialcurso] we can find the files provided by our lecturer as a sub-module of this repository. I will try to keep this as updated as possible, and maybe contribute to that repo. 




# Classes



## March 8th, 2017





> Homología
>
> Similaridad (computacional)



The similarity of genes does not imply homology (similar structure or function).



> Homework: "Dynamic Programming Matrix exercise"
>
> Also Smith-Waterman
>
> Complete the assignments at the end of the presentation.


## March 15th, 2017



### Global alignment

The next algorithm appears in [this file][align-algo].





![][dynamic]



The next code appears in [this file][global-code]

```python


# Program to compute the global and local pairwise alignment
# It implements the Smith-Waterman Algorithm for global
# alignments 

# Main function where the program starts
def main ():
    seq1 = "ACGT"
    seq2 = "ACCGCT"
    scoringPenalties = {"match" : 2, "mismatch" : -1, "gap" : -2}
    #scoringPenalties = {"match" : 10, "mismatch" : -5, "gap" : -7}
    print "\n>>> Global Alignment of %s and %s" % (seq1, seq2)
    matrix = global_align (seq1,seq2, scoringPenalties)
    printMatrix (seq1, seq2, matrix)

# Construct the alignment matrix filled with zeros 
def make_matrix(seq1, seq2):
    nColumns = len (seq1) + 1
    nFilas = len (seq2) + 1
    return [[0]*nColumns for i in xrange(nFilas)]

# Construct a global alignment
def global_align(seq1, seq2, scoringPenalties):
    # create a zero-filled matrix
    A = make_matrix (seq1, seq2)
    A[0][0] = 0

    # Initialize the first row with gap penalties 
    for i in range (1, len(seq2)+1):
      A[i][0] += scoringPenalties["gap"] + A[i-1][0]
    # Initialize the first column with gap penalties  
    for j in range (1, len(seq1)+1):
      A[0][j] += scoringPenalties ["gap"] + A[0][j-1]
    
    # fill in A in the right order
    for i in xrange(1, len(seq2)+1):
	for j in xrange(1, len(seq1)+1):
	    typeSubstitution = "match" if seq2[i-1] == seq1[j-1] else "mismatch"

	    diagonalScore   = A[i-1][j-1] + scoringPenalties [typeSubstitution]
	    horizontalScore = A[i][j-1] + scoringPenalties ["gap"]
	    verticalScore   = A[i-1][j] + scoringPenalties ["gap"]

	    A[i][j] = max (horizontalScore, verticalScore, diagonalScore)

    return A

# Print the matrix with the partial scores
def printMatrix(seq1, seq2, matrix):
	seq1 = "-" + seq1
	seq2 = "-" + seq2

	print "-",
	for j in seq1:
		print j.rjust(3),
	print ("")
	
	for i in range (0, len (matrix)):
		print seq2 [i],
		for j in range (0, len (matrix[i])):
			print str (matrix [i][j]).rjust(3),

		print ""


#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
if __name__ == "__main__":
	main ()


```



### Local alignment

- No initial gaps are penalized, so the first column and first row are filled with zeroes.
- The local alignment scoring function  always returns values greater or equal than zero.



### Motifs

Pattern that is suposed to have a biological function.



#### Regulatory motifs



Activators or regressors: Turns a gene on or off. 

The regulatory region of a gene typically starts 100-1000 bp upstream.

> Notation 3' - 5' , the direction of strand



The double helix is made of one positive strand and a negative strand.



Transcription factor identify the regulatory motif so transcription starts from there.

> notations for motifs (this is explained in the course's slides)
>
> - Logo notation
> - Consensus notation



### Construction of a PWM



The simplest way of modeling the backgroung is to assign an uniform distribution:

| Nucleotide | Probability |
| ---------- | ----------- |
| A          | 0.25        |
| C          | 0.25        |
| G          | 0.25        |
| T          | 0.25        |



With this, we build a weight vector for the background:

$ b = \left[  \array{ 0.25 \\ 0.25 \\ 0.25 \\ 0.25 }  \right]$

### Slide PWM



every column id divided by the respective background weight:



$M'_{i,j} = M_{i,j} \div b_i$



[materialcurso]: material_curso_garreta/

[dynamic]: images/dynamic.png
[align-algo]: material_curso_garreta/00-alignments/00-bioinfo-alignments-algorithmDP.pdf
[global-code]: material_curso_garreta/00-alignments/code/global_alignment.py

