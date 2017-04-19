#### Module 2

Welcome to my documentation of module 2, lectured by Dr. Luis Garreta.

## Folder structure

- In the folder [material_curso_garreta][materialcurso] we can find the files provided by our lecturer as a sub-module of this repository. I will try to keep this as updated as possible, and maybe contribute to that repo. 

- To update (clone) that submodule for the first time, go into the folder and...

  ```
  $ git submodule init
  $ git submodule update
  ```

  For additional updates:

  ```
  $ git fetch
  $ git merge origin/master
  ```

  ​




# Classes



## March 8th, 2017





> HomologÃ­a
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

> Review distances between strings



> What's this?
>
> - 3 - amero
> - 3 - mers
> - k - mers



___


> make slides from:
- ["*Significance of Alignments*"][lesk_fragment]. [Introduction to Bioinformatics][lesk]¨ (Arthur M. Lesk) p.267 (Biblioteca PUJ).
   - ["*Databases Searches*"][krane_fragment]. [Fundamental Concepts of Bioinformatics][krane] (Krane & Raymer) p.48. (PDF)


> **R/.** Slides are in [this file][01slides] 





## March 22nd, 2017


> How to predict when a dice is charged?


### Hidden Markov Models (HMM)

#### Gene finding and prediction

- Genes have 
- Start
- Finish
- Parts
  - exons
  - codond (flags)
  - start codons
  - middle codons
  - donor site
- Poly-A site
  - very distinctive


#### Protein profile HMMs

 > complete this

### HMM Architecture

#### Finite state machine
- Deterministic model
- Markov assumption:
- Next state depends only on all previous states
- Order $n$ model
- Next state depends only on $n$ previous states
- Simplest: Order 1


#### State transition matrix


$a_{ij} \rightarrow$ transition probability from state $i$ to state $j$

We always have that $ \sum_{j=1}^K {a_{ij}}$

e.g.:

$A= \left\{ \array{   & Sun  & Cloud & Rain \\ Sun & 0.5 & 0.25 & 0.25 \\ Cloud & 0.375 & 0.125 & 0.375 \\ Rain & 0.125 & 0.625 & 0.375 } \right\}$


#### Initial state

$\Pi = \left\{ \array{ Sun & Cloud & Rain \\ 1 & 0 & 0 }  \right\} $ 

#### Two sets of states

- Observable states (moisture in seaweed)
- Hidden states $\rightarrow$ what I want to predict (weather)


#### Emission matrix




$\array{ & \bf{Seaweed} \\  \bf{Weather} &  Dry & Dryish & Damp & Soggy \\ Sun & \\ Cloud \\ Rain}$



### The problems of HMMs



- The first problem is to determine if a sequence of dice results belongs to our model 
- **Evaluation question of HMMs**
- What is the state of the dice (fair or loaded) for any given dice result?
- **Decoding question of HMMS**
- How the state **transition** and the **emission** matrices change as more results are thrown by the croupier?
- **Learning question of HMMs**

### Hidden states and observable sequence

$\pi_t$: hidden state
$x_t$: emitted symbol during $\pi_t$

Path: $\pi_1, \pi_2,...,\pi_N$

### Memory-less

$P(\pi_{t+1}= k | \mbox{ 'all that has happened so far'} )$
$P(\pi_{t+1}= k | x_1, x_2, ..., x_t )$
$P(\pi_{t+1}= k | x_t $


> **CpG islands**
- Higer probability of finding a gene in subsequences containing a high ratio of C and G nucleotides.


### Problem 2: Decoding

- Find the best parse of a sequence



> **Homework**:
> Find ADN sequences with previously identified CpG islands.



## March 29th, 2017



Eucariots have start codons = **ATG**

introns are not transcripted into RNA i.e. introns in DNA will not make part of ARN, much less of protein.

**Gramatics of DNA**

5'UTR-(EXON-INTRON)*-EXON-3'UTR

EXONs have different structure than INTRON. This allows to search for them or predict their location based on this structure.



> Polymerase performs the search of specific parts in DNA in order to... ???



### Architecture:

1. Next hidden state matrix: $A=$ probabilities given a current state
2. Initial states: $\Pi=$ 
3. Emissions matrix: $B=$ probabilities of a hidden state given an observable state



### Markov's 3 problems

1. **Evaluation**: how likely is the sequence explained by my model?

2. **Decoding**: which portion of the sequence belongs to one state or other?

3. **Learning**: e_{\pi_i} 

   ​

$P(x, \pi) = \prod_{i=1}^N a_{\pi_{i-1},\pi_{i} }e_{\pi_i} (x_i)$ 



### Decoding





## April 7th, 2017

> Pathways in protein structure



Pathways are described by a relatively small amount of steps: 100-200 steps

On the other hand, molecular dynamics usually present millions of steps.

In pathways, multiple metrics are evaluated:

- Energy (difficult)
- External area of absorption
- Quantity of natural contacts.
- Hydrogen bonds
- etc., typically up to **15 features**

Problems?

- Large amount of features
  - How to reduce them? Maybe one variable is highly explained by other(s)
  - Some non-measurable metrics can be estimated using others.
- Large amount of molecular structures 
  - How can we propose a smaller representative group? 
  - How to get them?
    - The Anton computer has produced large amount on which we can work.



### Filtering and clustering





> NCDI: Non-redundant database



There are huge databases that, are very useful by themselves, but as users it is not feasible to download all data, due to bandwidth and storage constraints.

>  Anton, David Shaw: Terabytes of information are generated in that database just for a 20 to 100 amino-acids protein



> Free energy of a molecular structure
>
> Greedy: CD-Hit
>
> MUSCLE, UCLUST, USEARCH



> Single linkage
>
> Complete-linkage: Very sensitive to outliers because the distance for decisions of clustering two groups is the one of the farthest away elements. 
>
> Average-linkage (=UPGMA)



> Try in R:   plot(hclust(dist()), method="single" )



### K-means



> Lloyd's algorithm



> Squared-error-distorsion (SED)



> Bayesian clustering



### Representative datasets







> Hobbon's 2 algorithms
>
> 



#### Hobbon1

It is greedy in the sense that the very first element different enough from the others is considered to be the centroid of its cluster.

#### CD-Hit (Li, et al.)

It's a greedy method.

- At any given point, the optimal is assumed to be the local solution.

CDHIT is an open source program.

MUSCLE replaces CLUSTALW (not open source - very closed in fact)

> Application of each one?



> Comparisson of k-mers? $\rightarrow$ Look for the generalities of this algorithm



##### CDHIT basic algorithm

- [How CD-HIT Works?][CDHIT]



$X={x_1,...x_n}$, $x_i$ is a short sequence of nucleotides (7 to 15 base pairs)



The k-mer distance is measured using:

> Bush
>
> Smith



> Homework: Implement CD-HIT
>
> - The distance can be a hamming


## April 19th, 2017





[materialcurso]: material_curso_garreta/

[dynamic]: images/dynamic.png
[align-algo]: material_curso_garreta/00-alignments/00-bioinfo-alignments-algorithmDP.pdf
[global-code]: material_curso_garreta/00-alignments/code/global_alignment.py

[lesk]: biblio/lesk.pdf
[lesk_fragment]: 01slides/lesk_fragment.pdf
[krane]: biblio/krane.pdf
[krane_fragment]: 01slides/krane_fragment.pdf

[01slides]: 01slides/01slides.pdf
[CDHIT]: http://blog.nextgenetics.net/?e=26

