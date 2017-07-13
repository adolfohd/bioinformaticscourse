



# Module 3: Notes


## May 5th, 2017 ([2017-05-05 folder](20170505))

### Non coding RNA

Results are still very dependant on laboratory studies that have not been carried on.

Mother database: [gene bank][genebank]

its FTP holds genomes, proteoms, etc. We are looking if in this DB there is information about non-coding RNA.


Let's go to this FTP address:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/mapview/
```

## Homework 

1. Look in **rna.q** file for rows with the string "non coding"
2. Extract the **GeneID** of these rows
3. Look for these **GeneID**s in the **gene_seq.md** file, and extract the start and end positions of these genes in the genome.
4. Compute the gene's length
5. Visualize an histogram of the non-coding genes lengths.





# Meeting June 13th, 2017


## CSV to NWB ([2017-05-25 folder](20170525))

- Add width and height columns in the nodes declaration section
  - for each node, width and height are equal in order to draw a circle
  - this values are equal to the frequency of appearance of each node in the total of interactions.

- In the definition of interactions area:
  - take the column names from the source CSV file
  - assume the type of data of each column
    - **int** for **source** and **target** columns
    - **float** for numeric columns
    - **string** for non-numeric columns, or a for a column named "label"



[NWB files specification](http://nwb.slis.indiana.edu/software.html)
​		
## Clustering ([2017-06-02 folder](20170602))	

- Include a method that allow us to explain the cluster size
- Based on the code provided by J. Murillo:

```R
wss<- (nrow(scores[1:199, ]-1)*sum(apply(scores[1:199, ], 2, var)))
for (i in 2:15) 
	wss[i] <- sum(kmeans(scores[200:397, ], centers = i)$withinss)
plot(1:15, wss, type ="b", xlab = "Number of cluster", ylab = "withinss groups sum of error squares")
 
library(vegan)
fit <- cascadeKM(scale(scores[1:199,], center =  TRUE, scale = TRUE), 1, 10, iter = 1000)
plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best   
```

## Article Notes

The "rna.q" file contains ... (explain)

In order to filter all non-coding RNA (ncRNA) from the "rna.q" database file, we extract all RNA entries in which description the text "non-coding RNA" is found. Since the same *gene_id* can be found in multiple RNA entries, we only take the one with the longest number of bases.

In order to... (explain). For this, we look the previously selected ncRNA *gene_id*s into the "seq.gene.md" database file. This file contains (explain)

Afterwards, we compute gene lengths by subtracting their final and initial positions in their respective chromosomes ( $chr\_stop - chr\_start$).

Finally, we visualize these data by plotting a histogram of  lengths of all the ncRNA genes found. 

![][ncRNA_Chr]

This figure also divides each one of the histogram's bars into smaller bars with heights proportional to the quantity of genes that belong to each chromosome.

___

Blast results

___

The results of BLAST contain information about the significant interactions among pairs of ncNRA genes. 

In order to (explain), we count the number of interactions each gene has with others, and take averages of the features those interactions have.

Figure X shows an histogram of the number of interactions found for each gene.

![][interHist]

Figure XX shows the histogram of average gene lengths



![][lenHist]

### Clustering

(explain)



In order to choose an adequate number of clusters for this dataset, we will use two criteria: 

- Sums of squared errors


- Calinski-Harabasz criterion

Figure XX shows a plot of the sum of squared errors of all k-mean clusters with respect to their centroids, versus the number of clusters. 

![][sqError]

A rule of thumb to choose a "k" value (number of centroids) is to find an "elbow" in the plot. This is because beyond that point the decrease of sums of squared errors is considerably not as pronounced compared to the left of that "elbow". This value presents itself around a "k" value of 12 to 15.

Figure XX shows a plot for the Calinski-Harabasz criterion. According to the right plot, a "k" value of twelve (12) is adequate.





![][calinski]



According to both criteria, a "k" value of 12 is chosen. Figure XX shows the clustered dataset.

![][kmeans]

___

### Script to paste images

 ```bash
montage -mode concatenate -tile 4x3 cluster\ *.jpg clusters.jpg

 ```

___

Figure XX shows clusters.

![][clusters]





[genebank]:ftp://ftp.ncbi.nlm.nih.gov/genbank
[ncRNA_Chr]:20170505/images/histogram.jpg
[interHist]:20170602/images/histograms/interaction.histogram.jpg
[lenHist]:20170602/images/histograms/length.histogram.jpg
[calinski]:20170602/images/kalinskicriterion.jpg
[sqError]:20170602/images/sumofsquaresVSnumberofgroups.jpg
[kmeans]:20170602/images/clusters/clustering.jpg
[clusters]:20170602/images/clusters/clusters.jpg