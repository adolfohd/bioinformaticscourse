

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




#####

[genebank]:ftp://ftp.ncbi.nlm.nih.gov/genbank
