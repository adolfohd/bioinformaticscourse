

# Classes


## May 5th, 2017

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




#####

[genebank]:ftp://ftp.ncbi.nlm.nih.gov/genbank