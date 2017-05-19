import sys
import os
import re

rnaq = open("../data/humangenome/rna.q")
a = rnaq.readline()
noncoding = []
while a != "":
	b = re.search("non-coding", a) 
	if b is not None:
		a = re.search("(?<=GeneID:)\\d+\\t", a).group(0)
		noncoding.append(a)		
	a = rnaq.readline()
rnaq.close()
# Filter duplicates
noncoding = list(set(noncoding))

seq_gene = open("../data/humangenome/seq_gene.md")
a = seq_gene.readline()
while a != "":
	for geneid in noncoding:
		b = re.search("GeneID:"+geneid+"$", a)
		if b is not None:
			print(a)
	a = seq_gene.readline()
