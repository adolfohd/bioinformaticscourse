geneids=$(grep "non-coding" ../data/humangenome/rna.q)

pat=$(echo ${geneids[@]}|tr " " "|")

# grep -Eow "$pat" ../data/humangenome/seq_gene.md

awk '$12=="GENE"' ../data/humangenome/seq_gene.md | while read -r seq; do
        	"
	grep -Eow "$pat" <<< "$seq"

done
