grep "non-coding" ../data/humangenome/rna.q | while read -r line ; do
    
        geneid=$(grep -Po "(?<=GeneID:)\\d+" <<< "$line")
        # echo "$geneid"
        awk '$12=="GENE"' ../data/humangenome/seq_gene.md | while read -r seq; do
                query="GeneID:$geneid\$"
                grep "$query" <<< "$seq"
        done
done

