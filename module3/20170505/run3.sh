
echo "Searching Gene IDs"
args=$(grep "non-coding" ../data/humangenome/rna.q | grep -Po "(?<=GeneID:)\\d+")
echo "Gene IDs found"
echo "${args[@]:0:100}, etc"

echo "Searching in MD file"
while read -r line
do
	echo "$line"
    for i in ${args[@]}
    do
		query="GeneID:$geneid\$"
                grep "$query" <<< "$line"
#		grep "GeneID:$i\$" <<< "$line"
    done
done <"../data/humangenome/seq_gene.md"


        # echo "$geneid"

#        awk '$12=="GENE"' ../data/humangenome/seq_gene.md | while read -r seq; do
#                query="GeneID:$geneid\$"
#                grep "$query" <<< "$seq"
#        done
#done


