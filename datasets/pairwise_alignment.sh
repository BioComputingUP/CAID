files=(data/fastas/*)
echo $files
total=${#files[@]}
echo $total
i=0
for f1 in "${files[@]}"; do
    i=$(( i + 1))
    j=0
    for f2 in "${files[@]}"; do 
	j=$(( j + 1))
	if (( j > i )); then
	    name1=$(basename $f1 .fasta)
	    name2=$(basename $f2 .fasta)
	    echo $i $j $f1 $f2 $name1 $name2
	    needle $f1 $f2 -gapopen 10.0 -gapextend 0.5 data/pairwise/"$name1"_"$name2".aln
	fi
    done
done
