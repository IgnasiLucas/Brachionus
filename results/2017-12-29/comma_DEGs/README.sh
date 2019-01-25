#!bin/bash

ANNOTATION=/home/eva/Escritorio/diff_out/2017-12-20/comma_DEGs/B_plicat_annot_sinmRNA1.txt
WORKING=/home/eva/Escritorio/diff_out/2017-12-22/list_genes_venn

for file in $WORKING/*.txt ; do
        echo $file
        name2=$(basename "$file" | cut -d "_" -f 1,2)
        echo $name2 
        cat $file | grep "," | sed 's/$/\n-----/g'| sed 's/,/\n/g' > ./parts/$name2\_comma.txt
done
# cat B_plicat_annotation_170209.annot | sed 's/-mRNA-1     /       /' > B_plicat_annot_sinmRNA1.txt
for file2 in ./parts/*_comma.txt ; do
        echo $file2
        name3=$(basename "$file2" | cut -d "_" -f 1,2)
        echo $name3 
	grep -wFf $file2 $ANNOTATION > ./parts/$name3\_function.txt

	awk 'FNR == NR { lineno[$1] = NR; next}     {print lineno[$1], $0;}' ./parts/$name3\_comma.txt ./parts/$name3\_function.txt  | sort -k 1,1n | grep -v "\---NA\---" > ./$name3\_list_noNA.txt

done
echo "DONE"


