#!/bin/bash

ANNOTATION=/home/eva/Escritorio/diff_out/2017-12-07/B_plicat_annotation_170209.annot
WORKING=/home/eva/Escritorio/diff_out/2017-12-22

#remove first two columns from gne_diff file with all RPKM information

for i in ./list_genes_FC2/*_all.txt ; do 
	echo $i
	name=$(basename "$i" | cut -d "_" -f 1)
	echo $name
#	cat $i | cut -f 3-20 > $name\_fpkm_cut.txt
	cat $i | cut -f 3-20 > ./$name\_fpkm_cut.txt
done

#Add an annotation information to the file with only list of genes exclusively and shared

for file in $WORKING/list_genes_venn/* ; do
	echo $file
	name2=$(basename "$file" | cut -d "_" -f 1,2)
	echo $name2 
	grep -wFf $file ./B_plicat_annot_sinmRNA1.txt >./$name2\_annot.txt
#Change name to the gene in DEGs_annotation
	#cat ./$name2\_annot.txt | sed 's/-mRNA-1	/	/' > ./$name2\_annot_sinmRNA1.txt
done

echo "DONE"
