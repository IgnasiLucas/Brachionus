#!/bin/bash

#File 1: list genes exclusively and shared (esta en la carpeta 2017-12-22)
#File 4: list full RPKM (recien hechos quitando las 2 primeras columnas)

WORKING=/home/eva/Escritorio/diff_out/2017-12-22/list_genes_venn

grep -wFf $WORKING/I1vsI2_genes_exclusively.txt ./I1vsI2_fpkm_cut.txt > ./genes_fpkm/I1vsI2_fpkm_exclusively.txt
grep -wFf $WORKING/P1vsP2_genes_exclusively.txt ./P1vsP2_fpkm_cut.txt > ./genes_fpkm/P1vsP2_fpkm_exclusively.txt
grep -wFf $WORKING/I1vsP1_genes_exclusively.txt ./I1vsP1_fpkm_cut.txt > ./genes_fpkm/I1vsP1_fpkm_exclusively.txt
grep -wFf $WORKING/I2vsP2_genes_exclusively.txt ./I2vsP2_fpkm_cut.txt > ./genes_fpkm/I2vsP2_fpkm_exclusively.txt
grep -wFf $WORKING/I1vsP1_I2vsP2_genes_shared.txt ./I1vsP1_fpkm_cut.txt > ./genes_fpkm/I1vsP1_I2vsP2_fpkm_exclusively.txt
grep -wFf $WORKING/I1vsI2_P1vsP2_genes_shared.txt ./I1vsI2_fpkm_cut.txt > ./genes_fpkm/I1vsI2_P1vsP2_fpkm_exclusively.txt


awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/I1vsI2_fpkm_exclusively.txt ./I1vsI2_genes_annot.txt > ./full_final/I1vsI2_full.txt

echo "DONE_1"

awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/P1vsP2_fpkm_exclusively.txt ./P1vsP2_genes_annot.txt > ./full_final/P1vsP2_full.txt

echo "DONE_2"

awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/I1vsP1_fpkm_exclusively.txt ./I1vsP1_genes_annot.txt > ./full_final/I1vsP1_full.txt

echo "DONE_3"

awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/I2vsP2_fpkm_exclusively.txt ./I2vsP2_genes_annot.txt > ./full_final/I2vsP2_full.txt

echo "DONE_4"

awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/I1vsP1_I2vsP2_fpkm_exclusively.txt ./I1vsP1_I2vsP2_annot.txt > ./full_final/I1vsP1_I2vsP2_full.txt

echo "DONE_5"

awk 'BEGIN {
  FS = OFS = "\t"
  }
NR == FNR {
  # while reading the 1st file
  # store its records in the array f
  f[$1] = $0
  next
  }
$1 in f {
  # when match is found
  # print all values
  print f[$1], $0
  }' ./genes_fpkm/I1vsI2_P1vsP2_fpkm_exclusively.txt ./I1vsI2_P1vsP2_annot.txt > ./full_final/I1vsI2_P1vsP2_full.txt

echo "DONE_6"
