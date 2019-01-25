#bin/bash/

# FILE_1: LIST OF GENES
# FILE_2: LIST WITH FULL ANNOTATION
# FILE_3: LIST WITH full FPKM 

less file_3 |  cut -f 3-20 > fpkm_cut.txt

grep -wFf file_1 fpkm_cut.txt > fpkm_exclusive.txt

grep -wFf file_1 file_2 > annot_exclusive.txt

less annot_exclusive | sed 's/-mRNA-1	/	/'> annot_exclusive_sinmRNA1.txt

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
  }' fpkm_exclusive.txt annot_exclusive_sinmRNA1.txt > full_.txt


