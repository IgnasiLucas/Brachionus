#!/bin/bash
#
#				2019-10-19
#				==========
#
# Here the goal is to use interproscan to functionally annotate the transcripts previously
# analysed. Interproscan outputs much more information that I can use now. I will limit the
# functional analysis to the Gene Ontology terms.

INTERPROSCAN=~/bin/my_interproscan/interproscan-5.36-75.0/interproscan.sh
PROTEINS=../2019-07-10/transcripts.fa.transdecoder.pep
TRANSCRIPTS=../2019-07-10/transcripts.fa

if [ ! -e proteins.fa ]; then
   sed 's/\*//g' $PROTEINS > proteins.fa
fi

if [ ! -e functional_annotation.tsv ]; then
   $INTERPROSCAN --disable-precalc --goterms --input proteins.fa --output-file-base functiona_annotation --disable-residue-annot --cpu 50 --highmem
fi

if [ ! -e gene2GO.txt ] || [ ! -e trans2GO.txt ]; then
   # The $TRANSCRIPTS file, from previous folder, contains the information on what gene each
   # transcript belongs to. The GO terms of transcripts are in the 'functional_annotation.tsv' file.
   # Below, I give gawk both pieces of information, one as standard input, and one as a file to read.
   # Then I use an array of arrays to keep track of the GO terms associated to each transcript and gene.
   # To make sure that each GO term is registered only once in a gene's annotation, I record them as
   # keys of subarrays.
   #
   # Note also the use of a named pipe and the 'tee' command to split gawk's output in two different
   # files. I could also have asked gawk to print to two different files.
      # Apparently, I need to create and delete an element of the subarrays, to make gawk aware that
      # these are subarrays. Otherwise, it stops and complains that I'm using a scalar in array context.
   mkfifo pipe
   cat pipe | grep "^XLOC" | sort -nk 1,1 > gene2GO.txt &
   grep "^>" $TRANSCRIPTS | sed -E 's/^>|gene_id://g' | gawk '(FILENAME != "functional_annotation.tsv") {
      GENE[$1] = $2
      ANNOTATION[$1]["GO:1"] = 1
      ANNOTATION[$2]["GO:1"] = 1
      delete ANNOTATION[$1]["GO:1"]
      delete ANNOTATION[$2]["GO:1"]
   }((FILENAME == "functional_annotation.tsv") && (/GO:/)){
      TRANS = $1
      gsub(/\.p[0-9]+/, "", TRANS)
      split($NF, TERMS, /\|/)
      for (i in TERMS) {
         if (!(TERMS[i] in ANNOTATION[TRANS])) ANNOTATION[TRANS][TERMS[i]] = 1
         if (!(TERMS[i] in ANNOTATION[GENE[TRANS]])) ANNOTATION[GENE[TRANS]][TERMS[i]] = 1
      }
   }END{
      for (tag in ANNOTATION) {
         LIST = ""
         SEP  = ""
         for (term in ANNOTATION[tag]) {
            LIST = LIST SEP term
            SEP  = "|"
         }
         print tag "\t" LIST
      }
   }' - functional_annotation.tsv | tee pipe | grep "^TCONS" | sort -nk 1,1 > trans2GO.txt
   # Note the "-" after the awk code. Otherwise it does not read standard input, but only the file.
   rm pipe
fi

if [ ! -e summary.txt ]; then
   gawk 'BEGIN{
      WHATIS["TCON"] = "transcript"
      WHATIS["XLOC"] = "gene"
   }{
      NUM[WHATIS[substr($1,1,4)]]++
      split($2,A,/\|/)
      FREQ[WHATIS[substr($1,1,4)]][length(A)] += 1
   }END{
      print  "-----------------------------------------------"
      print  " Tag           Annotated  Terms/Tag      Total"
      print  "-----------------------------------------------"
      for (tag in NUM) {
         delete FREQ[tag][0]
         SUM = 0
         NUM_ANNOTATED = 0
         for (size in FREQ[tag]) {
            SUM += FREQ[tag][size] * size
            NUM_ANNOTATED += FREQ[tag][size]
         }
         printf " %-11s % 11d % 10.2f % 10d\n", tag, NUM_ANNOTATED, SUM/NUM_ANNOTATED, NUM[tag]
      }
      print  "-----------------------------------------------"
   }' *2GO.txt > summary.txt
fi

# CONCLUSIONS
# -----------
#
# There are 53789 genes and 77728 transcripts. But only 9691 (18%) genes got annotated with at
# least one GO term using InterProScan. Among the transcripts, 19735 (25.4%) received GO annotation.
#
# -----------------------------------------------
#  Tag           Annotated  Terms/Tag      Total
# -----------------------------------------------
#  gene               9691       2.57      53789
#  transcript        19735       2.55      77728
# -----------------------------------------------
#
# It seems insufficient for a functional enrichment analysis.
