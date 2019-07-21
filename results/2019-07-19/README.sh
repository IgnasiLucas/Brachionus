#!/bin/bash
#
#				2019-10-19
#				==========
#

INTERPROSCAN=~/bin/my_interproscan/interproscan-5.36-75.0/interproscan.sh
PROTEINS=../2019-07-10/transcripts.fa.transdecoder.pep

if [ ! -e proteins.fa ]; then
   sed 's/\*//g' $PROTEINS > proteins.fa
fi

if [ ! -e functional_annotation.tsv ]; then
   $INTERPROSCAN --disable-precalc --goterms --input proteins.fa --output-file-base functiona_annotation --disable-residue-annot --cpu 50 --highmem
fi

