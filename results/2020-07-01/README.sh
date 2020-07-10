#!/bin/bash
#
#				2020-07-01
#				==========
#
TRANSCRIPTS='../2019-07-10/transcripts.fa'

# transcript2gene is a space separated text file with two columns: transcript
# and gene names.
if [ ! -e transcript2gene.txt ]; then
   grep "^>" $TRANSCRIPTS | sed -E 's/^>|gene_id://g' > transcript2gene.txt
fi

if [ ! -e interpretation.html ]; then
   R --save -q -e "rmarkdown::render('interpretation.Rmd', output_file='interpretation.html')"
fi
