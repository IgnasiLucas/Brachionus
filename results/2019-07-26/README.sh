#!/bin/bash
#
#				2019-07-26
#				==========
#
# To perform the functional enrichment analysis, I want to use the R package topGO.
# TopGO is installed through conda in the environment 'topGO'. I record here the
# spec file and the yaml file of my topGO environment to help reproducibility.
#
# There are four enrichment analysis: genes between selective regimes, genes between
# hatching conditions, isoforms (or transcripts) between selective regimes and
# isoforms between hatching conditions.

TAGTYPE=( genes transcripts )
ANNOTATION=( ../2019-07-19/gene2GO.txt ../2019-07-19/trans2GO.txt )
LIST_DIRS=( ../2019-04-03/genes ../2019-04-03/isoforms )
CONDITION=( regime hatching )

for i in 0 1; do
   if [ ! -d ${TAGTYPE[$i]} ]; then mkdir ${TAGTYPE[$i]}; fi
   if [ ! -e ${TAGTYPE[$i]}/annotation.txt ]; then
      grep "GO:" ${ANNOTATION[$i]} > ${TAGTYPE[$i]}/annotation.txt
   fi
   for condition in ${CONDITION[@]}; do
      if [ ! -e ${TAGTYPE[$i]}/$condition.html ]; then
         # R uses white space as separator when parsing the expression. I can use "-e" more than once to
         # execute more than one expression. This is clearer than embeding the R expression in a printf expression
         # (see ../2019-04-03/README.sh) to pass arguments to the R expression itself. The 'Enrichment.Rmd' script
         # requires two additional arguments, the list of genes and the annotation file.
         R --save -q -e "OutputFile='${TAGTYPE[$i]}/$condition.html'" \
                  -e "rmarkdown::render('Enrichment.Rmd',output_file=OutputFile)" \
                  --args ${LIST_DIRS[$i]}/$condition.txt ${TAGTYPE[$i]}/annotation.txt
      fi
   done
done

# You can see a rendered version of the reports using the following links:
#
# Genes differentially expressed between selective regimes:
#    https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/genes/regime.html
#
# Genes differentially expressed between hatching conditions:
#    https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/genes/hatching.html
#
# Transcripts differentially expressed between selective regimes:
#    https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/transcripts/regime.html
#
# Transcripts differentially expressed between hatching conditions:
#    https://htmlpreview.github.io/?https://github.com/IgnasiLucas/Brachionus/blob/master/results/2019-07-26/transcripts/hatching.html

