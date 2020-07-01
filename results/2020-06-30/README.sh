#!/bin/bash
#
#				2020-06-30
#				==========
#

if [ ! -e enrichment.html ]; then
   R -q --no-save -e "rmarkdown::render('enrichment.Rmd', output_file='enrichment.html')"
fi
