#!/bin/bash
#				2020-06-29
#				==========
#

if [ ! -e figures.html ]; then
   R -q --no-save -e "rmarkdown::render('figures.Rmd', output_file = 'figures.html')"
fi
