#!/bin/bash
#
#				2020-01-23
#				----------
#

if [ ! -e interpretation.html ]; then
   R --save -q -e "rmarkdown::render('interpretation.Rmd', output_file='interpretation.html')"
fi
