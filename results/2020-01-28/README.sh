#!/bin/bash
#
#				2020-01-28
#				----------
#

if [ ! -e variation.html ]; then
   R --save -q -e "rmarkdown::render('variation.Rmd', output_file='variation.html')"
fi
