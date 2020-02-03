#!/bin/bash
#
#				2020-01-14
#				----------
#
# I separate genes from isoforms. Then, the enrichment analysis can be run on genes
# that are differentially expressed between selective regimes, and also on those
# differentially expressed between hatching condition. For the moment, I do not
# aim at functions associated with genes regulated by the interaction between
# selective regime and hatching condition.
#
# Plus, each analysis can be performed using two different orderings of the genes.
# The effect of either selective regime or hatching condition has been measured
# on 2020-01-08 in two different ways: with a partition of variance, and with a
# test for differential expression. Both options will be reported in the same
# document.

for tag in genes isoforms; do
   if [ ! -d $tag ]; then mkdir $tag; fi
   for var in regime hatching; do
      if [ ! -d $tag/$var ]; then mkdir $tag/$var; fi
      if [ ! -e $tag/$var/enrichment.html ]; then
         # I used to pass parameters to R when rendering an .Rmd file, in order to re-use it
         # with a different data set. But then, I did no know how to render the .Rmd file from
         # Rstudio. I resorted to using more than one, almost identical .Rmd files (see 2020-01-08).
         # But I just discovered I can 'knit with parameters' from Rstudio!
         # Note that for the parameters to be available in render(), to set to output file, I need
         # to call render within a function that uses those parameters as arguments.
         R --no-save -q -e "render_report <- function(tag, var){
                               rmarkdown::render('enrichment.Rmd',
                                  params = list(TAG = tag, VAR = var),
                                  output_file = paste(tag, var, 'enrichment.html', sep='/'))
                            }" \
                        -e "render_report('$tag', '$var')"
      fi
   done
done
