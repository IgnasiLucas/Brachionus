#!/bin/bash
#
#				2020-01-08
#				==========
#
# After a conversation with Edu, we decided to give another try to fitting a split-plot,
# mixed model. There is the *dream* method, which I have not tried yet.
#
# Notes on reproducibility
# ------------------------
#
# Lastly, I am using RStudio-server daily. It is very confortable, but it makes it more
# difficult to reproduce the actual environment. I cannot rely on conda environments to
# offer you the packages and versions I am using, because from within the active rstudio-
# server session I cannot switch conda environments. Suggestions online involve changing
# rstudio-server's default environment, which is not what I want.

if [ ! -d genes ]; then mkdir genes; fi
if [ ! -e genes/mixed.html ]; then
   R --no-save -q -e "rmarkdown::render('mixed_genes.Rmd', output_file='genes/mixed.html')"
fi

if [ ! -d isoforms ]; then mkdir isoforms; fi
if [ ! -e isoforms/mixed.html ]; then
   R --no-save -q -e "rmarkdown::render('mixed_isoforms.Rmd', output_file='isoforms/mixed.html')"
fi
