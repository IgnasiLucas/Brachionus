#!/bin/bash
#
#				2019-04-03
#				==========
#
# Here I will run edgeR to identify differentially expressed genes. The samples
# consist of 3 replicates of the 4 combinations of two two-level treatments:
# the selective regime (predictable or unpredictable) and the hatching time
# (either after a forced 28 days diapause, or right away). The distribution of
# samples among treatments was the following:
#
#  +------------------------------------+------------------------------------+
#  |    Predictable selective regime    |  Unpredictable selective regime    |
#  +------------------------------------+------------------------------------+
#  | Forced diapause | Without diapause | Forced diapause | Without diapause |
#  +-----------------+------------------+-----------------+------------------+
#  |     3C_S11      |      3A_S9       |     1C_S1       |      1A_S8       |
#  |     5C_S4       |      5A_S2       |     2C_S5       |      2A_S7       |
#  |     6C_S10      |      6A_S3       |     4C_S12      |      4A_S6       |
#  +-----------------+------------------+-----------------+------------------+
#
# Note that samples are paired between diapause treatments. After the selection
# experiment, "diapausing eggs obtained for each laboratory population were divided
# into two groups to be subjected to two different diapause conditions: (1) non-forced
# diapause (NFD) and (2) forced diapause (FD)." Thus, when testing for the effect
# of diapause (or hatching) condition, it would not be efficient to compare the
# average level of expression among samples forced to hatch with the average level
# among samples forced to diapause.
#
# One way to model gene expression in these samples is
# described in section 3.5 of the edgeR manual. Since the effect of the population
# only matters within a selective regime, we can nest it, and re-use the labels
# for the populations. Below, I order the populations in a more natural way:
#
#  +------------------+------------+----------+-----------+
#  | Selective regime | Population | Hatching | Sample Id |
#  +------------------+------------+----------+-----------+
#  |      Random      |      1     | hatching |   1A_S8   |
#  |      Random      |      1     | diapause |   1C_S1   |
#  |      Random      |      2     | hatching |   2A_S7   |
#  |      Random      |      2     | diapause |   2C_S5   |
#  |      Random      |      3     | hatching |   4A_S6   |
#  |      Random      |      3     | diapause |   4C_S12  |
#  |     Periodic     |      1     | hatching |   3A_S9   |
#  |     Periodic     |      1     | diapause |   3C_S11  |
#  |     Periodic     |      2     | hatching |   5A_S2   |
#  |     Periodic     |      2     | diapause |   5C_S4   |
#  |     Periodic     |      3     | hatching |   6A_S3   |
#  |     Periodic     |      3     | diapause |   6C_S10  |
#  +------------------+------------+----------+-----------+
#
# A model such as this:
#
#       ~regime + regime:population + hatching
#
# would estimate the orthogonal effects of regime and hatching, but not their
# interaction. Maybe I should first test for the significance of the interaction.
# Well, because the same model must be fitted to every gene, the way to decide
# what model is better is to look at the number of genes with a significant interaction
# term, I guess.

GENE_COUNTS=../2019-03-29/genes.PostCount.txt
ISOFORM_COUNTS=../2019-03-29/isoforms.PostCount.txt


if [ ! -d genes ]; then mkdir genes; fi
if [ ! -e genes/report.html ]; then
   R --save -e "rmarkdown::render('RunEdgeR.Rmd', output_file='genes/report.html')" --args $GENE_COUNTS genes
   mv .RData genes/
fi

if [ ! -d isoforms ]; then mkdir isoforms; fi
if [ ! -e isoforms/report.html ]; then
   #Rscript --no-save RunEdgeR.R $ISOFORM_COUNTS isoforms
   R --save -e "rmarkdown::render('RunEdgeR.Rmd', output_file='isoforms/report.html')" --args $ISOFORM_COUNTS isoforms
   mv .RData isoforms/
fi
