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

for folder in genes isoforms; do
   if [ ! -d $folder ]; then mkdir $folder; fi
   for report in Preliminar Hatching Interactions Regime; do
      if [ ! -e $folder/$report.html ]; then
         R --save -e $(printf "rmarkdown::render('%s.Rmd',output_file='%s.html',output_dir='%s')" $report $report $folder) --args ../2019-03-29/$folder.PostCount.txt $folder
      fi
   done
   R --save -e "rm(list=ls())"
done

# CONCLUSIONS
# -----------
#
# Very few genes respond consistently to hatching condition, which is no unexpected. If I recall
# correctly, hatching condition (either immedate hatching or induced diapause) was expected to
# modulate the set of genes the expression of which could be affected by selective regime.
#
# In this respect, models 5 and 6 fail to identify any gene or isoform with a significant interaction
# term between selective regime and hatching condition.
#
# In general, genes affected by regime seem to have their expression changed in the same direction
# whether when immediately hatching or not. And significant genes under both conditions overlap extensively,
# according to model 3.
#
# However, model 3 does not account for a possible batch effect of populations. Every pair of immediately
# hatching or diapausing populations was derived from a common population at the end of the selection
# experiment. Thus, there were 6 such pairs, derived from 6 selection experiments. Selection lasted
# 8 cycles of either regular or unpredictable dissecation-hidration periods. It is possible that genetic
# drift, together with selection, differentiated those populations to the point of creating a batch
# effect.
#
# Model 4 is considered better than any other model. It includes the population effect, in addition to
# the main effects of hatching condition and selective regime. According to model 4, no more than 65
# genes (56 isoforms) are up- or down-regulated by the selective regime.
