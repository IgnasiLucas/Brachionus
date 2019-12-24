#!/bin/bash
#
#				2019-12-23
#				==========
#
# Let's summarize the results up to now. Six original populations were
# evolved in two different regimes: three under a regular environment,
# and three under a random one. After the evolution experiment, every
# evolved population was split in two: one was treated for immediate
# hatching of resting eggs, and the other was forced to undergo diapause.
# Gene expression was measured in every population. Genes with too low
# expression level in more than 8 libraries were dropped. The minimum
# expression level required was such that would produce at least 5 counts
# in the smallest library.
#
# The filtering step, being less aggressive than suggested by the edgeR
# manual, reduced the number of genes from 53789 to 18497. In the case of
# isoforms, the filtering left 32925 out of 77728. This step may explain
# the main differences between the EdgeR and the Cufflinks analyses. In
# any case, it is clear that the filtering is necessary and cannot be softer.
#
# Among the different statistical models, the ones that account for the
# population effect produce a better fit, suggesting that not all populations
# responded equally to the treatments.
#
# The hatching condition by itself affected the expression of very few genes.
# Only two genes (XLOC_052427 and XLOC_045087), and one transcript of one of
# them (TCONS_00075780) are significantly (down)regulated in forced diapause
# with respect to the immediate hatching condition. One of those genes, XLOC_052427,
# codes for a protein with homologues in a few other species. The other doesn't.
#
# What genes are regulated by selective regime is less clear, because the number
# varies wildly among models. It seems, however, that the large numbers of
# regulated genes suggested by models 2, 3 and 5 are due to their failure to
# account for the dependence between samples descended from the same original
# population. The experiment is better modeled as a split-plot design, where
# the selective-regime levels are randomly assigned to 6 original populations
# (whole-plot factor), while diapause treatments are applied to sub-populations
# from each original population (split-plot factor). This is a mixed model,
# because the error term shared among descendants of the same original population
# is a random error. Unfortunately, the glm functionality of EdgeR does not
# support random effects. Neither does DSeq2. A new package called "dream"
# (Hoffman and Roussos, 2018, http://dx.doi.org/10.1101/432567) may be useful
# fot this. For the moment, I adopted the common practice of ignoring the
# random nature of the original population effect.
#
# Thus, nested within each selective regime there are three populations of
# origin. Models 4 and 6 were specified as follows:
#
#   Model 4:  ~ Regime + Diapause + Regime:Nested_population
#   Model 6:  ~ Regime + Diapause + Regime:Diapause + Regime:Nested_population
#
# Model 6 identified only 2 genes significantly regulated by the selective
# regime, which seems very little. This model includes some interactions: zero
# genes are counted with a significant interaction term between selective regime
# and hatching condition (FDR = 0.25). But even at FDR = 0.01, 3 genes show up
# with at least one of the 4 terms of interaction between selective regime and
# population significative. I need to inspect this more closely.
#

if [ ! -e randomness.html ]; then
   R --save -q -e "rmarkdown::render('randomness.Rmd', ouptut_file='randomness.html')"
fi
