# -------------------------
#  LIBRARIES AND FUNCTIONS
# -------------------------
library(edgeR)
library(ggplot2)

# The following function is taken from Kamil Slowikowski: https://gist.github.com/slowkow/9041570
# It creates a quantile-quantile plot for p-values, assumed to be uniformly distributed,
# and using ggplot2. Confidence intervals assume independence among tests. A qqplotr package
# exists as an extension to ggplot2 for Q-Q plots. However, it is not available through conda
# yet.

gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)
}

# ---------------
#  READ THE DATA
# ---------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
   stop("At least one argument must be supplied: input file.", call.=FALSE)
} else if (length(args) == 1) {
   args[2] <- "."
}

input_file_name <- args[1]
output_dir <- args[2]
counts <- read.table(input_file_name, row.names=1, header=TRUE)

# Three factors: population, diapause or hatching condition, and selective regime.
# I prefer the immedaite hatching (instead of forced diapause) and the regular or
# periodic selective regime as the reference levels of the corresponding factors.
population  <- factor(rep(1:6, each=2))
diapause <- factor(rep(c("hatch","diapause"), 6))
diapause <- relevel(diapause, "hatch")
regime <- factor(c("random","random","random","random","regular","regular",
                 "random","random","regular","regular","regular","regular"))
regime <- relevel(regime, "regular")
TwoFactors <- factor(paste(regime, diapause, sep='.'))
y <- DGEList(counts=counts, group=population)
y <- calcNormFactors(y)
png(filename=paste(output_dir, 'MDS.png', sep='/'))
   plotMDS(y$counts)
dev.off()

# -----------
#  FILTERING
# -----------
#
# The recommended filtering consists on removing genes that have such a low expression
# level that do not have at least 5-10 counts on as many samples as groups there are. The
# rationale is that those genes cannot have enough reads in at least one sample of each
# group. However, instead of filtering on the raw counts, the edgeR manual suggests to
# translate whatever minimum number of reads per library in terms of counts per million.
# The reason is that the filtering should take into account differences in library sizes
# between samples. It probably makes sense, because otherwise, I would be applying an
# effectively higher (harsher) minimum threshold to the gene expression level in samples
# with smaller libraries. That would probably bias results, since we would be enriching
# the dataset in genes more highly expressed in some samples than in others.
#
# Say that library sizes range between 250000 and 2000000. Then, applying a minimum expression
# level of 20 counts per million is the only way to make sure that even the small libraries
# are required to contain at least 5 counts, while at the same time setting the cutoff at the
# same biologically relevant level. In any case, as long as the minimum expression
# level is not enforced in all samples, the dataset retained will contain genes that do not
# reach the minimum in some samples. Since counts per million are (expectedly and actually)
# similarly distributed among samples, it should not be the case that any sample is enriched
# in low-expression genes.
#

threshold <- 5.0 / (min(y$samples$lib.size) / 1000000)
keep <- rowSums(cpm(y) > threshold) >= 6
y <- y[keep, ,keep.lib.sizes=FALSE]

# ----------------
#  ESTIMATING BCV
# ----------------

design1 <- model.matrix(~population + diapause)
design2 <- model.matrix(~0 + TwoFactors)
colnames(design2) <- levels(TwoFactors)
design3 <- model.matrix(~regime + diapause + regime:diapause)
y1 <- estimateDisp(y, design1)
y2 <- estimateDisp(y, design2)
y3 <- estimateDisp(y, design3)

# ---------------
#  MODEL FITTING
# ---------------

fit1 <- glmQLFit(y1, design1)
fit2 <- glmQLFit(y2, design2)
fit3 <- glmQLFit(y3, design3)

# -----------
#  CONTRASTS
# -----------
#
# 1. Genes that respond to the induction of diapause
#    -----------------------------------------------
#
# First, we use design 1 to find genes that respond to the enforcement of diapause, taking
# into account the variation among the 6 starting populations. This also takes care of the
# variation between selective regimes.

qlf_diapause_1 <- glmQLFTest(fit1, coef=7)
zplot <- gg_qqplot(qlf_diapause_1$table$PValue) + theme_grey(base_size=24)
ggsave(paste(output_dir,"qq_diapause_1.png",sep="/"), zplot)

diapause_1 <- topTags(qlf_diapause_1, n=dim(y)[1])$table
write.table(diapause_1, paste(output_dir,'diapause_design1.txt',sep='/'), quote=FALSE, sep='\t')
diapause_1$logPValue <- -log10(diapause_1$PValue)
zplot <- ggplot(data=diapause_1, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir,"volcano_diapause_1.png",sep="/"), zplot)

# We could also search for genes that respond to the induction of diapause using design 2,
# which ignores the population blocks, but takes selective regime into account. By ignoring
# the variance due to populations, this test is expected to be less sensitive. I see two ways 
# to do this. One is to look for genes that are up- or down-regulated by induction of diapause
# in both selective regimes. That is, I compare first random.diapause and random.hatch, and 
# also regular.diapause and regular.hatch, and then I identify the genes that are found
# significant in both comparisons.

qlf_random.diapause_random.hatch <- glmQLFTest(fit2, contrast=c(1,-1,0,0))
qlf_regular.diapause_regular.hatch <- glmQLFTest(fit2, contrast=c(0,0,1,-1))
qlf_diapause_2_common <- merge(
   topTags(qlf_random.diapause_random.hatch, n=dim(y)[1])$table,
   topTags(qlf_regular.diapause_regular.hatch, n=dim(y)[1])$table,
   by='row.names', all=FALSE, sort=FALSE)
#qlf_diapause_2_common <- qlf_diapause_2_common[qlf_diapause_2_common$logFC.x * qlf_diapause_2_common$logFC.y > 0,]
write.table(qlf_diapause_2_common, paste(output_dir, 'diapause_design2_common.txt', sep="/"), quote=FALSE, sep='\t')

# But, I would like a single list of p values for all genes, instead of having to select the few
# that exceed an arbitrary threshold in both comparisons. This is different from the classic
# problem of combining p values in a meta-analysis, where all p values test the same hypothesis,
# and one significant rejection is enough. I am testing two different hypotheses for each gene,
# and want to measure the probability of observing both p-values under the null hypothesis that
# at least one of the underlying hypotheses is true. I can't figure out how to do that, but I
# think to use the largest of the two p values is a sensitive choice.

qlf_diapause_2_common_pvalues <- pmax(
   qlf_random.diapause_random.hatch$table$PValue,
   qlf_regular.diapause_regular.hatch$table$PValue
)
zplot <- gg_qqplot(qlf_diapause_2_common_pvalues) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_diapause_2_common.png", sep="/"), zplot)

# The second way is to compare the average diapause effect between the two selective regimes
# with the average effect of the immediate hatching. 

qlf_diapause_2_mean <- glmQLFTest(fit2, contrast=c(0.5, -0.5, 0.5, -0.5))
zplot <- gg_qqplot(qlf_diapause_2_mean$table$PValue) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_diapause_2_mean.png", sep="/"), zplot)

diapause_2_mean <- topTags(qlf_diapause_2_mean, n=dim(y)[1])$table
write.table(diapause_2_mean, paste(output_dir, 'diapause_design2_mean.txt', sep="/"), quote=FALSE, sep='\t')
diapause_2_mean$logPValue <- -log10(diapause_2_mean$PValue)
zplot <- ggplot(data=diapause_2_mean, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir, "volcano_diapause_2_mean.png", sep="/"), zplot)
# Finally, I can use design 3 to identify genes that respond to the induction of diapause.
# This model accounts for the effect of the selective regime, the hatching condition, and
# their interaction. Coefficient 3 expresses the additive effect of the induction of diapause.
# That is, the amount of change in expression level due to diapause that is common to both
# selective regimes.

qlf_diapause_3 <- glmQLFTest(fit3, coef=3)
zplot <- gg_qqplot(qlf_diapause_3$table$PValue) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_diapause_3.png", sep="/"), zplot)
diapause_3 <- topTags(qlf_diapause_3, n=dim(y)[1])$table
write.table(diapause_3, paste(output_dir, 'diapause_design3.txt', sep='/'), quote=FALSE, sep='\t')
diapause_3$logPValue <- -log10(diapause_3$PValue)
zplot <- ggplot(data=diapause_3, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir, "volcano_diapause_3.png", sep='/'), zplot)

# It may be worth taking a look at how the p values correlate among the three differnt ways
# of identifying genes affected by the induction of diapause.
png(filename = paste(output_dir, 'diapause_designs_pvalues.png', sep='/'))
   pairs(cbind(qlf_diapause_1$table$PValue,
               qlf_diapause_2_common_pvalues,
               qlf_diapause_2_mean$table$PValue,
               qlf_diapause_3$table$PValue),
         pch='.',
         labels=c("design 1", "d.2 common", "d.2 mean", "design 3"))
dev.off()

# 2. Genes that respond to the randomness of the selective regime
#
# To search for genes regulated by the selective regime, I could use design 1 to compare
# populations 1, 2, and 4 with populations 3, 5, and 6. But coeficients in that model are
# expressed relative to population 1, and I find it difficult to express that contrast. It
# is more natural to use designs 2 and 3.
#
# In design 2, I face the same dilemma as before. I can search for genes that respond to the
# selective regime first among populations that underwent forced diapause and then also
# among populations that didn't, and identify the genes with a significant respond in both
# cases. And I could also identify genes that on average are expressed differently among
# replicates in random environments than in replicates from regular environments. I think
# the first option is better, but I'll run both:

qlf_random.diapause_regular.diapause <- glmQLFTest(fit2, contrast=c(1,0,-1,0))
qlf_random.hatch_regular.hatch       <- glmQLFTest(fit2, contrast=c(0,1,0,-1))

# As before, I will use the largest p value, for plotting purposes. Although I will need the
# two lists of genes to select the common ones. 

qlf_randomness_2_common_pvalues <- pmax(
   qlf_random.diapause_regular.diapause$table$PValue,
   qlf_random.hatch_regular.hatch$table$PValue
)
qq_plot2 <- gg_qqplot(qlf_randomness_2_common_pvalues) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_randomness_2_common.png", sep='/'), qq_plot2)

qlf_randomness_2_common <- merge(
   topTags(qlf_random.diapause_regular.diapause, n=dim(y)[1])$table,
   topTags(qlf_random.hatch_regular.hatch, n=dim(y)[1])$table,
   by='row.names', all=FALSE, sort=FALSE)
write.table(qlf_randomness_2_common, paste(output_dir, 'randomness_design2_common.txt', sep='/'), quote=FALSE, sep='\t')

qlf_randomness_2_mean <- glmQLFTest(fit2, contrast=c(0.5, 0.5, -0.5, -0.5))
qq_plot2 <- gg_qqplot(qlf_randomness_2_mean$table$PValue) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_randomness_2_mean.png", sep="/"), qq_plot2)
randomness_2_mean <- topTags(qlf_randomness_2_mean, n=dim(y)[1])$table
write.table(randomness_2_mean, paste(output_dir, 'randomness_design2_mean.txt', sep='/'), quote=FALSE, sep='\t')
randomness_2_mean$logPValue <- -log10(randomness_2_mean$PValue)
zplot <- ggplot(data=randomness_2_mean, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir, 'volcano_randomness_2_mean.png', sep='/'), zplot)

# And finally, the relative and additive effect of the random environment is represented
# by coefficient 2 in design 3. This effect does not include the interaction with the
# hatching condition. Thus, genes with a significant coeficient 2 are systematically
# affected by randomness, irrespective of hatching condition.

qlf_randomness_3 <- glmQLFTest(fit3, coef=2)
qq_plot3 <- gg_qqplot(qlf_randomness_3$table$PValue) + theme_grey(base_size=24)
ggsave(paste(output_dir, "qq_randomness_3.png", sep='/'), qq_plot3)
randomness_3 <- topTags(qlf_randomness_3, n=dim(y)[1])$table
write.table(randomness_3, paste(output_dir, 'randomness_design3.txt', sep='/'), quote=FALSE, sep='\t')
randomness_3$logPValue <- -log10(randomness_3$PValue)
zplot <- ggplot(data=randomness_3, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir, "volcano_randomness_3.png", sep='/'), zplot)

png(filename = paste(output_dir, 'randomness_designs_pvalues.png', sep='/'))
   pairs(cbind(qlf_randomness_2_common_pvalues,
               qlf_randomness_2_mean$table$PValue,
               qlf_randomness_3$table$PValue),
         pch='.',
         labels=c("design 1", "d.2 common", "d.2 mean", "design 3"))
dev.off()

rm(zplot, threshold, counts, args, input_file_name, y)
