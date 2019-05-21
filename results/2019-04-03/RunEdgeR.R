# -------------------------
#  LIBRARIES AND FUNCTIONS
# -------------------------
library(edgeR)
library(ggplot2)

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
# I prefer the immediate hatching (instead of forced diapause) and the regular or
# periodic selective regime as the reference levels of the corresponding factors.

population  <- factor(c("1","1","2","2","1","1","3","3","2","2","3","3"))
diapause <- factor(rep(c("hatch","diapause"), 6))
diapause <- relevel(diapause, "hatch")
regime <- factor(c("random","random","random","random","regular","regular",
                 "random","random","regular","regular","regular","regular"))
regime <- relevel(regime, "regular")

y <- DGEList(counts=counts, group=regime)
y <- calcNormFactors(y, method="TMM")

png(filename=paste(output_dir, 'library_size.png', sep='/'))
   barplot(y$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
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
keep <- rowSums(cpm(y) > threshold) >= 4
y <- y[keep, ,keep.lib.sizes=FALSE]

png(filename=paste(output_dir, 'MDS.png', sep='/'))
   plotMDS(y$counts)
dev.off()

for (i in 1:12) {
   png(filename=paste(output_dir, "/MDnorm_", colnames(y)[i], ".png", sep=""))
   plotMD(cpm(y, log=TRUE), column=i)
   dev.off()
}

# ----------------
#  ESTIMATING BCV
# ----------------

design1 <- model.matrix(~regime + diapause + regime:population)
design2 <- model.matrix(~regime + diapause + regime:population + regime:diapause)

y1 <- estimateDisp(y, design1)
y2 <- estimateDisp(y, design2)

png(filename=paste(output_dir, 'BCV_1.png', sep='/'))
   plotBCV(y1)
dev.off()
png(filename=paste(output_dir, 'BCV_2.png', sep='/'))
   plotBCV(y2)
dev.off()

# ---------------
#  MODEL FITTING
# ---------------

fit1 <- glmQLFit(y1, design1)
fit2 <- glmQLFit(y2, design2)

png(filename=paste(output_dir, 'dispersion_1.png', sep='/'))
   plotQLDisp(fit1)
dev.off()
png(filename=paste(output_dir, 'dispersion_2.png', sep='/'))
   plotQLDisp(fit2)
dev.off()

# -----------
#  CONTRASTS
# -----------
#
# 1. Is the interaction term between selective regime and hatching condition necessary?
#    ----------------------------------------------------------------------------------

qlf_2 <- glmQLFTest(fit2, coef=8)
genes_with_interaction_term <- topTags(qlf_2, n=dim(y)[1])$table
genes_with_interaction_term$logPValue <- -log10(genes_with_interaction_term$PValue)
p <- ggplot(data=genes_with_interaction_term, mapping=aes(x=logFC, y=logPValue, color=FDR)) +
   geom_point() + theme_grey(base_size=20) + ylab(expression(paste("-log"[10], plain("(p value)"))))
ggsave(paste(output_dir, "volcano_genes_with_interaction_term.png", sep="/"), p)

NumGenes <- sum(genes_with_interaction_term$FDR <= 0.10)

if (NumGenes > 2) {
   fit <- fit2
} else {
   fit <- fit1
}

# That is, if there were more than two genes with a significant interaction between the selective
# regime and the hatching condition (at FDR = 0.1), then I would keep the more complex model,
# including that interaction. Otherwise (which is actually the case), I keep working with the
# simpler model with no interaction. Note that in both models coefficients 2 and 3 are the random
# regime and the forced diapause effects.
#
# 2. 
