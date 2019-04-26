library(edgeR)

# ---------------
#  READ THE DATA
# ---------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
   stop("At least one argument must be supplied: input file.", call.=FALSE)
} else if (length(args) == 1){
   args[2] = 'z1.out'
}

input_file_name <- args[1]
counts <- read.table(input_file_name, row.names=1, header=TRUE)

# Three factors: population, diapause or hatching condition, and selective regime.
# I prefer the immedaite hatching (instead of forced diapause) and the regular or
# periodic selective regime as the reference levels of the corresponding factors.
population  <- factor(rep(1:6, each=2))
diapause <- factor(rep(c("hatch","diapause"), 6))
diapause <- relevel(diapause, "hatch")
regime <- factor(c("random","random","random","random","regular","regular",
                 "random","random","regular","regular","regular","regular"))
regime <- relevel(regime, regular)
TwoFactors <- factor(paste(regime, diapause, sep='.'))
y <- DGEList(counts=counts, group=population)
y <- calcNormFactors(y)
png(filename='MDS.png')
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

design1 <- model.matrix(~group + diapause)
design2 <- model.matrix(~0 + TwoFactors)
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

qlf_diapause <- glmQLFTest(fit1, coef=7)
qlf_diapause_inefficient <- glmQLFTest(fit3, coef=3)
