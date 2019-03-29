library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
outputFile <- args[1]
inputFile  <- args[2]

fpkm <- read.table(inputFile)
fpkm_positive <- data.frame(V1=fpkm[fpkm$V1 > 0,])
p <- ggplot(data=fpkm_positive, mapping=aes(x=V1)) + geom_histogram(bins=100) + scale_x_log10()
ggsave(outputFile, plot=p)
