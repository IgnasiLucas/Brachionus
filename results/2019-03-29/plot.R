library(reshape2)
library(ggplot2)

Expected  <- read.table('genes.ExpCount.txt', header=TRUE)
Posterior <- read.table('genes.PostCount.txt', header=TRUE)
# This is a convenient way to change the shape of the data frames to long, and
# merge them.
counts <- merge(melt(Expected), melt(Posterior), by=c("genes","variable"))
names(counts) <- c("genes", "sample", "expected", "posterior")

g <- ggplot(data=counts, mapping=aes(x=expected, y=posterior)) + geom_point() + geom_smooth() + 
     facet_wrap(~sample, ncol=4)
ggsave('Expected_posterior.png', plot=g)
