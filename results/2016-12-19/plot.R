# This takes the arguments that I pass with --args and saves them
# in a list. There's only one, the sample name that I will use in
# the title of the plot.

args <- commandArgs(TRUE)
Sample <- args[1]

# The input files zdel and zins are lists of numbers, but there first
# row needs to be skipped.

indels <- read.table('zindels', skip=1, col.names=c('size'))

png(filename = 'indels.png')
indels_hist <- hist(indels$size, plot=FALSE)
plot(indels_hist$mids, indels_hist$count, log='y', type='h', lwd=20, lend=2, main=Sample, xlab='Size (bp)', ylab='Frequency')
dev.off()

# Never forget the 'dev.off()' statement to close an image file.
