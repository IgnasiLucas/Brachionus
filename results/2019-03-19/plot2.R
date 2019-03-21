library(ggplot2)

A1 <- read.table('z1A_S8')
A2 <- read.table('z2A_S7')
A3 <- read.table('z3A_S9')
A4 <- read.table('z4A_S6')
A5 <- read.table('z5A_S2')
A6 <- read.table('z6A_S3')
C1 <- read.table('z1C_S1')
C2 <- read.table('z2C_S5')
C3 <- read.table('z3C_S11')
C4 <- read.table('z4C_S12')
C5 <- read.table('z5C_S4')
C6 <- read.table('z6C_S10')

p <- ggplot() + 
     geom_boxplot(data=A1, mapping=aes(x=1, y=log(V1))) +
     geom_boxplot(data=C1, mapping=aes(x=2, y=log(V1))) +
     geom_boxplot(data=A2, mapping=aes(x=3, y=log(V1))) +
     geom_boxplot(data=C2, mapping=aes(x=4, y=log(V1))) +
     geom_boxplot(data=A3, mapping=aes(x=5, y=log(V1))) +
     geom_boxplot(data=C3, mapping=aes(x=6, y=log(V1))) +
     geom_boxplot(data=A4, mapping=aes(x=7, y=log(V1))) + 
     geom_boxplot(data=C4, mapping=aes(x=8, y=log(V1))) + 
     geom_boxplot(data=A5, mapping=aes(x=9, y=log(V1))) +
     geom_boxplot(data=C5, mapping=aes(x=10, y=log(V1))) + 
     geom_boxplot(data=A6, mapping=aes(x=11, y=log(V1))) +
     geom_boxplot(data=C6, mapping=aes(x=12, y=log(V1))) +
     labs(x='Samples', y='Intron lengths')

ggsave('introns.png', plot=p)

