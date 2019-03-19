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
     boxplot_geom(data=A1, mapping=aes(x=1, y=V1)) +
     boxplot_geom(data=C1, mapping=aes(x=2, y=V1)) +
     boxplot_geom(data=A2, mapping=aes(x=3, y=V1)) +
     boxplot_geom(data=C2, mapping=aes(x=4, y=V1)) +
     boxplot_geom(data=A3, mapping=aes(x=5, y=V1)) +
     boxplot_geom(data=C3, mapping=aes(x=6, y=V1)) +
     boxplot_geom(data=A4, mapping=aes(x=7, y=V1)) + 
     boxplot_geom(data=C4, mapping=aes(x=8, y=V1)) + 
     boxplot_geom(data=A5, mapping=aes(x=9, y=V1)) +
     boxplot_geom(data=C5, mapping=aes(x=10, y=V1)) + 
     boxplot_geom(data=A6, mapping=aes(x=11, y=V1)) +
     boxplot_geom(data=C6, mapping=aes(x=12, y=V1)) +
     labs(x='Samples', y='Intron lengths')

ggsave('introns.png', plot=p)

