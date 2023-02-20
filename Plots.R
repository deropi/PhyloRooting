library(dplyr)
library(ggplot2)
library(reshape2)
library("ggpubr")


## Fig 2c
df <- read.table('Opisthokonta.df', header=TRUE, sep = ',')
dfroot<- df %>% group_by(Tree) %>% 
  filter(AD==min(AD), GeneFam=="CSC")
dfroot.2 <- as.data.frame(dcast(df, Tree ~ Partition_id, value.var="AD"))


## Fig 2b
CSC <- df[(df$GeneFam=='CSC'),]
dfroot.2 <- as.data.frame(dcast(CSC, Tree ~ Partition_id, value.var="AD"))
ggplot(dfroot.2, aes(x=`1`, y=`2`)) + geom_hex(binwidth = c(0.035, 0.035)) + coord_fixed()+ xlim(0,1)+ylim(0,1)+ geom_abline()+xlab('AD Partition 1')+ylab('AD Partition 2')#+ geom_point()+ geom_abline()


#Fig 4
df <- read.table('Proteobacteria.df', header=TRUE, sep = ',')
ggplot(df , aes(x=as.numeric(AD), color=as.character(Partition_id)))+ stat_ecdf()+theme_bw() + xlab('AD') + ylab('Cumulative distribution') + facet_wrap(vars(GeneFam))

