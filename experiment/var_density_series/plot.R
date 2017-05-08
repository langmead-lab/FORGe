#!/usr/bin/env Rscript

library('dplyr')
library('ggplot2')

m <- read.table('tab.csv', sep=',', header=T)
m$aligned <- bitwAnd(m$flag, 4) == 0
mg <- m %>% group_by(type, stride, alts) %>% summarise(frac_al=sum(aligned)/n(), frac_cor=sum(cor == 1)/n(), frac_al_cor=sum(cor == 1)/sum(aligned))
ggplot(mg, aes(x=1/stride, color=alts, linetype=type)) + geom_line(aes(y=frac_al)) + labs(x='Density', y='Fraction')
