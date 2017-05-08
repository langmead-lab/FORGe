#!/usr/bin/env Rscript

library('dplyr')
library('ggplot2')
library('gridExtra')

m <- read.table('overall.csv', sep=',', header=T)
m$aligned <- bitwAnd(m$flag, 4) == 0
mg <- m %>% group_by(type, stride, pct) %>% summarise(frac_al=sum(aligned)/n(), frac_cor=sum(cor == 1)/n(), frac_al_cor=sum(cor == 1)/sum(aligned))
mg$stride <- factor(mg$stride)

myplot <- function() {
       frac_al     <- ggplot(mg, aes(x=pct, color=stride, linetype=type)) + geom_line(aes(y=frac_al    )) + geom_point(aes(y=frac_al    )) + labs(x='% vars included', y='Fraction', title='Fraction aligned')
       frac_al_cor <- ggplot(mg, aes(x=pct, color=stride, linetype=type)) + geom_line(aes(y=frac_al_cor)) + geom_point(aes(y=frac_al_cor)) + labs(x='% vars included', y='Fraction', title='Fraction correct of aligned')
       frac_cor    <- ggplot(mg, aes(x=pct, color=stride, linetype=type)) + geom_line(aes(y=frac_cor   )) + geom_point(aes(y=frac_cor   )) + labs(x='% vars included', y='Fraction', title='Fraction correct')
       grid.arrange(frac_al, frac_al_cor, frac_cor, ncol=3)
}

# PDF
pdf('all3.pdf', height=3, width=10)
myplot()
dev.off()

# PNG
png('all3.png', height=3, width=10, units='in', res=73)
myplot()
dev.off()
