library(ggplot2)

gt<-read.table("W55.Somatic.gt", header = T)

samples<-colnames(gt)[3:ncol(gt)]

ggplot(gt)+geom_tile(aes=)