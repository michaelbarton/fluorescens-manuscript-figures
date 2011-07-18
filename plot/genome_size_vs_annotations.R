#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)



sizes <- read.csv('data/genome_size.csv')
sizes <- subset(sizes, source == "genome")

p <- ggplot(sizes,aes(y=gene_count/1000,x=genome_size/1000000,colour=Pseudomonas))
p <- p + geom_point(size=10)
p <- p + geom_text(aes(x = 6.25, y = 5.6, vjust = 1.5, hjust = -0.5, label = "R124"), colour = "black")
p <- p + scale_y_continuous("Gene Count (000s)")
p <- p + scale_x_continuous("Size (MBp)")
p <- p + theme_bw()

png('out/genome_size_vs_annotations.png',width=600)
print(p)
graphics.off()
