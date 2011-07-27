#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)

p <- ggplot(
  read.csv("data/pili/plot_data.csv"),
  aes(x=x,y=y,color=name))

p <- p + geom_point()

p <- p + scale_x_continuous("Second Dimensional Scaling",limits=c(-10,10))
p <- p + scale_y_continuous("First Dimensional Scaling",limits=c(-10,10))
p <- p + scale_colour_manual("Genes",values =c("black",brewer.pal(9,"Set1")))

print(p)
