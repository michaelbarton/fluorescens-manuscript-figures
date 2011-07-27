#!/usr/bin/env Rscript

library(reshape)
library(ggplot2)
library(vegan)

log.score <- function(score){
  -1 *
  log(
    score +
    min(score[score > 0] / 100))
}

create.matrix.from.adjacency <- function(adjacency){
  as.matrix(
    cast(
      melt(adjacency,id=c("match","query")),
    query ~ match,fill=0))
}

create.plot.data <- function(scaling){
  data.frame(
    x = scaling$points[,1],
    y = scaling$points[,2]
  )
}

adjacency <- read.csv("data/pili/adjacency.csv")
adjacency$score <- log.score(adjacency$score)

distances <- dist(
               create.matrix.from.adjacency(
                 adjacency))

# Adjust 0 values to very small value
# See https://stat.ethz.ch/pipermail/r-help/2006-April/103909.html
distances[distances == 0] <- 0.000001

scaling <- metaMDSiter(distances)

p <- ggplot(create.plot.data(scaling),aes(x=x,y=y))
p <- p + geom_point()
p <- p + theme_bw()
p <- p + scale_x_continuous("First Scaling Dimension")
p <- p + scale_y_continuous("Second Scaling Dimension")
