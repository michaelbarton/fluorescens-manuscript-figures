#!/usr/bin/env Rscript

library(reshape)
library(ggplot2)
library(vegan)
library(yaml)

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

generate.MDS <- function(score.matrix){
  names     <- row.names(score.matrix)
  distances <- dist(score.matrix)

  # Adjust 0 values to very small value
  # See https://stat.ethz.ch/pipermail/r-help/2006-April/103909.html
  distances[distances == 0] <- 0.000001
  s <- metaMDSiter(distances)$points

  row.names(s) <- names
  colnames(s)  <- c("x","y")
  s
}

adjacency <- read.csv("data/pili/adjacency.csv")
adjacency$score <- log.score(adjacency$score)

write.csv(file="data/pili/scaling.csv",
  generate.MDS(
    create.matrix.from.adjacency(
      adjacency)))
