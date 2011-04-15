#!/usr/bin/env Rscript

library('plyr')
library('ggplot2')

source('lib/r/drop.levels.r')
source('lib/r/find.replace.R')

options(stringsAsFactors = FALSE)

data <- read.csv('data/alignment/nucmer.coords',sep="\t",header = FALSE)
names(data)[length(names(data)) - 1] <- "reference"
names(data)[length(names(data))]     <- "assembly"

# Select for scaffold 8 matches
data <- drop.levels( data[grep("fluorescens",data$reference),] )

# Select only those matches with at least 10 hits
data <- ddply(data,.(reference),function(df){
  if(dim(df)[1] > 10){
    return(df)
  }
})

data$reference <- as.factor(find.replace(data$reference,
  c("fluorescens_pf01_genome","fluorescens_pf5_genome","fluorescens_sbw25_genome"),
  c("Pf-01","Pf-5","SBW-25")
  ))

png("out/figure_2.png",width = 600)

p <- ggplot(data,aes(x = V1/1000000, xend = V2/1000000,
  y = V3/1000000, yend = V4/1000000,color=reference))
p <- p + geom_segment(size = 4)
p <- p + scale_x_continuous("Reference Genome Position (MBp)")
p <- p + scale_y_continuous("R124 Genome Position (MBp)")
p <- p + theme_bw()

print(p)
graphics.off()
