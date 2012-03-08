#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)

sizes <- read.csv('data/genome_size.csv')
sizes <- subset(sizes, source == "genome")

sizes$gene_count <- sizes$gene_count
sizes$genome_size <- sizes$genome_size/1000

median_density <- with(subset(sizes,strain != "R124"),{
  median(gene_count / genome_size)
})

p <- ggplot(sizes,aes(y=gene_count,x=genome_size,colour=Pseudomonas))

# Add slope line
p <- p + geom_abline(
  intercept = 0,
  slope     = median_density,
  lty       = 2,
  color     = "grey30",
  size      = 1
  )

# Annotate with median gene density
p <- p + geom_text(aes(
    x      = 5200,
    y      = 4500,
    vjust  = 0,
    hjust  = 0,
    label  = paste("median =",signif(median_density,digits=3),"Genes/KBp")),
  colour = "grey30")

# Plot size data
p <- p + geom_point(size=10)

# Highlight R124 genome
p <- p + geom_text(aes(
    x      = 6250,
    y      = 5600,
    vjust  = 1.5,
    hjust  = -0.5,
    label  = "R124"),
  colour = "grey30")

p <- p + scale_y_continuous("Gene Count")
p <- p + scale_x_continuous("Genome Size (KBp)")
p <- p + theme_bw()

postscript('out/genome_size_vs_annotations.eps',width=600)
print(p)
graphics.off()
