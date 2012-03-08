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
  c("Pf0-1","Pf-5","SBW-25")
  ))

base <- ggplot(data) + theme_bw()

p1 <- base + geom_segment(aes(y = V1/1000000, yend = V2/1000000,
  x = V3/1000000, xend = V4/1000000,color=reference),size = 6)
p1 <- p1 + scale_y_continuous("P. fluorescens Reference Genome Position (MBp)")
p1 <- p1 + xlab(NULL)
p1 <- p1 + geom_text(aes(x,y,label = text),
                     data.frame(x=0.25, y=7, text="A"))

p2 <- base + geom_density(size = 1, adjust = 1/3, aes(x = (V4 + V3)/2000000,color=reference))
p2 <- p2 + scale_x_continuous("P. fluorescens R124 Genome Position (MBp)")
p2 <- p2 + scale_y_continuous("Sequence Alignment Density",breaks = c(0))
p2 <- p2 + geom_text(aes(x,y,label = text),
                     data.frame(x=0.25, y=0.25, text="B"))

legend <- p2 + opts(keep = "legend_box")
p1 <- p1 + opts(legend.position = "none")
p2 <- p2 + opts(legend.position = "none")

Layout <- grid.layout(
  nrow    = 2,
  ncol    = 2,
  widths  = unit(c(2, 0.4), c("null", "null")),
  heights = unit(c(5, 2),   c("null", "null"))
)

vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}

subplot <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

postscript("out/genome_dot_plot.eps", width = 600)

vplayout()
print(p1,     vp = subplot(1, 1))
print(p2,     vp = subplot(2, 1))
print(legend, vp = subplot(1:2, 2))

graphics.off()
