#!/usr/bin/env Rscript

library(seqinr)
library(plyr)
library(ggplot2)

annotations <- read.csv('data/genome/annotation/gene_list.csv')

plasmid <- grep('scaffold00008',annotations$Scaffold.Name)
plasmid <- plasmid

genome <- union(grep('scaffold',annotations$Scaffold.Name),
  grep('contig00001',annotations$Scaffold.Name))
genome <- setdiff(genome, plasmid)

genes    <- length(genome)
assembly <- length(read.fasta(file = 'data/genome/assembly/assembly.fna')[[1]])

files <- dir('data/reference/gene',pattern='*',recursive=TRUE)
files <- files[grep('genome',files)]
sizes <- adply(files,c(1),.parallel=T,function(x){
  genes  <- paste('data/reference/gene',x,sep="/")
  genome <- paste('data/reference/genomes/',(strsplit(x,'\\.')[[1]][1]),".gb",sep="")

  data.frame(
   genes  = length(read.fasta(file = genes, as.string = TRUE, seqtype = "DNA")),
   genome = as.integer(strsplit(readLines(file(genome,"r"), 1),'\\s+')[[1]][3]),
   Genome = "Pseudomonas"
  )
})

sizes <- rbind(sizes,
  data.frame(X1 = 0, genes = genes, genome=assembly,Genome="R124"))[,2:4]

p <- ggplot(sizes,aes(y=genes/1000,x=genome/1000000,colour=Genome))
p <- p + geom_point(size=10)
p <- p + scale_y_continuous("Gene Count (000s)")
p <- p + scale_x_continuous("Size (MBp)")
p <- p + theme_bw()

png('out/genome_size_vs_annotations.png',width=600)
print(p)
graphics.off()
