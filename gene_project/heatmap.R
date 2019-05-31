# Author: Akshara Balachandra
# Date: May 30, 2019

library(bio3d)

# read alignment
alignment <- read.fasta('./clustal_alignment.fst')

#rd <- rmsd(alignment$ali, fit = TRUE)
#hc <- hclust(as.dist(rd))
#hclustplot(hc, labels = alignment$id)

str(alignment)
