#' 01_look_see_download.r
#'
#' ## Date Written:
#'   2015-08-22
#'
#' ## Description:
#'   Find out how may transducin-like structures there are in the 
#'   protein databank (PDB) then:  
#'   - download these to disc  
#'   - extract sequences  
#'   - cluster by sequence identity  
#'   - superpose structures and cluster by RMSD  
#'   - annotate and assign sequence identity and RMSD cluster groups
#'   - compare consistency of clustering patterns with figure for collaborator.    
#'
#' ## Usage and versions:  
#'   Developed in interactive R mode then rendered with rmarkdown (i.e.)  
#'   `library(rmarkdown);  render("01_look_see_download.r")`  
#'
#'   See sessionInfo() output at end.  
#'

# Load required packages
library(bio3d)
library(dendextend)


# Use PDB structure as search seed input
input <- "1GIT"

# Extract sequence from structure and BLAST search PDB
pdb <- read.pdb(input)
aa1 <- pdbseq(pdb)
blast <- blast.pdb( pdbseq(pdb) )

# Summarize hit stat profile. 
raw.hits <- plot(blast)

# Save plot to pdf
#pdf("../results/pdb_blast_plot.pdf")
plot(blast)
#dev.off()


# Annotate significant hits
annot <- pdb.annotate(raw.hits$pdb.id, 
                      anno.terms=c("resolution", "source", 
                                   "ligandId", "citation"))

# Print table...
annot

# Save annotated hit table 
write.table(annot, sep="\t", quote=FALSE, file="../results/pdb_blast_hits.tbl")

# Extract and align sequences of all hits
raw.files <- get.pdb(raw.hits$pdb.id, split=TRUE, path="../data/")
aln <- pdbaln( raw.files )

# Print alignment
aln

# Save alignment FASTA file
write.fasta(aln, file="../results/pdb_alignment.fa")

# Produce a nice colored HTML file for this alignment
aln2html(aln, append=FALSE, file="../results/pdb_alignment.html")
aln2html(aln, colorscheme="ent", file="../results/pdb_alignment.html")


# Sequence identity analysis
ide <- seqidentity(aln)
hc <- hclust(as.dist(1-ide), method="ward.D")

#plot(hc)

# Save plot
#pdf("../results/sequence_clustering.pdf")
hclustplot(hc, k=3, labels=basename.pdb(aln$id), cex=0.8,
             ylab="1-Identity", main="Sequence Identity Clustering", fillbox=FALSE)
#dev.off()


#' ### Structure analysis  

# Fit structures
xyz <- pdbfit(aln)

# RMSD (structural distance)
rd <- rmsd(xyz)

# Clustering
hc.rmsd <- hclust(as.dist(rd), method="ward.D")
#plot(hc.rmsd)



hclustplot(hc.rmsd, k=3, labels=basename.pdb(aln$id), cex=0.8,
             ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)


# Save both plots to PDF file
pdf("../results/clustering.pdf")
hclustplot(hc, k=3, labels=basename.pdb(aln$id), cex=0.6,
             ylab="1-Identity", main="Sequence Identity Cluster Dendrogram", fillbox=FALSE)

hclustplot(hc.rmsd, k=3, labels=basename.pdb(aln$id), cex=0.6,
             ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
dev.off()


# Shorten labels for tanglegram plot
hc$labels=basename.pdb(aln$id)
hc.rmsd$labels=hc$labels
tanglegram(hc, hc.rmsd, lab.cex=0.7)


# Save plot
#pdf("../results/cluster_compare.pdf")
#tanglegram(hc, hc.rmsd)
#dev.off()

# Save all data
save.image(file="../data/01_store.RData")


# Report on software versions
sessionInfo()
