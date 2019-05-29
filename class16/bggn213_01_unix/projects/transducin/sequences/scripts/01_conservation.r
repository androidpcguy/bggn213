#' 01_conservation.r  
#'
#' ## Date Written:  
#'   2015-08-22  
#'
#' ## Descrption:  
#'   FInd conserved residue clusters in GalphaI PFAM alignment.   
#'   First:  
#'   - Read alignment   
#'   - Calculate conservation  
#'   - Calculate consensus sequence  
#'   - Plot conservation score along with conserved consnsus position labels  
#'
#' ## Usage and versions:
#'   Developed in interactive R mode then rendered with rmarkdown (i.e.)
#'   `library(rmarkdown);  render("01_conservation.r")`
#'
#'   See sessionInfo() output at end.
#'

library(bio3d)

# Read PFAM alignment
aln <- read.fasta("../data/PF00503_seed.txt")

# Calculate and plot conservation score per alignment position
sim <- conserv(aln)
plot(sim, typ="h", ylab="Conservation", xlab="Alignment position")


# Zero-out low scoring positions for plot below
sim[sim < 0.6] = 0


## Calculate consensus sequence
con <- consensus(aln, cutoff = 0.8)
paste(con$seq, collapse="")

# Determine which positions score more 0.75
conserved.inds <- which(sim > 0.75)

# Plot setup development
plot(sim, typ="h", ylab="Conservation", xlab="Alignment position", col="gray20")
text(conserved.inds, sim[conserved.inds], labels=con$seq[conserved.inds], col="blue")
rug(conserved.inds, col="blue")


# Plot results to PDF
pdf("../results/conservation.pdf", width=10)
plot(sim, typ="h", ylab="Conservation", xlab="Alignment position", col="gray20")
text(conserved.inds, sim[conserved.inds], labels=con$seq[conserved.inds], col="blue")
rug(conserved.inds, col="blue")
dev.off()


#' ## Report on software versions
sessionInfo()
