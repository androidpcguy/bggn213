#' 02_make_figure.r
#'
#' ## Date Written:
#'   2015-08-22
#'
#' ## Description:
#'   Produce a nicer summary figure of sequence and structure clustering.  
#'   See '01_look_see_download.r' for data setup and exploratory analysis.  
#'
#' ## Usage and versions:  
#'   Developed in interactive R mode then rendered with rmarkdown (i.e.)  
#'   `library(rmarkdown);  render("02_make_figure.r")`  
#'
#'   See sessionInfo() output at end.    
#'

# Load required packages
library(bio3d)
library(dendextend)

# Load processed data
load("../data/01_store.RData")

# Figure setup
tanglegram(hc, hc.rmsd, lab.cex=0.6, 
	main_left ="Sequence Clustering", main_right="Structure Clustering"
	, common_subtrees_color_branches=TRUE, lwd=2)


## Lets plot a larger PDF version
pdf("../results/tanglegram_figure.pdf", width=9, height=7)

tanglegram(hc, hc.rmsd, lab.cex=0.6, 
	main_left ="Sequence Clustering", main_right="Structure Clustering"
	, common_subtrees_color_branches=TRUE, lwd=2)

dev.off()


# Report on software versions
sessionInfo()

