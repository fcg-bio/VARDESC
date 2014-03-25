#!/opt/well/R/R-3.0.2/bin/Rscript


# Load libraries and functions
{
	library(RColorBrewer, warn.conflicts = F, quietly = T)
	library(ggplot2, warn.conflicts = F, quietly = T)
	library(reshape, warn.conflicts = F, quietly = T)
	library(fmsb, warn.conflicts = F, quietly = T)
	library(gplots, warn.conflicts = F, quietly = T)
}


# Usage
{
usage <- "
	var-annot.R : R script to create the annotation plot from tsv
	usage : var-annot.R [options] <input file>
	Input format: Tab separated file. Additional columns between brackets
		1. Sample
		2. Gene
		[3. Group] : Grouping variable for sample
		[4. MutationCat] : String with categories for mutation. Ex: variant effect
	Options:
		help : print usage
		prefix : prefix for output files. Default = vardesc
		minMut : minimum number of samples that containts mutations per each gene. Genes with smaller number will be discarded from the plot. Default = 2
		sclust : By default samples are displayed in columns sorted by alphanumeric order. There are three strings accepted for this parameter:
			1. sclust=alphanumeric : Default
			2. sclust=group : If Group column is included in the input table, cluster the samples by group
			3. sclust=hclust : Hierarchical clustering (hclust). Dendogram will be computed and columns reordered based on columns means
	"
}


# Getting Arguments indicated by '=' and printing usage
{

arguments <- list (file=NULL, prefix="vardesc", minMut=2, sclust="alphanumeric")	
argsRaw<-commandArgs(trailingOnly = T)

if (length(argsRaw) == 0)  stop(usage)
if (argsRaw[1] == "help") stop(usage)

# Get input file name
arguments$file = argsRaw[length(argsRaw)]
argsRaw <- argsRaw[-length(argsRaw)]

if(!file.exists(arguments$file)) stop("Input file doesn't exists",usage)


# Read the other arguments
for(i in 1:length(argsRaw)){
	argsRaw2<-unlist(strsplit(argsRaw[i],"="))
	arguments[argsRaw2[1]]<-argsRaw2[2]
}


# Check arguments
if (! arguments$sclust %in% c("alphanumeric","group","hclust") ) stop("sclust input parameter not recognised",usage)
arguments$minMut <- as.numeric(arguments$minMut)


}


# Reading file
{
	#	arguments <- list (file=NULL, bySample=NULL, minMut=2, sclust="alphanumeric")
	#	arguments <- list (file=NULL, bySample=NA, minMut=2, sclust="alphanumeric", prefix = "test")
	#	arguments$file<-"annot.SNV.PRCC.tsv"

	vartab <- read.delim(arguments$file, stringsAsFactors = F)
	if(is.null(vartab$Sample)) stop("Input table doesn't contain Sample colum",usage)
	if(is.null(vartab$Group) & arguments$sclust == "group") stop("group option for sclust could not be used because input table doesn't contain Group variable",usage)

	#	vartab$Group[vartab$Sample %in% c("P21","P02","P14","P04","P05","P06","P07","P08","P09","P10","P11")] <- "Type1"
	#	vartab$Group[vartab$Sample %in% c("P12","P13","P03","P15","P16","P17","P18","P01","P19","P20","P22","P23")] <- "Type2"
}


# Quality control of the table and Modification
{
	# Delete space
	for (i in names(vartab)) vartab[,i] <- gsub(" ", "", vartab[,i])
}


# PLOT : Heatmap Gene ~ Sample
{
# Table Gene ~ Sample
htabm <- table(vartab$Gene, vartab$Sample)

# Select Genes with number of different samples mutated > arguments$minMut; sort the matrix htabmu by number of mutations
htabmu <- htabm
htabmu[htabmu > 1] <- 1
nMutGene <- rowSums(htabmu)
htabGenes <- names(sort(nMutGene[nMutGene >= arguments$minMut], decreasing = T))
htabm <- htabmu

# Subset of htabm with selected genes
htabm <- htabm[htabGenes,]

# Grouping Variables
if(!is.null(vartab$Group)) {
	groupSample <- unique(vartab[,c("Sample","Group")])
	group_char <- groupSample$Group[unlist(lapply(dimnames(htabm)[[2]], function(x) { which(x == groupSample$Sample)}))]
	group_colors <- as.character(as.numeric(as.factor(group_char)))
	group_levels <- levels(as.factor(group_char))
}

# Clustering options
{
Dendo <- "none" # Default definition = alphanumeric
ColSort <- FALSE	# Default definition = alphanumeric
if (arguments$sclust == "hclust" ) {Dendo <- "column"; ColSort <- TRUE}
if (arguments$sclust == "group" ) {
	# 1. Sort groupSample by group level	
	groupSampleSorted <- data.frame(matrix(nrow = 0, ncol = 2, dimnames = list(c(),c("Sample","Group"))))
	for (i in group_levels) {
		sorted1 <- groupSample[groupSample$Group == i, ]	# subset
		sorted2 <- sorted1[order(sorted1$Sample),]			# sort alphanumerically internally for the subset
		groupSampleSorted <- rbind(groupSampleSorted, sorted2) # append to genereral data.frame
	}
	groupSample <- groupSampleSorted
	
	# 2. Sort htabm columns by group level
	htabm <- htabm[, groupSample$Sample]
	
	# 3. Recalculate group variables
	group_char <- groupSample$Group
	group_colors <- as.character(as.numeric(as.factor(group_char)))
	group_levels <- levels(as.factor(group_char))
}

}

# Color scheme
heatCol <- brewer.pal(3, "PuBuGn") # palette for heatmap

# Build Plot
{
	Height = 0.1372549 * dim(htabm)[1]
	Width = 0.20 * dim(htabm)[2]
	
	pdfName <- paste(arguments$prefix,"annot_heatmap.pdf", sep=".")
	pdf(pdfName, title = "Mutation annotation heatmap", height = Height, width = Width)

	if(is.null(vartab$Group)) {
		heatmap.2(htabm, col=heatCol, dendrogram=Dendo, scale="none", xlab="Samples", ylab="Genes", trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepcolor='white', key=F, cexRow = 0.6, cexCol = 0.8, Rowv = F, Colv = ColSort )
	}
	if(!is.null(vartab$Group)) {
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		heatmap.2(htabm, col=heatCol, dendrogram=Dendo, scale="none", xlab="Samples", ylab="Genes", trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepcolor='white', key=F, cexRow = 0.6, cexCol = 0.8, Rowv = F, Colv = ColSort , ColSideColors=group_colors)
		legend("topleft", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
	}
	
	dev.off()
	
	svgName <- paste(arguments$prefix,"annot_heatmap.svg", sep=".")
	svg(svgName, height = Height, width = Width)
	
	if(is.null(vartab$Group)) {
		heatmap.2(htabm, col=heatCol, dendrogram=Dendo, scale="none", xlab="Samples", ylab="Genes", trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepcolor='white', key=F, cexRow = 0.6, cexCol = 0.8, Rowv = F, Colv = ColSort )
	}
	if(!is.null(vartab$Group)) {
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		heatmap.2(htabm, col=heatCol, dendrogram=Dendo, scale="none", xlab="Samples", ylab="Genes", trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepcolor='white', key=F, cexRow = 0.6, cexCol = 0.8, Rowv = F, Colv = ColSort , ColSideColors=group_colors)
		legend("topleft", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
	}
	
	dev.off()
	
	
}

}




