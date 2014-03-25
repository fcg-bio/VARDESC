#!/opt/well/R/R-3.0.2/bin/Rscript

# Usage
{
usage <- "
	var-context.R : R script to create mutation context plots from tsv
	usage : var-context.R [options] <input file>
	example : var-context.R prefix=outfile BSgenome=BSgenome.Hsapiens.UCSC.hg19 fasta=genome.fa bySample file.tsv
	
	Input format: Tab separated file with the following names (the order is not important). Additional columns between brackets
		1. Chromosome
		[2. Start]
		[3. End]
		4. Ref
		5. Alt
		[6. Sample] : only necessary if bySample option is selected
		[7. Group] : Grouping variable for sample. If defined, the heatmap will contain the a color coloumn with the grouping info
		8. Upstream : String with upstream (5') region
		9. Downstream : String with downstream (3') region
		
	Options:
		help : print usage
		prefix : prefix for output files. Default = vardesc
		bySample : plot by sample
		BSgenome : name of the genome as described in th R library BSgenome. By default BSgenome.Hsapiens.UCSC.hg19. Other genomes must be installed
		fasta : fasta file containing the genome to analyse. In case that BSgenome is not installed or Exome analysis. Default = NULL
		
Example of exome analysis :
#1. Create Targeted Fasta from a bed, gff or vcf file
bedtools getfasta -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>
#2. Run var-context
var-context.R prefix=outfile bySample fasta=genome.fa file.tsv

	"
}


# Getting Arguments indicated by '=' and printing usage
{

arguments <- list (file=NULL, bySample=NULL, prefix="vardesc", BSgenome="BSgenome.Hsapiens.UCSC.hg19", fasta=NULL)	
argsRaw<-commandArgs(trailingOnly = T)

if (length(argsRaw) == 0)  stop(usage)
if (argsRaw[1] == "help") stop(usage)

# Get input file name
arguments$file = argsRaw[length(argsRaw)]
argsRaw <- argsRaw[-length(argsRaw)]

if(!file.exists(arguments$file)) stop("Input file doesn't exists",usage)
if(!is.null(arguments$fasta) & !file.exists(arguments$file)) stop("Input fasta file doesn't exists",usage)
if(is.null(arguments$BSgenome) & is.null(arguments$fasta)) stop("Please, provide either BSgenome or fasta options",usage)


# Read the other arguments
for(i in 1:length(argsRaw)){
	argsRaw2<-unlist(strsplit(argsRaw[i],"="))
	arguments[argsRaw2[1]]<-argsRaw2[2]
}

#	arguments <- list (file=NULL, bySample=NA, prefix="test", BSgenome="BSgenome.Hsapiens.UCSC.hg19")

}


# Load libraries and functions
{
#source("/users/fcastro/R/functions/functions.R")
#source("~/R/functions/functions.R")
##source("/well/tomlinson/fcastro/devel/VARDESC/var-desc-functions.R")
#source("var-desc-functions.R")

library(RColorBrewer, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(reshape, warn.conflicts = F, quietly = T)
library(fmsb, warn.conflicts = F, quietly = T)
library(gplots, warn.conflicts = F, quietly = T)
library(pheatmap, warn.conflicts = F, quietly = T)
library(Biostrings, warn.conflicts = F, quietly = T)
library(BSgenome, warn.conflicts = F, quietly = T)
if(!is.null(arguments$BSgenome)) library(arguments$BSgenome, warn.conflicts = F, quietly = T ,character.only = TRUE)

}


# Reading files
{
#	arguments <- list (file=NULL, bySample=NULL, prefix="test")
#	arguments <- list (file=NULL, bySample=NA, prefix="test", fasta="test.bed.fa")
#	arguments <- list (file=NULL, bySample=NA, prefix="test")
#	arguments$file<-"context.SNV.PRCC.tsv"

vartab <- read.table(arguments$file, stringsAsFactors = F, header = T)

#	vartab$Group[vartab$Sample %in% c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11")] <- "Type1"
#	vartab$Group[vartab$Sample %in% c("P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23")] <- "Type2"
}


# Loading Genomes
{

if(is.null(arguments$fasta)) {
	BSgenome.package <- grep("BSgenome.", search(), value = T)
	GenomeName <- unlist(strsplit(BSgenome.package, "\\."))[2]
	Genome <- get(GenomeName)
}
if(!is.null(arguments$fasta)) {
	Genome <- readDNAStringSet(arguments$fasta, "fasta")
}

}


# Get Number of Samples
{
	if(is.null(arguments$bySample)) {
		nsamples <- 1
	}

	if(!is.null(arguments$bySample)) {
		nsamples <- nlevels(factor(vartab$Sample))
	}

}


# Quality control of the table and Modification
{

	# Delete space
	for (i in names(vartab)) vartab[,i] <- gsub(" ", "", vartab[,i])

	# Select only SNPs
	snptab <- vartab
	snptab <- snptab[nchar(snptab$Ref)==1 & snptab$Ref!="-", ]
	snptab <- snptab[nchar(snptab$Alt)==1 & snptab$Alt!="-", ]
	
	# Detect number of bases in upstream and downstream
	nUp <- unique(nchar(snptab$Upstream))
	nDown <- unique(nchar(snptab$Downstream))
	
	if(length(nUp) > 1 | length(nDown) > 1 ) stop("Different length of context string")
	
	# Modify Ref, Alt, Upstream and Downstream variables according to reference strand used (Biostrings)
	snptab$cRef <- snptab$Ref
	snptab$cRef[snptab$Ref %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(snptab$Ref[snptab$Ref %in% c("A","G")])))
	snptab$cAlt <- snptab$Alt
	snptab$cAlt[snptab$Ref %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(snptab$Alt[snptab$Ref %in% c("A","G")])))
	
	snptab$cUpstream <- snptab$Upstream
	snptab$cUpstream[snptab$Ref %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(snptab$Upstream[snptab$Ref %in% c("A","G")])))
	snptab$cDownstream <- snptab$Downstream
	snptab$cDownstream[snptab$Ref %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(snptab$Downstream[snptab$Ref %in% c("A","G")])))
		
	# Create String of reference sequence
	snptab$RefRegion <- paste(snptab$cUpstream, snptab$cRef, snptab$cDownstream, sep ="")
	snptab$RefRegion[snptab$Ref %in% c("A","G")] <- paste(snptab$cDownstream[snptab$Ref %in% c("A","G")], snptab$cRef[snptab$Ref %in% c("A","G")], snptab$cUpstream[snptab$Ref %in% c("A","G")], sep ="")
	snptab$AltRegion <- paste(snptab$cUpstream, snptab$cAlt, snptab$cDownstream, sep ="")
	snptab$AltRegion[snptab$Ref %in% c("A","G")] <- paste(snptab$cDownstream[snptab$Ref %in% c("A","G")], snptab$cAlt[snptab$Ref %in% c("A","G")], snptab$cUpstream[snptab$Ref %in% c("A","G")], sep ="")	

	# Create Variant columns
	snptab$Variant <- paste(snptab$cRef, snptab$cAlt, sep =">") 
	snptab$VariantRegion <- paste(snptab$Variant, snptab$RefRegion, sep =":") 
	
	# Dinucleotide
	Np <- rep(NA, dim(snptab)[1])
	Np[snptab$Ref %in% c("C","T")] <- substr(snptab$cUpstream[snptab$Ref %in% c("C","T")], nchar(snptab$cUpstream[1]), nchar(snptab$cUpstream[1]))
	Np[snptab$Ref %in% c("A","G")] <- substr(snptab$cUpstream[snptab$Ref %in% c("A","G")], 1, 1)
	snptab$Dinucleotide <- paste(Np, snptab$cRef, sep="")

	
#	# Add Category All into Sample column is bySample option is specified
#	if(!is.null(arguments$bySample)) {
#		alltab <- snptab
#		alltab$Sample <- "All"
#		snptab <- rbind(snptab,alltab)
#	}

}


# Calculate oligonucleotide distribution in the reference genome
{

context.length <- 1 + nchar(vartab$Upstream[1]) + nchar(vartab$Downstream[1])

# Calculate the frequency of each context string in the reference genome + dinucleotide composition

	# BSgenome
	if (is.null(arguments$fasta)) {
		Chrs <- names(Genome)[grep("random|chrUn|upstream|_", names(Genome), invert = T)]
		for (chr in Chrs) {
			Seq <- Genome[[chr]]
			onf <- oligonucleotideFrequency(Seq, context.length)
			dnf <- dinucleotideFrequency(Seq)

			data.onf <- data.frame(RefRegion = names(onf), RefRegionCounts = onf)
			data.dnf <- data.frame(RefRegion = names(dnf), RefRegionCounts = dnf)

			if(exists("ndist")) {
				nodist <- merge(nodist, data.onf, by="RefRegion", all=T)
				nodist$RefRegionCounts <- nodist$RefRegionCounts.x + nodist$RefRegionCounts.y
				nodist <- nodist[,-c(2,3)]
				
				ndist <- merge(ndist, data.dnf, by="RefRegion", all=T)
				ndist$RefRegionCounts <- ndist$RefRegionCounts.x + ndist$RefRegionCounts.y
				ndist <- ndist[,-c(2,3)]
			}
			if(!exists("ndist")) {nodist <- data.onf; ndist <- data.dnf}
		}
	}

	# Fasta
	if (!is.null(arguments$fasta)) {
		onf <- colSums(oligonucleotideFrequency(Genome, context.length))
		nodist <- data.frame(RefRegion=names(onf), RefRegionCounts=onf)
		
		dnf <- colSums(dinucleotideFrequency(Genome))
		ndist <- data.frame(RefRegion=names(dnf), RefRegionCounts=dnf)
	}

# Make reverse complementary for those positions = C or T in the middle
varPos <- round(context.length/2)	# get the index for the position in the middle = variant position
nodist$varPos <- substr(nodist$RefRegion, varPos, varPos)
nodist$RefRegion <- as.character(nodist$RefRegion)
nodist$RefRegion[nodist$varPos %in% c("A","G")] <- as.character(reverseComplement(DNAStringSet(nodist$RefRegion[nodist$varPos %in% c("A","G")])))

# Recalculate proportion + normalization factor
Counts <- c()
for( StrSeq in levels(factor(nodist$RefRegion))) {
	Counts <- append(Counts, sum(nodist$RefRegionCounts[nodist$RefRegion == StrSeq]))
}
Counts <- as.numeric(Counts)
nodist <- data.frame(RefRegion = levels(factor(nodist$RefRegion)), RefRegionCounts = Counts, RefRegionProp = Counts/ sum(Counts))
nodist$varPos <- substr(nodist$RefRegion, varPos, varPos)
evenProp <- 1 / dim(nodist)[1]
nodist$RefRegionNorm <- nodist$RefRegionProp / evenProp 

}


# Adding to snptab the oligonucleotide distribution in the reference genome
{
	snptab <- merge(snptab, nodist, by = "RefRegion", all.x = T)
	snptab$RefRegionCounts[is.na(snptab$RefRegionCounts)] <- 0
	snptab$RefRegionProp[is.na(snptab$RefRegionProp)] <- 0
}


# PLOT : Barplot mutation context 
{
Height <- 3
if (nsamples != 1) Height <- nsamples

tmp <- unique(snptab[,c("VariantRegion","Variant","RefRegion", "RefRegionCounts", "RefRegionProp", "RefRegionNorm")])

if(is.null(arguments$bySample)) {
	varprop <- data.frame(prop.table(table(snptab$VariantRegion)))
	names(varprop) <- c("VariantRegion", "Frequency")
	varprop <- merge(varprop,tmp)
	
	# Vertical Grid coordinates for Variant type
	vGridCoord <- as.numeric(cumsum(table(varprop$Variant)))
	vGridCoord <- vGridCoord[-length(vGridCoord)] + 0.5

	if(!is.null(snptab$Group)) {
		varprop <- data.frame(prop.table(table(snptab$VariantRegion, snptab$Group)))
		names(varprop) <- c("VariantRegion", "Group", "Frequency")
		varprop <- merge(varprop,tmp)
		# Vertical Grid coordinates for Variant type
		vGridCoord <- as.numeric(cumsum(table(varprop$Variant[varprop$Group==levels(varprop$Group)[1]])))
		vGridCoord <- vGridCoord[-length(vGridCoord)] + 0.5
	}
	
	# Normalize the Frequency
	varprop$normFrequency <- varprop$Frequency / varprop$RefRegionNorm
	
	# Build the plot
	out <- ggplot(data=varprop
		, aes(x=VariantRegion, y=Frequency, fill=Variant, width=0.6, main="Patterns of substitutions")) +
		geom_bar(stat="identity") + 
		theme_bw() + 
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			, axis.text.x  = element_text(angle=90, vjust=0.5, size=4)
		) + 
		scale_fill_brewer(palette = "Set1") + 
		ylab("Mutation Proportion") + 
		xlab("Samples") + 
		geom_vline(xintercept=vGridCoord, linetype="dashed", colour="grey")
	
	if(!is.null(snptab$Group)) out <- out + facet_grid(Group ~ .)
	
	# Build the plot Normalized proportion
	outNorm <- ggplot(data=varprop
		, aes(x=VariantRegion, y=normFrequency, fill=Variant, width=0.6, main="Patterns of substitutions")) +
		geom_bar(stat="identity") + 
		theme_bw() + 
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			, axis.text.x  = element_text(angle=90, vjust=0.5, size=4)
		) + 
		scale_fill_brewer(palette = "Set1") + 
		ylab("Normalized Mutation Proportion") + 
		xlab("Samples") + 
		geom_vline(xintercept=vGridCoord, linetype="dashed", colour="grey")
		
	if(!is.null(snptab$Group)) outNorm <- outNorm + facet_grid(Group ~ .)
}

if(!is.null(arguments$bySample)) {
	varprop <- prop.table(table(snptab$VariantRegion, snptab$Sample), 2)
	varprop <- data.frame(varprop)
	names(varprop) <- c("VariantRegion", "Sample", "Frequency")
	varprop <- merge(varprop,tmp)

	# Normalize the Frequency
	varprop$normFrequency <- varprop$Frequency / varprop$RefRegionNorm
	
	# Vertical Grid coordinates for Variant type
	vGridCoord <- as.numeric(cumsum(table(varprop$Variant[varprop$Sample==varprop$Sample[1]])))
	vGridCoord <- vGridCoord[-length(vGridCoord)] + 0.5

	# Build the plot	
	out <- ggplot(data=varprop
		, aes(x=VariantRegion, y=Frequency, fill=Variant, width=0.6, main="Patterns of substitutions")) +
		geom_bar(stat="identity") + 
		theme_bw() + 
		facet_grid(Sample ~ .) +
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			, axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
		) + 
		scale_fill_brewer(palette = "Set1") + 
		ylab("Mutation Proportion") + 
		xlab("Samples")  + 
		geom_vline(xintercept=vGridCoord, linetype="dashed", colour="grey")

	# Build the plot Normalized proportion
	outNorm <- ggplot(data=varprop
		, aes(x=VariantRegion, y=normFrequency, fill=Variant, width=0.6, main="Patterns of substitutions")) +
		geom_bar(stat="identity") + 
		theme_bw() + 
		facet_grid(Sample ~ .) +
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			, axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
		) + 
		scale_fill_brewer(palette = "Set1") + 
		ylab("Normalized Mutation Proportion") + 
		xlab("Samples")  + 
		geom_vline(xintercept=vGridCoord, linetype="dashed", colour="grey")
	
}

pdfName <- paste(arguments$prefix,"context_barplot.pdf", sep=".")
pdf(pdfName, width= 12, height= Height, title = "Mutation context barplot")
	print(out)
dev.off()
pdfName <- paste(arguments$prefix,"context_normalized_barplot.pdf", sep=".")
pdf(pdfName, width= 12, height= Height, title = "Normalized Mutation context barplot")
	print(outNorm)
dev.off()

svgName <- paste(arguments$prefix,"context_barplot.svg", sep=".")
svg(svgName, width= 12, height= Height)
	print(out)
dev.off()
svgName <- paste(arguments$prefix,"context_normalized_barplot.svg", sep=".")
svg(svgName, width= 12, height= Height)
	print(outNorm)
dev.off()

rm(tmp)
}


# TABLE : Statistical differences in mutation spoutectrum by sample vs background
{
	write.table(varprop, file = paste(arguments$prefix,"context_table.txt", sep="."), 
		quote = F, sep = "\t", row.names = F)

}


# PLOT : Heatmap mutation context 
{
if(!is.null(arguments$bySample)) {

	Height <- 0.15*nlevels(factor(snptab$VariantRegion))
	Width1 <- 0.15*nlevels(factor(snptab$VariantRegion)) / 3
	Width2 <- 0.5*nlevels(factor(snptab$Sample))
	Width <- max(Width1, Width2)
	
	tmp <- unique(snptab[,c("VariantRegion","Variant","RefRegion", "RefRegionCounts", "RefRegionProp", "RefRegionNorm")])

	varprop <- prop.table(table(snptab$VariantRegion, snptab$Sample), 2)

	htabm <- data.matrix(varprop)

	# Normalize htabm
	htabmNorm <- htabm
	for (i in row.names(htabm)) {
		normFactor <- tmp[tmp$VariantRegion==i,"RefRegionNorm"]
		htabmNorm[row.names(htabmNorm) == i, ] <- htabm[row.names(htabm) == i, ] / normFactor
	}
	rm(tmp)

	# Group color scheme
	if(!is.null(vartab$Group)) {
		groupSample <- unique(vartab[,c("Sample","Group")])
		group_char <- groupSample$Group[unlist(lapply(dimnames(varprop)[[2]], function(x) { which(x == groupSample$Sample)}))]
		group_colors <- as.character(as.numeric(as.factor(group_char)))
		group_levels <- levels(as.factor(group_char))
	}

	# Color scheme

#	heatCol <- brewer.pal(8, "PuBuGn") # palette for heatmap
	heatCol <- rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlGnBu")))(100))

	# Create the plot
	pdfName <- paste(arguments$prefix,"context_heatmap.pdf", sep=".")
	pdf(pdfName, width= Width, height= Height, title = "Mutation context barplot")
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		if(!is.null(vartab$Group)) {
			heatmap.2(htabm, col=heatCol, dendrogram="column", scale="column", ColSideColors=group_colors, xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
			legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
		}

		if(is.null(vartab$Group)) { 
			pheatmap(htabm, border_color=NA, col=heatCol)
#			heatmap.2(htabm, col=heatCol, dendrogram="column", scale="column", xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
		}
	dev.off()

	svgName <- paste(arguments$prefix,"context_heatmap.svg", sep=".")
	svg(svgName, width= Width, height= Height)
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		if(!is.null(vartab$Group)) {
			heatmap.2(htabm, col=heatCol, dendrogram="column", scale="column", ColSideColors=group_colors, xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
			legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
		}

		if(is.null(vartab$Group)) {
			pheatmap(htabm, border_color=NA, col=heatCol)
#			heatmap.2(htabm, col=heatCol, dendrogram="column", scale="column", xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
		}
	dev.off()

	# Create the Normalized plot
	pdfName <- paste(arguments$prefix,"context_heatmap.pdf", sep=".")
	pdf(pdfName, width= Width, height= Height, title = "Mutation context barplot")
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		if(!is.null(vartab$Group)) {
			heatmap.2(htabmNorm, col=heatCol, dendrogram="column", scale="column", ColSideColors=group_colors, xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabmNorm)[2]), rowsep=c(1:dim(htabmNorm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
			legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
		}

		if(is.null(vartab$Group)) {
			pheatmap(htabmNorm, border_color=NA, col=heatCol)
#			heatmap.2(htabmNorm, col=heatCol, dendrogram="column", scale="column", xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabmNorm)[2]), rowsep=c(1:dim(htabmNorm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
		}
	dev.off()

	svgName <- paste(arguments$prefix,"context_heatmap.svg", sep=".")
	svg(svgName, width= Width, height= Height)
		palette(brewer.pal(8,"Set1")) # palette for grouping scheme
		if(!is.null(vartab$Group)) {
			heatmap.2(htabmNorm, col=heatCol, dendrogram="column", scale="column", ColSideColors=group_colors, xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabmNorm)[2]), rowsep=c(1:dim(htabmNorm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
			legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
		}
		if(is.null(vartab$Group)) {
			pheatmap(htabmNorm, border_color=NA)
#			heatmap.2(htabmNorm, col=heatCol, dendrogram="column", scale="column", xlab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabmNorm)[2]), rowsep=c(1:dim(htabmNorm)[1]), sepwidth=c(0.00005, 0.005), sepcolor='white', keysize=1)
		}
	dev.off()


}

}


# PLOT : CpG Dinucleotide
{
# GpC proportion in the genome
gCpG <- sum(ndist$RefRegionCounts[grep("CG", ndist$RefRegion)])
gNpG <- sum(ndist$RefRegionCounts[grep("CG|AG|TG", ndist$RefRegion)])
gCpGprop <- gCpG / (gCpG + gNpG)
	
if(is.null(arguments$bySample) & is.null(vartab$Group)) {

	tmp <- snptab
	tmp$CpG[tmp$Ref=="G"] <- "NpG"
	tmp$CpG[tmp$Ref=="G" & tmp$Dinucleotide == "GC"] <- "CpG"
	
	varprop <- data.frame(table(tmp$Variant[tmp$Ref=="G"], tmp$CpG[tmp$Ref=="G"]))
	varprop <- cast(varprop, Var1~Var2, value = "Freq")
	names(varprop)[1] <- c("Variant")
	varprop$NpG <- varprop$CpG + varprop$NpG
	varprop$CpGprop <- varprop$CpG / varprop$NpG
	
	varprop$pvalue <- NA
	for (i in 1:dim(varprop)[1]) {
		x <- varprop[i,c("NpG","CpG")]
		m <- matrix(c(gNpG, gCpG, x$NpG, x$CpG), nrow=2, byrow=T, dimnames=list(c("Genome","Query"),c("NpG","CpG")))
		varprop$pvalue[i] <- chisq.test(m)$p.value
	}
	pvalR <- round(varprop$pvalue,4)
	pvalR <- paste ("P", pvalR, sep="=")
	pvalR[varprop$pvalue < 0.0001] <- "P<0.0001"
	varprop$pvalue <- pvalR
	
	varprop$Variant <- as.character(varprop$Variant)
	varprop <- rbind(c("Genome", gCpG, gNpG, gCpGprop, ""), varprop)
	varprop$Variant[varprop$Variant == "C>A"] <- "G>T\nC>A"
	varprop$Variant[varprop$Variant == "C>T"] <- "G>A\nC>T"
	varprop$Variant[varprop$Variant == "C>G"] <- "G>C\nC>G"
	varprop$Variant <- factor(varprop$Variant, levels=varprop$Variant)
	varprop$CpGprop <- as.numeric(varprop$CpGprop)

	out <-  ggplot(data=varprop, aes(x=Variant, y=CpGprop, fill=Variant, width=0.6)) +
		geom_bar(stat="identity")  + 
		ylab("Fraction of CpG per NpG") +
		theme_bw() + 
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			,legend.position="none"
		) + 
		scale_fill_brewer(palette = "Set1") +
		geom_text(data=varprop, 
            aes(Variant, CpGprop, label=varprop$pvalue),
            size = 4, vjust=0) + ggtitle("Fraction of guanine mutations\noccurring at CpG dinucleotides")
            
   # Create the plot
	pdfName <- paste(arguments$prefix,"context_dinucleotide.pdf", sep=".")
	pdf(pdfName, title = "Fraction of the three classes of guanine mutations occurring at CpG dinucleotides")
		print(out)
	dev.off()
	
	svgName <- paste(arguments$prefix,"context_dinucleotide.svg", sep=".")
	svg(svgName)
		print(out)
	dev.off()
	
	
}

if(is.null(arguments$bySample) & !is.null(vartab$Group)) {

	tmp <- snptab
	tmp$CpG[tmp$Ref=="G"] <- "NpG"
	tmp$CpG[tmp$Ref=="G" & tmp$Dinucleotide == "GC"] <- "CpG"
	
	varprop <- data.frame(table(tmp$Variant[tmp$Ref=="G"], tmp$CpG[tmp$Ref=="G"], tmp$Group[tmp$Ref=="G"] ))
	names(varprop) <- c("Variant","Dn","Group", "Freq")
	for(i in levels(factor(varprop$Group))) {
		if(exists("tmpCast")) { 
			x <- cast(varprop[varprop$Group==i,], Variant~Dn, value = "Freq")
			x$Variant <- as.character(x$Variant)
			x <- rbind(x, c("Genome", gCpG, gNpG))
			x$Group <- i
			tmpCast <- rbind(tmpCast, x)
		}
		if(!exists("tmpCast")) {
			tmpCast <- cast(varprop[varprop$Group==i,], Variant~Dn, value = "Freq")
			tmpCast$Variant <- as.character(tmpCast$Variant)
			tmpCast <- rbind(tmpCast, c("Genome", gCpG, gNpG))
			tmpCast$Group <- i
		}

	}
	varprop <- tmpCast
	rm(tmpCast)

	varprop$CpG <- as.numeric(varprop$CpG)
	varprop$NpG <- as.numeric(varprop$NpG)
	varprop$NpG <- varprop$CpG + varprop$NpG
	varprop$CpGprop <- varprop$CpG / varprop$NpG
	varprop$CpGprop[is.na(varprop$CpGprop)] <- 0
	
	varprop$pvalue <- NA
	for (i in 1:dim(varprop)[1]) {
		x <- varprop[i,c("NpG","CpG")]
		m <- matrix(c(gNpG, gCpG, x$NpG, x$CpG), nrow=2, byrow=T, dimnames=list(c("Genome","Query"),c("NpG","CpG")))
		varprop$pvalue[i] <- chisq.test(m)$p.value
	}
	pvalR <- round(varprop$pvalue,4)
	pvalR <- paste ("P", pvalR, sep="=")
	pvalR[varprop$pvalue < 0.0001] <- "P<0.0001"
	pvalR[varprop$pvalue == 0] <- ""
	varprop$pvalue <- pvalR
	
	varprop$Variant[varprop$Variant == "C>A"] <- "G>T\nC>A"
	varprop$Variant[varprop$Variant == "C>T"] <- "G>A\nC>T"
	varprop$Variant[varprop$Variant == "C>G"] <- "G>C\nC>G"
	varprop$Variant <- factor(varprop$Variant, levels=c("Genome", "G>T\nC>A", "G>C\nC>G", "G>A\nC>T"))

	out <-  ggplot(data=varprop, aes(x=Variant, y=CpGprop, fill=Variant, width=0.6)) +
		geom_bar(stat="identity")  + 
		facet_grid(Group ~ .) +
		ylab("Fraction of CpG per NpG") +
		theme_bw() + 
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			,legend.position="none"
		) + 
		scale_fill_brewer(palette = "Set1") +
		geom_text(data=varprop, 
            aes(Variant, CpGprop, label=varprop$pvalue),
            size = 4, vjust=0) + ggtitle("Fraction of guanine mutations\noccurring at CpG dinucleotides")
   
   # Create the plot
	pdfName <- paste(arguments$prefix,"context_dinucleotide.pdf", sep=".")
	pdf(pdfName, title = "Fraction of the three classes of guanine mutations occurring at CpG dinucleotides")
		print(out)
	dev.off()
	
	svgName <- paste(arguments$prefix,"context_dinucleotide.svg", sep=".")
	svg(svgName)
		print(out)
	dev.off()
	
}

if(!is.null(arguments$bySample)) {

	tmp <- snptab
	tmp$CpG[tmp$Ref=="G"] <- "NpG"
	tmp$CpG[tmp$Ref=="G" & tmp$Dinucleotide == "GC"] <- "CpG"
	
	varprop <- data.frame(table(tmp$Variant[tmp$Ref=="G"], tmp$CpG[tmp$Ref=="G"], tmp$Sample[tmp$Ref=="G"] ))
	names(varprop) <- c("Variant","Dn","Sample", "Freq")
	for(i in levels(varprop$Sample)) {
		if(exists("tmpCast")) { 
			x <- cast(varprop[varprop$Sample==i,], Variant~Dn, value = "Freq")
			x$Variant <- as.character(x$Variant)
			x <- rbind(x, c("Genome", gCpG, gNpG))
			x$Sample <- i
			tmpCast <- rbind(tmpCast, x)
		}
		if(!exists("tmpCast")) {
			tmpCast <- cast(varprop[varprop$Sample==i,], Variant~Dn, value = "Freq")
			tmpCast$Variant <- as.character(tmpCast$Variant)
			tmpCast <- rbind(tmpCast, c("Genome", gCpG, gNpG))
			tmpCast$Sample <- i
		}
	}
	varprop <- tmpCast
	rm(tmpCast)

	varprop$CpG <- as.numeric(varprop$CpG)
	varprop$NpG <- as.numeric(varprop$NpG)
	varprop$NpG <- varprop$CpG + varprop$NpG
	varprop$CpGprop <- varprop$CpG / varprop$NpG
	varprop$CpGprop[is.na(varprop$CpGprop)] <- 0
	
	varprop$pvalue <- NA
	for (i in 1:dim(varprop)[1]) {
		x <- varprop[i,c("NpG","CpG")]
		m <- matrix(c(gNpG, gCpG, x$NpG, x$CpG), nrow=2, byrow=T, dimnames=list(c("Genome","Query"),c("NpG","CpG")))
		varprop$pvalue[i] <- chisq.test(m)$p.value
	}
	pvalR <- round(varprop$pvalue,4)
	pvalR <- paste ("P", pvalR, sep="=")
	pvalR[varprop$pvalue < 0.0001] <- "P<0.0001"
	pvalR[varprop$pvalue == 0] <- ""
	varprop$pvalue <- pvalR
	
	varprop$Variant[varprop$Variant == "C>A"] <- "G>T\nC>A"
	varprop$Variant[varprop$Variant == "C>T"] <- "G>A\nC>T"
	varprop$Variant[varprop$Variant == "C>G"] <- "G>C\nC>G"
	varprop$Variant <- factor(varprop$Variant, levels=c("Genome", "G>T\nC>A", "G>C\nC>G", "G>A\nC>T"))

	out <-  ggplot(data=varprop, aes(x=Variant, y=CpGprop, fill=Variant, width=0.6)) +
		geom_bar(stat="identity")  + 
		facet_wrap( ~ Sample, ncol=2) +
		ylab("Fraction of CpG per NpG") +
		theme_bw() + 
		theme(
			plot.title = element_text(hjust=0)
			, axis.title.x = element_blank()
			,legend.position="none"
		) + 
		scale_fill_brewer(palette = "Set1") +
		geom_text(data=varprop, 
            aes(Variant, CpGprop, label=varprop$pvalue),
            size = 4, vjust=0) + ggtitle("Fraction of guanine mutations\noccurring at CpG dinucleotides")
   
   # Create the plot
	Height <- 1*nlevels(factor(varprop$Sample))

	pdfName <- paste(arguments$prefix,"context_dinucleotide.pdf", sep=".")
	pdf(pdfName, title = "Fraction of the three classes of guanine mutations occurring at CpG dinucleotides", height= Height)
		print(out)
	dev.off()
	
	svgName <- paste(arguments$prefix,"context_dinucleotide.svg", sep=".")
	svg(svgName)
		print(out)
	dev.off()
	
}

}
