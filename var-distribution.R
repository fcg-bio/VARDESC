#!/opt/well/R/R-3.0.2/bin/Rscript

# Load libraries
{
library(RColorBrewer, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(reshape, warn.conflicts = F, quietly = T)
library(fmsb, warn.conflicts = F, quietly = T)
library(lattice, warn.conflicts = F, quietly = T)
}


# Load functions
{

#http://onertipaday.blogspot.co.uk/2007/09/plotting-two-or-more-overlapping.html
plot.multi.dens <- function(s, xlab = xlab, main = main)
{
	junk.x = NULL
	junk.y = NULL
	for(i in 1:length(s))
	{
		junk.x = c(junk.x, density(s[[i]])$x)
		junk.y = c(junk.y, density(s[[i]])$y)
	}
	xr <- range(junk.x)
	yr <- range(junk.y)
	plot(density(s[[1]]), xlim = xr, ylim = yr, main = main, xlab = xlab)
	for(i in 1:length(s))
	{
		lines(density(s[[i]]), xlim = xr, ylim = yr, col = i)
	}
}

}


# Usage
{
usage <- "
	var-distribution.R : R script to create mutation distribution across CHR plots from tsv
	usage : var-distribution.R [options] <input file>
	example : var-distribution.R chrCoord=file centromere=file prefix=outfile bySample file.tsv
	
	Input format: Tab separated file with the following names (the order is not important)
		1. Chromosome
		2. Start
		3. End
		4. Ref
		5. Alt
		
	Options:
		help : print usage
		prefix : prefix for output files. Default = vardesc
		bySample : plot by sample 
		centromere : provide bed file with centromere coordinates tabular separated. Example:
			chr1\t121535434\t124535434
		chrCoord : provide bed file with chromosomes coordinates tabular separated. Example:
			chr1\t1\t249250621
		
	"
}


# Getting Arguments indicated by '=' and printing usage
{

arguments <- list (file=NULL, bySample=NULL, prefix="vardesc", centromere=NULL, chrCoord=NULL)	
argsRaw<-commandArgs(trailingOnly = T)

options(warning.length=8000)
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

if (is.null(arguments$centromere)) stop(usage)
if (is.null(arguments$chrCoord)) stop(usage)

if(!file.exists(arguments$centromere)) stop("Centromere coordinates file doesn't exists",usage)
if(!file.exists(arguments$chrCoord)) stop("Chromosomes coordinates file doesn't exists",usage)


}


# Reading files
{
#	arguments <- list (file=NULL, bySample=NULL, prefix="test", centromere="hg19_noChr.centromere.bed", chrCoord="hg19_noChr.coordinates.bed")
#	arguments <- list (file=NULL, bySample=NA, prefix="test", centromere="hg19_noChr.centromere.bed", chrCoord="hg19_noChr.coordinates.bed")
#	arguments$file<-"test.SNV.PRCC.tsv"

vartab <- read.delim(arguments$file, stringsAsFactors = F)


# Reading chromosomes coordinates
chrcoord <- read.delim(arguments$chrCoord, stringsAsFactors = F, col.names=c("chr","start","end"), header=F)

# Reading centromere position file
cm <- read.delim(arguments$centromere, stringsAsFactors = F, col.names=c("chr","start","end"), header=F)

}


# Quality control of the table
{

	# Delete space
	for (i in names(vartab)) vartab[,i] <- gsub(" ", "", vartab[,i])
	
	# Add Category All into Sample column is bySample option is specified
	if(!is.null(arguments$bySample)) {
		alltab <- vartab
		alltab$Sample <- "All"
		vartab <- rbind(vartab,alltab)
	}
}


# Annotate distance from Centromere
{
	cm.distance <- function(x,y){
		res <- NA
		if(x %in% cm$chr) {
		s <- as.numeric(cm[cm$chr==x,"start"])
		e <- as.numeric(cm[cm$chr==x,"end"])
		y <- as.numeric(y)
		res <- s-y
		if(abs(s-y) > abs(e-y)) {
		res <- e-y
		}
		}
		res
	}

	vartab$cmdist <- unlist(apply(vartab,1, function(x) cm.distance(x["Chromosome"],x["Start"])))
}


# Annotate relative distance to Centromere
{

	chrlength <- as.list(chrcoord$end)
	names(chrlength) <- chrcoord$chr

	cm.relative <- function(chr,y){
		res <- NA
		if(chr %in% cm$chr) {
		chrl <- chrlength[[chr]]
		s <- as.numeric(cm[cm$chr==chr,"start"])
		e <- as.numeric(cm[cm$chr==chr,"end"])
		y <- as.numeric(y)
		
		if(abs(s-y) < abs(e-y)) {
			res <- (y / s) -1
		}
		if(abs(s-y) > abs(e-y)) {
			res <- (y-e) / (chrl-e)
		}
		}
		res
	}

	vartab$cmrel <- unlist(apply(vartab,1, function(x) cm.relative(x["Chromosome"],x["Start"])))

}


# PLOT: Distribution of variants position relative to centromere position
{

if(is.null(arguments$bySample))  
{
out <- ggplot(vartab, aes(x=cmrel)) + 
	geom_density (alpha=.4, na.rm = T ) +
	xlab("Chromosome position relative to centromere") +
	theme(
	plot.background = element_blank()
	,panel.grid.major = element_blank()
	,panel.grid.minor = element_blank()
	,panel.border = element_blank()
	,panel.background = element_blank()
	,axis.line = element_line(colour='grey')
	)	
}

if(!is.null(arguments$bySample))  
{
out <- ggplot(vartab, aes(x=cmrel, colour=Sample)) + 
	geom_density (alpha=.4, na.rm = T, show_guide=FALSE) +
	stat_density(aes(x=cmrel, colour=Sample), geom="line",position="identity") +
	xlab("Chromosome position relative to centromere") +
	theme(
	plot.background = element_blank()
	,panel.grid.major = element_blank()
	,panel.grid.minor = element_blank()
	,panel.border = element_blank()
	,panel.background = element_blank()
	,axis.line = element_line(colour='grey')
	)
}

pdfName <- paste(arguments$prefix,"centromere_dist.pdf", sep=".")
pdf(pdfName, width= 15, height= 6, title = "Chromosome position relative to centromere")
	print(out)
dev.off()
svgName <- paste(arguments$prefix,"centromere_dist.svg", sep=".")
pdf(svgName, width= 15, height= 6)
	print(out)
dev.off()

}


# PLOT: Distribution of variants across Chromosome - position relative to centromere position
{


if(is.null(arguments$bySample))  
{

pdfName <- paste(arguments$prefix,"chromosome_pos_dist.pdf", sep=".")
pdf(pdfName, width= 19, height= 14, title = "Distribution across chromosome")
	nChr <- length(levels(as.factor(vartab$Chromosome)))
	par(mfrow=c(ceiling(nChr/5),5))
	for (chromosome in levels(as.factor(vartab$Chromosome))){
		d <- density(vartab$cmrel[vartab$Chromosome == chromosome])
		plot(d, main = chromosome, xlab = "Chromosome position")
		polygon(d, col="red", border="red") 
	}
dev.off()

svgName <- paste(arguments$prefix,"chromosome_pos_dist.pdf", sep=".")
svg(svgName, width= 19, height= 14)
	nChr <- length(levels(as.factor(vartab$Chromosome)))
	par(mfrow=c(ceiling(nChr/5),5))
	for (chromosome in levels(as.factor(vartab$Chromosome))){
		d <- density(vartab$cmrel[vartab$Chromosome == chromosome])
		plot(d, main = chromosome, xlab = "Chromosome position")
		polygon(d, col="red", border="red") 
	}
dev.off()

}



if(!is.null(arguments$bySample))  
{
nChr <- nlevels(as.factor(vartab$Chromosome))

pdfName <- paste(arguments$prefix,"chromosome_pos_dist.pdf", sep=".")
pdf(pdfName, width= 19, height= 14, title = "Distribution across chromosome")
	par(mfrow=c(ceiling(nChr/5),5))
	for (chromosome in levels(as.factor(vartab$Chromosome))){
		subdata <- vartab[vartab$Chromosome==chromosome,]
		list_cmrel <- list()
		for (Sample in levels(as.factor(subdata$Sample))) {
			sub_cmrel <- subdata[subdata$Sample==Sample, "cmrel"]
			if(length(sub_cmrel) < 2 ) sub_cmrel <- c(sub_cmrel, 0)
			list_cmrel <- append(list_cmrel, list(sub_cmrel))
		}
		
		plot.multi.dens(list_cmrel, xlab = "Relative chromosome position", main=chromosome)
		if(chromosome == levels(as.factor(vartab$Chromosome))[1]) legend("topleft", levels(as.factor(vartab$Sample)), fill=2+(0:nlevels(as.factor(vartab$Sample))), bty = "n", cex = .5)
	}
dev.off()

svgName <- paste(arguments$prefix,"chromosome_pos_dist.svg", sep=".")
pdf(svgName, width= 19, height= 14)
	par(mfrow=c(ceiling(nChr/5),5))
	for (chromosome in levels(as.factor(vartab$Chromosome))){
		subdata <- vartab[vartab$Chromosome==chromosome,]
		list_cmrel <- list()
		for (Sample in levels(as.factor(subdata$Sample))) {
			sub_cmrel <- subdata[subdata$Sample==Sample, "cmrel"]
			if(length(sub_cmrel) < 2 ) sub_cmrel <- c(sub_cmrel, 0)
			list_cmrel <- append(list_cmrel, list(sub_cmrel))
		}
		
		plot.multi.dens(list_cmrel, xlab = "Relative chromosome position", main=chromosome)
		if(chromosome == levels(as.factor(vartab$Chromosome))[1]) legend("topleft", levels(as.factor(vartab$Sample)), fill=2+(0:nlevels(as.factor(vartab$Sample))), bty = "n", cex = .5)
	}
dev.off()



}


}


# PLOT: Distribution of variants by Chromosome
{

if(is.null(arguments$bySample))  
{

chrDist <- data.frame(table(vartab$Chromosome))
names(chrDist) <- c("Chromosome", "Counts")

out <- ggplot(data=chrDist, aes(x=Chromosome, y=Counts)) +
    geom_bar(colour="black", stat="identity") +
    guides(fill=FALSE) + 
    theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	  ) + ggtitle("Distribution of mutations by chromosome")
    

}

if(!is.null(arguments$bySample))  
{

chrDist <- data.frame(table(vartab$Chromosome, vartab$Sample))
names(chrDist) <- c("Chromosome", "Sample", "Counts")
chrDist <- chrDist[chrDist$Sample != "All", ]

out<-ggplot(data=chrDist, aes(x=Sample, y=Counts, fill=Sample)) +
    geom_bar(colour="black", stat="identity") +
    facet_wrap( ~ Chromosome, ncol=4) +
    guides(fill=FALSE) + 
    theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	   ,axis.text.x  = element_text(angle=90, vjust=0.5)
	  ) + ggtitle("Distribution of mutations by chromosome")
    

}

pdfName <- paste(arguments$prefix,"chromosome_dist.pdf", sep=".")
pdf(pdfName, width= 19, height= 14, title = "Distribution by chromosome")
	print(out)
dev.off()

svgName <- paste(arguments$prefix,"chromosome_dist.svg", sep=".")
svg(svgName, width= 19, height= 14)
	print(out)
dev.off()

}


# PLOT: Rate of variants by Chromosome
{

if(is.null(arguments$bySample))  
{

chrDist <- data.frame(table(vartab$Chromosome))
names(chrDist) <- c("Chromosome", "Counts")
chrDist <- merge(chrDist, chrcoord, by.x="Chromosome", by.y="chr")
chrDist$Rate <- 1000000*(chrDist$Counts / chrDist$end)

out <- ggplot(data=chrDist, aes(x=Chromosome, y=Rate)) +
    geom_bar(colour="black", stat="identity") +
    guides(fill=FALSE) + 
    theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	  ) + 
	  ylab("Number of mutations per million bases") + 
	  ggtitle("Distribution of mutations by chromosome")
    

}

if(!is.null(arguments$bySample))  
{

chrDist <- data.frame(table(vartab$Chromosome, vartab$Sample))
names(chrDist) <- c("Chromosome", "Sample", "Counts")
chrDist <- chrDist[chrDist$Sample != "All", ]
chrDist <- merge(chrDist, chrcoord, by.x="Chromosome", by.y="chr")
chrDist$Rate <- 1000000*(chrDist$Counts / chrDist$end)


out <- ggplot(data=chrDist, aes(x=Sample, y=Counts, fill=Sample)) +
    geom_bar(colour="black", stat="identity") +
    facet_wrap( ~ Chromosome, ncol=4) +
    guides(fill=FALSE) + 
    theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	   ,axis.text.x  = element_text(angle=90, vjust=0.5)
	  ) + ggtitle("Distribution of mutations by chromosome")
    

}


pdfName <- paste(arguments$prefix,"chromosome_rate.pdf", sep=".")
pdf(pdfName, width= 19, height= 14, title = "Rate of variants by chromosome")
	print(out)
dev.off()
svgName <- paste(arguments$prefix,"chromosome_rate.scg", sep=".")
pdf(svgName, width= 19, height= 14)
	print(out)
dev.off()

}






