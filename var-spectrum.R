#!/opt/well/R/R-3.0.2/bin/Rscript


# Load libraries
{
library(RColorBrewer, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(reshape, warn.conflicts = F, quietly = T)
library(fmsb, warn.conflicts = F, quietly = T)
library(gplots, warn.conflicts = F, quietly = T)
library(fdrtool, warn.conflicts = F, quietly = T)
library(pheatmap, warn.conflicts = F, quietly = T)
library(plyr, warn.conflicts = F, quietly = T)
}


# Load functions
{
# Multiplot
	multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	  require(grid)

	  # Make a list from the ... arguments and plotlist
	  plots <- c(list(...), plotlist)

	  numPlots = length(plots)

	  # If layout is NULL, then use 'cols' to determine layout
	  if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
						ncol = cols, nrow = ceiling(numPlots/cols))
	  }

	 if (numPlots==1) {
		print(plots[[1]])

	  } else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots) {
		  # Get the i,j matrix positions of the regions that contain this subplot
		  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

		  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
										  layout.pos.col = matchidx$col))
		}
	  }
	}

# Polar histogram

#' polarHistogram builds many  histogram and arranges them around a circle to save space.
#' C. Ladroue, after Kettunen, J. et al. Genome-wide association study identifies multiple loci influencing human serum metabolite levels. Nat Genet advance online publication (2012). URL http://dx.doi.org/10.1038/ng.1073.
#' CC licence: http://creativecommons.org/licenses/by-nc-sa/3.0/ 
#' Attribution - You must attribute the work in the manner specified by the author or licensor (but not in any way that suggests that they endorse you or your use of the work).
#' Noncommercial - You may not use this work for commercial purposes.
#' Share Alike - If you alter, transform, or build upon this work, you may distribute the resulting work only under the same or similar license to this one. 
#'
#' v. 21022012
#'
#' The data frame is expected to have at least four columns: family, item, score and value. 
#' The three first columns are categorical variables, the fourth column contains non-negative values.
#' See testPolarHistogram.R for an example.
#' Each bar represents the proportion of scores for one item. Items can be grouped by families.
#' The resulting graph can be busy and might be better off saved as a pdf, with ggsave("myGraph.pdf").
#'
#' @author Christophe Ladroue
#' @param df a data frame containing the data
#' @param binSize width of the bin. Should probably be left as 1, as other parameters are relative to it.
#' @param spaceItem space between bins
#' @param spaceFamily space between families
#' @param innerRadius radius of inner circle
#' @param outerRadius radius of outer circle. Should probably be left as 1, as other parameters are relative to it.
#' @param guides a vector with percentages to use for the white guide lines
#' @param alphaStart offset from 12 o'clock in radians
#' @param circleProportion proportion of the circle to cover
#' @param direction whether the increasing count goes from or to the centre.
#' @param familyLabels logical. Whether to show family names
#' @return a ggplot object
#' @export
#' @examples
#' See testPolarHistogram.R

polarHistogram<-function(
  df,
  binSize=1,
  spaceItem=0.2,
  spaceFamily=1.2,
  innerRadius=0.3,
  outerRadius=1,
  guides=c(10,20,40,80),
  alphaStart=-0.3,
  circleProportion=0.8,
  direction="inwards",
  familyLabels=FALSE){
  
  # ordering
  df<-arrange(df,family,item,score)
  
  # summing up to one
  # TO DO: replace NA by 0 because cumsum doesn't ignore NA's.
  df<-ddply(df,.(family,item),transform,value=cumsum(value/(sum(value))))
  
  # getting previous value
  df<-ddply(df,.(family,item),transform,previous=c(0,head(value,length(value)-1)))
  
  # family and item indices. There must be a better way to do this
  df2<-ddply(df,.(family,item),summarise,indexItem=1)
  df2$indexItem<-cumsum(df2$indexItem)
  df3<-ddply(df,.(family),summarise,indexFamily=1)
  df3$indexFamily<-cumsum(df3$indexFamily)
    
  df<-merge(df,df2,by=c("family",'item'))
  df<-merge(df,df3,by="family")
    
  df<-arrange(df,family,item,score)
    
  # define the bins
  # linear projection  
  affine<-switch(direction,
                 'inwards'= function(y) (outerRadius-innerRadius)*y+innerRadius,
                 'outwards'=function(y) (outerRadius-innerRadius)*(1-y)+innerRadius,
                 stop(paste("Unknown direction")))
  
  df<-within(df,{
             xmin<-(indexItem-1)*binSize+(indexItem-1)*spaceItem+(indexFamily-1)*(spaceFamily-spaceItem)
             xmax<-xmin+binSize
             ymin<-affine(1-previous)
             ymax<-affine(1-value)
             }
             )
  
  # build the guides
  guidesDF<-data.frame(
    xmin=rep(df$xmin,length(guides)),
    y=rep(1-guides/100,1,each=nrow(df)))
  
  guidesDF<-within(guidesDF,{
    xend<-xmin+binSize
    y<-affine(y)
  })
      
  
  # Building the ggplot object
  
  totalLength<-tail(df$xmin+binSize+spaceFamily,1)/circleProportion-0

  # histograms
  p<-ggplot(df)+geom_rect(
    aes(
      xmin=xmin,
      xmax=xmax,
      ymin=ymin,
      ymax=ymax,
      fill=score)
    )
  
  # item labels
  readableAngle<-function(x){
    angle<-x*(-360/totalLength)-alphaStart*180/pi+90
    angle+ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,180,0)
  }
  readableJustification<-function(x){
    angle<-x*(-360/totalLength)-alphaStart*180/pi+90
    ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,1,0)
  }
  
  dfItemLabels<-ddply(df,.(item),summarize,xmin=xmin[1])
  dfItemLabels<-within(dfItemLabels,{
    x<-xmin+binSize/2
    angle<-readableAngle(xmin+binSize/2)
    hjust<-readableJustification(xmin+binSize/2)
    })

  p<-p+geom_text(
    aes(
      x=x,
      label=item,
      angle=angle,
      hjust=hjust),
    y=1.02,
    size=3,
    vjust=0.5,
    data=dfItemLabels)
  
  # guides  
  p<-p+geom_segment(
    aes(
      x=xmin,
      xend=xend,
      y=y,
      yend=y),
    colour="white",
    data=guidesDF)
  
  # label for guides
  guideLabels<-data.frame(
    x=0,
    y=affine(1-guides/100),
    label=paste(guides,"% ",sep='')
    )
  
  p<-p+geom_text(
    aes(x=x,y=y,label=label),
    data=guideLabels,
    angle=-alphaStart*180/pi,
    hjust=1,
    size=4)

  # family labels
  if(familyLabels){
#     familyLabelsDF<-ddply(df,.(family),summarise,x=mean(xmin+binSize),angle=mean(xmin+binSize)*(-360/totalLength)-alphaStart*180/pi)
    familyLabelsDF<-aggregate(xmin~family,data=df,FUN=function(s) mean(s+binSize))
    familyLabelsDF<-within(familyLabelsDF,{
      x<-xmin
      angle<-xmin*(-360/totalLength)-alphaStart*180/pi
    })

    p<-p+geom_text(
      aes(
        x=x,
        label=family,
        angle=angle),
    data=familyLabelsDF,
    y=1.2)
  }  
#   # empty background and remove guide lines, ticks and labels
  p<-p+theme(
    panel.background=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    legend.title=element_blank()
    )

  # x and y limits
  p<-p+xlim(0,tail(df$xmin+binSize+spaceFamily,1)/circleProportion)
  p<-p+ylim(0,outerRadius+0.2)

  # project to polar coordinates
  p<-p+coord_polar(start=alphaStart)

  # nice colour scale
  p<-p+scale_fill_brewer(palette='Spectral',type='div')

  p
}

}


# Usage
{
usage <- "
	var-spectrum.R : R script to create mutation spectrum plots from tsv
	usage : var-spectrum.R [options] <input file>
	example : var-spectrum.R prefix=outfile bySample file.tsv
	
	Input format: Tab separated file with the following names (the order is not important). Additional columns between brackets
		1. Chromosome
		2. Start
		3. End
		4. Ref
		5. Alt
		[6. Sample] : only necessary if bySample option is selected
		[7. Group] : Group name for sample
		
	Options:
		help : print usage
		prefix : prefix for output files. Default = vardesc
		bySample : plot by sample
		
	"
}


# Getting Arguments indicated by '=' and printing usage
{

arguments <- list (file=NULL, bySample=NULL, prefix="vardesc")	
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


}


# Reading files
{
#	arguments <- list (file=NULL, bySample=NULL, prefix="test")
#	arguments <- list (file=NULL, bySample=NA, prefix="test")
#	arguments$file<-"test.SNV.PRCC.tsv"
#	arguments$file<-"TCGA_CRCexome_var-spectrum.tsv"
#	arguments$file<-"TCGA_MELANOMA.snp.tsv"
# arguments <- list (file="/media/DATA/WTCHG/projects/PRCC/PRCC_vardesc_FCastro/PRCC.vardesc.tsv", bySample=T, prefix="test")
  
vartab <- read.table(arguments$file, stringsAsFactors = F, header = T)


#	vartab$Group[vartab$Sample %in% c("P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11")] <- "Type1"
#	vartab$Group[vartab$Sample %in% c("P12","P13","P14","P15","P16","P17","P18","P19","P20","P21","P22","P23")] <- "Type2"
}


# Quality control of the table and Modification
{

	# Delete space
	for (i in names(vartab)) vartab[,i] <- gsub(" ", "", vartab[,i])

	# Select only SNPs
	snptab <- vartab
	snptab <- snptab[nchar(snptab$Ref)==1 & snptab$Ref!="-", ]
	snptab <- snptab[nchar(snptab$Alt)==1 & snptab$Alt!="-", ]

	# Create Variant Variable
	snptab$Variant <- paste(snptab$Ref, snptab$Alt, sep =">")
	
	# Create Grouped Variant Variable
	snptab$gVariant[snptab$Variant == "C>A" | snptab$Variant == "G>T"] <- "C>A"
	snptab$gVariant[snptab$Variant == "C>G" | snptab$Variant == "G>C"] <- "C>G"
	snptab$gVariant[snptab$Variant == "C>T" | snptab$Variant == "G>A"] <- "C>T"
	snptab$gVariant[snptab$Variant == "T>A" | snptab$Variant == "A>T"] <- "T>A"
	snptab$gVariant[snptab$Variant == "T>C" | snptab$Variant == "A>G"] <- "T>C"
	snptab$gVariant[snptab$Variant == "T>G" | snptab$Variant == "A>C"] <- "T>G"
	
#	# Add Category All into Sample column is bySample option is specified
#	if(!is.null(arguments$bySample)) {
#		alltab <- snptab
#		alltab$Sample <- "All"
#		snptab <- rbind(snptab,alltab)
#	}
}


# PLOT: Barplot number of mutation by sample : Only if bySample != NULL
{
if(!is.null(arguments$bySample)) {
	mutCounts <- table(snptab$Sample)
	mutCounts <- data.frame(mutCounts[names(mutCounts) != "All"])
	mutCounts$Sample <- row.names(mutCounts)
	names(mutCounts) <- c("Counts","Sample")
	mutCounts$Sample <- factor(mutCounts$Sample, levels = mutCounts$Sample[order(mutCounts$Counts)])
	out <- ggplot(data=mutCounts, aes(x=Sample, y=Counts, fill=Sample)) + geom_bar(stat="identity") + theme(legend.position="none") +
		theme(
			plot.background = element_blank()
		   ,panel.grid.major = element_blank()
		   ,panel.grid.minor = element_blank()
		   ,panel.border = element_blank()
		   ,panel.background = element_blank()
		  )  + ggtitle("Mutation counts") + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
	
	pdfName <- paste(arguments$prefix,"mutCounts.pdf", sep=".")
	pdf(file = pdfName, width= 15, height= 9, title = "Mutation Counts")
		print(out)
	dev.off()
	svgName <- paste(arguments$prefix,"mutCounts.svg", sep=".")
	svg(file = svgName, width= 15, height= 9)
		print(out)
	dev.off()
}

if(!is.null(vartab$Group) && !is.null(arguments$bySample)) {
	mutCounts <- table(snptab$Sample)
	mutCounts <- data.frame(mutCounts[names(mutCounts) != "All"])
	mutCounts$Sample <- row.names(mutCounts)
	names(mutCounts) <- c("Counts","Sample")
	groupUnique <- unique(snptab[,c("Sample","Group")])
	mutCounts <- merge(mutCounts, groupUnique)
	mutCounts$Sample <- factor(mutCounts$Sample, levels = mutCounts$Sample[order(mutCounts$Counts)])
	out <- ggplot(data=mutCounts, aes(x=Sample, y=Counts, fill=Group)) + geom_bar(stat="identity") +
		theme(
			plot.background = element_blank()
		   ,panel.grid.major = element_blank()
		   ,panel.grid.minor = element_blank()
		   ,panel.border = element_blank()
		   ,panel.background = element_blank()
		  )  + ggtitle("Mutation counts") + theme(axis.text.x  = element_text(angle=90, vjust=0.5)) + 
		  scale_fill_brewer(type="diverging", palette="RdYlBu")
	
	pdfName <- paste(arguments$prefix,"mutCounts_group.pdf", sep=".")
	pdf(file = pdfName, width= 15, height= 9, title = "Mutation Counts")
		print(out)
	dev.off()
	svgName <- paste(arguments$prefix,"mutCounts_group.svg", sep=".")
	svg(file = svgName, width= 15, height= 9)
		print(out)
	dev.off()
}

}


# Create tables for mutation spectrum
{
if(is.null(arguments$bySample)) {
	stab <- data.frame(prop.table(table(snptab$Variant)))
	names(stab) <- c("Variant", "Freq")
	nsamples <- 1
	
	# Mutations grouped by type
	gtab <- data.frame(prop.table(table(snptab$gVariant)))
	names(gtab) <- c("Variant", "Freq")
	nsamples <- 1
}

if(!is.null(arguments$bySample)) {
	stab <- data.frame(prop.table(table(snptab$Sample, snptab$Variant), margin=1))
	names(stab) <- c("Sample", "Variant", "Freq")
	nsamples <- length(levels(stab$Sample))
	
	# Mutations grouped by type
	gtab <- data.frame(prop.table(table(snptab$Sample, snptab$gVariant), margin=1))
	names(gtab) <- c("Sample", "Variant", "Freq")
	nsamples <- length(levels(gtab$Sample))
}

if(!is.null(vartab$Group)) {
	sGtab <- data.frame(prop.table(table(snptab$Group, snptab$Variant), margin=1))
	names(sGtab) <- c("Group", "Variant", "Freq")
	ngroups <- length(levels(sGtab$Group))
	
	# Mutations grouped by type
	gGtab <- data.frame(prop.table(table(snptab$Group, snptab$gVariant), margin=1))
	names(gGtab) <- c("Group", "Variant", "Freq")
	ngroups <- length(levels(gGtab$Group))
}

}


# TABLE : Statistical differences in mutation spectrum by sample vs background
{
if(!is.null(arguments$bySample))  
{
	tSamVar <- table(snptab$Sample, snptab$gVariant)
	background <- matrix(colSums(tSamVar), nrow=1, dimnames = list(c("Background"), names(colSums(tSamVar))))
	
	pSamVar <- prop.table(tSamVar,1)
	weighted <-	apply(pSamVar, 2, mean) * sum(background)
	weighted <- matrix(weighted, nrow=1, dimnames = list(c("Weighted"), names(colSums(pSamVar))))
	
	out.tab <- matrix(nrow=dim(tSamVar)[1], ncol=dim(tSamVar)[2]+3, dimnames = list(c(row.names(tSamVar)), c(names(colSums(tSamVar)), "counts", "p.value", "p.value.weighted")))
	for ( i in 1:dim(tSamVar)[1]) {
		x <- rbind(background, tSamVar[i,])
		p.val <- chisq.test(x)$p.value
		x <- rbind(weighted, tSamVar[i,])
		p.val.weighted <- chisq.test(x)$p.value
		out.tab[i,] <- c(round(prop.table(tSamVar[i,]), 3), sum(tSamVar[i,]), p.val, p.val.weighted)
	}
	out.tab <- cbind(row.names(out.tab), out.tab)
	colnames(out.tab)[1] <- "Sample"
		
	out.tab <- cbind(out.tab, fdrtool(as.numeric(out.tab[,"p.value"]), statistic="pvalue", plot=F)$qval)
	out.tab <- cbind(out.tab, fdrtool(as.numeric(out.tab[,"p.value.weighted"]), statistic="pvalue", plot=F)$qval)
	dimnames(out.tab)[[2]][c(dim(out.tab)[2]-1,dim(out.tab)[2])] <- c("q.value", "q.value.weighted")

	dStatsSpectrum <- data.frame(out.tab)
	write.table(out.tab, file = paste(arguments$prefix,"spectrum_table.txt", sep="."), 
		quote = F, sep = "\t", row.names = F)
		
	dStatsSpectrum$p.value <- as.numeric(as.character(dStatsSpectrum$p.value))
	dStatsSpectrum$p.value <- -log10(dStatsSpectrum$p.value)
	dStatsSpectrum <- dStatsSpectrum[order(dStatsSpectrum$p.value),]
	dStatsSpectrum$Sample <- factor(dStatsSpectrum$Sample, levels=dStatsSpectrum$Sample)
	
	dStatsSpectrum$p.value.weighted <- as.numeric(as.character(dStatsSpectrum$p.value.weighted))
	dStatsSpectrum$p.value.weighted <- -log10(dStatsSpectrum$p.value.weighted)
	
	dStatsSpectrum$q.value <- as.numeric(as.character(dStatsSpectrum$q.value))
	dStatsSpectrum$q.value <- -log10(dStatsSpectrum$q.value)
	
	dStatsSpectrum$q.value.weighted <- as.numeric(as.character(dStatsSpectrum$q.value.weighted))
	dStatsSpectrum$q.value.weighted <- -log10(dStatsSpectrum$q.value.weighted)
}
}


# PLOT : Barplot mutation spectrum using Non-Grouped mutations
{
Height <- 9
if (nsamples != 1) Height <- nsamples * 1.6

if(is.null(arguments$bySample)) {
	out <- ggplot(data=stab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + theme(legend.position="none") +
		theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	   ) + ggtitle("Mutation spectrum")
  }

if(!is.null(arguments$bySample) ) {
	out <- ggplot(data=stab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + facet_grid(Sample ~ .) + 
		theme_bw() +
		theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,legend.position="none"
		) + ggtitle("Mutation spectrum")
  }

pdfName <- paste(arguments$prefix,"spectrum_nongrouped.pdf", sep=".")
pdf(pdfName, width= 7, height= Height, title = "Mutation Spectrum")
	print(out)
dev.off()	
svgName <- paste(arguments$prefix,"spectrum_nongrouped.svg", sep=".")
svg(svgName, width= 7, height= Height)
	print(out)
dev.off()	

if(!is.null(vartab$Group)) {
	if (ngroups != 1) Height <- ngroups * 1.6
	out <- ggplot(data=sGtab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + facet_grid(Group ~ .) + 
		theme_bw() +
		theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,legend.position="none"
		) + ggtitle("Mutation spectrum")
		
pdfName <- paste(arguments$prefix,"spectrum_nongrouped_group.pdf", sep=".")
pdf(pdfName, width= 7, height= Height, title = "Mutation Spectrum")
	print(out)
dev.off()	
svgName <- paste(arguments$prefix,"spectrum_nongrouped_group.svg", sep=".")
svg(svgName, width= 7, height= Height)
	print(out)
dev.off()	

  }
  
  
}


# PLOT : Barplot mutation spectrum Grouped mutations
{
Height <- 9
if (nsamples != 1) Height <- nsamples * 1.6

if(is.null(arguments$bySample)) { 
	out <- ggplot(data=gtab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + theme(legend.position="none") +
		theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	  )  + ggtitle("Mutation spectrum")
  }
if(!is.null(arguments$bySample)) { 
	out <- ggplot(data=gtab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + facet_grid(Sample ~ .) +
		theme_bw() +
		theme(
		plot.background = element_blank()
		,panel.grid.major = element_blank()
		,panel.grid.minor = element_blank()
		,legend.position="none"
	  ) + ggtitle("Mutation spectrum")
	  
  }
  
pdfName <- paste(arguments$prefix,"spectrum.pdf", sep=".")
pdf(pdfName, width= 7, height= Height, title = "Mutation Spectrum")
print(out)
dev.off()

svgName <- paste(arguments$prefix,"spectrum.svg", sep=".")
pdf(svgName, width= 7, height= Height)
print(out)
dev.off()

if(!is.null(vartab$Group)) {
	if (ngroups != 1) Height <- ngroups * 1.6
	out <- ggplot(data=gGtab, aes(x=Variant, y=Freq, fill=Variant)) + geom_bar(stat="identity") + coord_flip() + scale_fill_brewer(palette = "Spectral") + facet_grid(Group ~ .) +
		theme_bw() +
		theme(
		plot.background = element_blank()
		,panel.grid.major = element_blank()
		,panel.grid.minor = element_blank()
		,legend.position="none"
	  ) + ggtitle("Mutation spectrum")
	  
pdfName <- paste(arguments$prefix,"spectrum_group.pdf", sep=".")
pdf(pdfName, width= 7, height= Height, title = "Mutation Spectrum")
print(out)
dev.off()

svgName <- paste(arguments$prefix,"spectrum_group.svg", sep=".")
pdf(svgName, width= 7, height= Height)
print(out)
dev.off()
  
  }
 
  
}


# PLOT : Stacked barplot mutation spectrum Grouped mutations (sorted by Significance)
{
if(!is.null(arguments$bySample))  
{
# Sort data according to P-Value : 
gtabS <- gtab
orderedSampleFactors <- dStatsSpectrum$Sample[order(dStatsSpectrum$p.value)]
gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

## Sort data according to Ts proportion
#gtabS <- gtab
#gtabS$Ts <- NA
#for (i in levels(gtabS$Sample)){
#	tsSum <- sum(gtabS[gtabS$Sample==i & (gtabS$Variant=="C>T" | gtabS$Variant=="T>C"), "Freq"])
#	gtabS[gtabS$Sample==i, "Ts"] <- tsSum
#}
#orderedSampleFactors <- unique(gtabS[order(gtabS$Ts),"Sample"])
#gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

# Build the plot

# Build the plot
xTextSize <- 1
if(nlevels(gtabS$Sample) < 260) xTextSize <- 2
if(nlevels(gtabS$Sample) < 200) xTextSize <- 3
if(nlevels(gtabS$Sample) < 160) xTextSize <- 5
if(nlevels(gtabS$Sample) < 140) xTextSize <- 6
if(nlevels(gtabS$Sample) < 100) xTextSize <- 8


out <- ggplot(data=gtabS, aes(x=Sample, y=Freq, fill=Variant)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Spectral") + theme(
    plot.background = element_blank()
   ,panel.grid.major = element_blank()
   ,panel.grid.minor = element_blank()
   ,panel.border = element_blank()
   ,panel.background = element_blank()
  )  + 
  ggtitle("Mutation spectrum (sorted by profile significance)") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=xTextSize))

stacked.barplot <- out

pdfName <- paste(arguments$prefix,"spectrum_stacked.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(out)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_stacked.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(out)
dev.off()


}

}


# PLOT : Stacked barplot mutation spectrum Grouped mutations + Significance differences in mutation spoutectrum by sample vs background
{
if(!is.null(arguments$bySample))  
{
	hlines <- data.frame(cutoff="p-val<0.05", yval=-log10(0.05))

	out <-  ggplot(data=dStatsSpectrum, aes(x=Sample, y=p.value)) + 
		geom_point() +
		scale_y_continuous(name="-log10(p-value)") + 
		geom_hline(yintercept=-log10(0.05), colour="#990000", linetype="dashed") + 
	  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=xTextSize))
		
	stacked.barplot <- stacked.barplot + theme(legend.position="bottom")

pdfName <- paste(arguments$prefix,"spectrum_stacked_significance.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	multiplot(stacked.barplot, out, cols=1)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_stacked_significance.svg", sep=".")
svg(svgName, width= 12, height= 9)
	multiplot(stacked.barplot, out, cols=1)
dev.off()

}
}


# PLOT : Stacked barplot mutation spectrum Grouped mutations (sorted by Number of mutations)
{
if(!is.null(arguments$bySample))  
{
# Sort data according number of mutations : 
gtabS <- gtab
CountsTab <- data.frame(table(snptab$Sample))
names(CountsTab) <- c("Sample","Freq")
CountsTab <- CountsTab[order(CountsTab$Freq),]
CountsTab$Sample<-factor(as.character(CountsTab$Sample), levels = as.character(CountsTab$Sample))
orderedSampleFactors <- levels(CountsTab$Sample)
gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

## Sort data according to Ts proportion
#gtabS <- gtab
#gtabS$Ts <- NA
#for (i in levels(gtabS$Sample)){
#	tsSum <- sum(gtabS[gtabS$Sample==i & (gtabS$Variant=="C>T" | gtabS$Variant=="T>C"), "Freq"])
#	gtabS[gtabS$Sample==i, "Ts"] <- tsSum
#}
#orderedSampleFactors <- unique(gtabS[order(gtabS$Ts),"Sample"])
#gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

# Build the plot
xTextSize <- 1
if(nlevels(gtabS$Sample) < 260) xTextSize <- 2
if(nlevels(gtabS$Sample) < 200) xTextSize <- 3
if(nlevels(gtabS$Sample) < 160) xTextSize <- 5
if(nlevels(gtabS$Sample) < 140) xTextSize <- 6
if(nlevels(gtabS$Sample) < 100) xTextSize <- 8

out <- ggplot(data=gtabS, aes(x=Sample, y=Freq, fill=Variant)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Spectral") + theme(
    plot.background = element_blank()
   ,panel.grid.major = element_blank()
   ,panel.grid.minor = element_blank()
   ,panel.border = element_blank()
   ,panel.background = element_blank()
  )  + 
  ggtitle("Mutation spectrum (sorted by number of variants per sample)") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=xTextSize))

stacked.barplot <- out

pdfName <- paste(arguments$prefix,"spectrum_stacked_sortCounts.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(out)
dev.off()

svgName <- paste(arguments$prefix,"spectrum_stacked_sortCounts.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(out)
dev.off()


}

}


# PLOT : Stacked barplot mutation spectrum Grouped mutations + Significance differences in mutation spoutectrum by sample vs background
{
if(!is.null(arguments$bySample))  
{

	out <-  ggplot(data=CountsTab, aes(x=Sample, y=Freq)) + 
		geom_point() +
		scale_y_continuous(name="Number of Mutations") + 
	  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=xTextSize))
		
	stacked.barplot <- stacked.barplot + theme(legend.position="bottom")

pdfName <- paste(arguments$prefix,"spectrum_stacked_counts.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	multiplot(stacked.barplot, out, cols=1)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_stacked_counts.svg", sep=".")
svg(svgName, width= 12, height= 9)
	multiplot(stacked.barplot, out, cols=1)
dev.off()

}
}



# PLOT : Polar Histogram mutation spectrum Grouped mutations (sorted by Significance)
{

if(!is.null(arguments$bySample))  
{

	# Sort data according to P-Value : 
	gtabS <- gtab
	orderedSampleFactors <- dStatsSpectrum$Sample[order(dStatsSpectrum$p.value)]
	gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

#	# Sort data according to Ts proportion
#	gtabS <- gtab
#	gtabS$Ts <- NA
#	for (i in levels(gtabS$Sample)){
#		tsSum <- sum(gtabS[gtabS$Sample==i & (gtabS$Variant=="C>T" | gtabS$Variant=="T>C"), "Freq"])
#		gtabS[gtabS$Sample==i, "Ts"] <- tsSum
#	}
#	orderedSampleFactors <- unique(gtabS[order(gtabS$Ts),"Sample"])
#	gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

	# Build the plot
	dpolar <- gtabS
	names(dpolar) <- c("item","score","value")
	dpolar$family <- 1
	p<-polarHistogram(dpolar,familyLabel=FALSE, guides=c(25,50,75), circleProportion=0.95)

pdfName <- paste(arguments$prefix,"spectrum_polar_significance.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(p)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_polar_significance.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(p)
dev.off()
}

if(!is.null(vartab$Group) & !is.null(arguments$bySample)) {
	# Build the plot
	groupUnique <- unique(snptab[,c("Sample","Group")])
	dpolar <- merge(gtab, groupUnique)
	names(dpolar) <- c("item","score","value", "family")
	p<-polarHistogram(dpolar,familyLabel=TRUE, guides=c(25,50,75), circleProportion=0.95, spaceFamily=8, spaceItem=0.001)


pdfName <- paste(arguments$prefix,"spectrum_polar__significance_group.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(p)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_polar_significance_group.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(p)
dev.off()

}

}


# PLOT : Polar Histogram mutation spectrum Grouped mutations (sorted by Number of mutations)
{

if(!is.null(arguments$bySample))  
{

	# Sort data according to Mutations Counts : 
	gtabS <- gtab
	CountsTab <- data.frame(table(snptab$Sample))
	names(CountsTab) <- c("Sample","Freq")
	CountsTab <- CountsTab[order(CountsTab$Freq),]
	CountsTab$Sample<-factor(as.character(CountsTab$Sample), levels = as.character(CountsTab$Sample))
	orderedSampleFactors <- levels(CountsTab$Sample)
	gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)
	
#	# Sort data according to Ts proportion
#	gtabS <- gtab
#	gtabS$Ts <- NA
#	for (i in levels(gtabS$Sample)){
#		tsSum <- sum(gtabS[gtabS$Sample==i & (gtabS$Variant=="C>T" | gtabS$Variant=="T>C"), "Freq"])
#		gtabS[gtabS$Sample==i, "Ts"] <- tsSum
#	}
#	orderedSampleFactors <- unique(gtabS[order(gtabS$Ts),"Sample"])
#	gtabS$Sample <- factor(as.character(gtabS$Sample), levels = orderedSampleFactors)

	# Build the plot
	dpolar <- gtabS
	names(dpolar) <- c("item","score","value")
	dpolar$family <- 1
	p<-polarHistogram(dpolar,familyLabel=FALSE, guides=c(25,50,75), circleProportion=0.95)

pdfName <- paste(arguments$prefix,"spectrum_polar_counts.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(p)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_polar_counts.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(p)
dev.off()
}

if(!is.null(vartab$Group) & !is.null(arguments$bySample)) {
	# Build the plot
	groupUnique <- unique(snptab[,c("Sample","Group")])
	dpolar <- merge(gtab, groupUnique)
	names(dpolar) <- c("item","score","value", "family")
	p<-polarHistogram(dpolar,familyLabel=TRUE, guides=c(25,50,75), circleProportion=0.95, spaceFamily=8, spaceItem=0.001)


pdfName <- paste(arguments$prefix,"spectrum_polar__counts_group.pdf", sep=".")
pdf(pdfName, width= 12, height= 9, title = "Mutation Spectrum")
	print(p)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_polar_counts_group.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(p)
dev.off()

}

}

# PLOT : Transition / Transversion rate
{

if(!is.null(arguments$bySample))  
{

tstv <- snptab
tstv$TsTv <- "Transvertion"
tstv$TsTv[tstv$Variant=="C>T"] <- "Transition"
tstv$TsTv[tstv$Variant=="T>C"] <- "Transition"
tstv$TsTv[tstv$Variant=="A>G"] <- "Transition"
tstv$TsTv[tstv$Variant=="G>A"] <- "Transition"

tstvtab <- table(tstv$Sample,tstv$TsTv)
tstvdf <- data.frame(matrix(tstvtab, ncol=2, dimnames=list(c(), dimnames(tstvtab)[[2]])))
tstvdf$Sample <- row.names(tstvtab)
tstvdf$TsTv <- tstvdf$Transition / tstvdf$Transvertion
tstvdf <- tstvdf[tstvdf$Sample != "All", ]

out <- ggplot(data=tstvdf, aes(x=Sample, y=TsTv)) +
    geom_bar(fill="red",stat="identity" ) +
    guides(fill=FALSE) + 
    theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_blank()
	   ,panel.grid.minor = element_blank()
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	   ,axis.text.x  = element_text(angle=90, vjust=0.5)
	  ) + 
	  ggtitle("Transition/Transversion Ratio")

pdfName <- paste(arguments$prefix,"TsTv.pdf", sep=".")
pdf(pdfName)
	print(out)    
dev.off()  
svgName <- paste(arguments$prefix,"TsTv.svg", sep=".")
svg(svgName, width= 12, height= 9)
	print(out)    
dev.off()  

}

}


# PLOT : RadarPlot mutation spectrum
{

if(!is.null(arguments$bySample))  
{

# GGPLOT
	out <- ggplot(gtab, aes(x = Variant, y = Freq, colour = Sample, group = Sample)) +
	  geom_line() +
	  coord_polar(theta = "x", direction = -1) +
	  theme(
		plot.background = element_blank()
	   ,panel.grid.major = element_line(colour='grey', linetype = 'dotted')
	   ,panel.border = element_blank()
	   ,panel.background = element_blank()
	  ) 
pdfName <- paste(arguments$prefix,"spectrum_spiderChart.pdf", sep=".")
pdf(pdfName, width= 9, height= 9, title = "Mutation Spectrum")
	print(out)
dev.off()
svgName <- paste(arguments$prefix,"spectrum_spiderChart.svg", sep=".")
svg(svgName, width= 9, height= 9)
	print(out)
dev.off()


  
# RADARCHART


	radardf <- cast(gtab, Sample~Variant, value = "Freq")
	row.names(radardf) <- radardf$Sample
	radardf$Sample <- NULL
	radardf <- rbind(c(0,0,0,0),radardf) # Min
	radardf <- rbind(c(1,1,1,1),radardf) # Max

	if (nsamples <= 11) colplot <- brewer.pal(nsamples,"Spectral")
	if (nsamples > 11) colplot <- brewer.pal(11,"Spectral")

pdfName <- paste(arguments$prefix,"spectrum_radarChart.pdf", sep=".")
pdf(pdfName, width= 12, height= 12, title = "Mutation Spectrum")
	radarchart(radardf, axistype=1, pty=32, plty=1, plwd=3, axislabcol="grey", 
		na.itp=FALSE, pcol=colplot, title="Mutation spectrum")
	legend("topleft", row.names(radardf)[c(-1,-2)], col = colplot, lwd = 3, bty = "n", cex = .8)
dev.off()

svgName <- paste(arguments$prefix,"spectrum_radarChart.svg", sep=".")
svg(svgName, width= 12, height= 12)

	radarchart(radardf, axistype=1, pty=32, plty=1, plwd=3, axislabcol="grey", 
		na.itp=FALSE, pcol=colplot, title="Mutation spectrum")

	legend("topleft", row.names(radardf)[c(-1,-2)], col = colplot, lwd = 3, bty = "n", cex = .8)
dev.off()

} 

}


# PLOT : HeatMap mutation spectrum
{

if(!is.null(arguments$bySample))  
{
# Use heatmap.3 if more than one group is specified. Source and example at ~/R/functions/heatmap.3.R
# Use pheatmap if group is not specified
	
# Create Data with heatmap layout
htab <- cast(gtab, Sample~Variant)
htab <- htab[htab$Sample != "All", ]
if(!is.null(vartab$Group)) {
	groupSample <- unique(vartab[,c("Sample","Group")])
	htab <- merge(htab, groupSample)
	group_colors <- as.character(as.numeric(as.factor(htab$Group)))
	group_levels <- levels(as.factor(htab$Group))
	htab$Group <- NULL
}
rownames(htab) <- htab$Sample

# Data Matrix
htab$Sample <- NULL
htabm <- t(data.matrix(htab))
def.palette <- palette()

# Color scheme
heatCol <- brewer.pal(5, "PuBu") # palette for heatmap

# Plot size
Height <- 0.05*dim(htabm)[2]
Width <- 0.08*dim(htabm)[2]

# Column text size
xTextSize <- 0.5
if(dim(htabm)[2] < 300) xTextSize <- 1
if(dim(htabm)[2] < 260) xTextSize <- 2
if(dim(htabm)[2] < 200) xTextSize <- 3
if(dim(htabm)[2] < 160) xTextSize <- 4
if(dim(htabm)[2] < 140) xTextSize <- 5


# Plotting data
pdfName <- paste(arguments$prefix,"spectrum_heatmap.pdf", sep=".")
pdf(pdfName, width= Width, height= Height, title = "Mutation Spectrum")
	palette(brewer.pal(8,"Set1")) # palette for grouping scheme
	if(!is.null(vartab$Group)) {
		heatmap.2(t(htabm), col=heatCol, dendrogram="row", scale="row", RowSideColors=group_colors, ylab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.01, 0.01), sepcolor='white', keysize=1)
		
		legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
	}

#	if(is.null(vartab$Group)) heatmap.2(t(htabm), col=heatCol, dendrogram="row", scale="row", ylab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.01, 0.01), sepcolor='white', keysize=1)

	if(is.null(vartab$Group)) pheatmap(htabm, border_color=NA, col=heatCol, fontsize_col = xTextSize, fontsize=xTextSize)

		
	palette(def.palette)

dev.off()

svgName <- paste(arguments$prefix,"spectrum_heatmap.svg", sep=".")
svg(svgName, width= Width, height= Height)
	palette(brewer.pal(8,"Set1")) # palette for grouping scheme
	if(!is.null(vartab$Group)) {
		heatmap.2(t(htabm), col=heatCol, dendrogram="row", scale="row", RowSideColors=group_colors, ylab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.01, 0.01), sepcolor='white', keysize=1)
		
		legend("topright", group_levels, fill = as.character(as.numeric(as.factor(group_levels))), bty = 'n')
	}

#	if(is.null(vartab$Group)) heatmap.2(t(htabm), col=heatCol, dendrogram="row", scale="row", ylab="Samples", na.rm=T, trace="none", colsep=c(1:dim(htabm)[2]), rowsep=c(1:dim(htabm)[1]), sepwidth=c(0.01, 0.01), sepcolor='white', keysize=1)

	if(is.null(vartab$Group)) pheatmap(htabm, border_color=NA, col=heatCol, fontsize_col = xTextSize)
	palette(def.palette)

dev.off()


} 

}



