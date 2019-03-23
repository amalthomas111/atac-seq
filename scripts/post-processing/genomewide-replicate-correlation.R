#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
# RScript to find the correlation of replicates using non overlapping windows
# Author: A.T

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0)
{
  stop("Usage: Rscript genowidereplicate-correlation.R <XMLconfig. file> ", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input XML file not found!\
Usage: Rscript plotinsertsize.R <XML file>", call.=FALSE)
}
#libraries
suppressPackageStartupMessages({
	library(BSgenome.Hsapiens.UCSC.hg19)
	library(Rsamtools)
	library(GenomicAlignments)
	library(gplots)
	library(RColorBrewer)
	library(XML)
})

#parse from XML
data = xmlToList(xmlParse(args[1]))
bamfiles = vector()
bam.names = vector()
for(f in data$BamFiles){
	bamfiles = c(bamfiles,f)
}
for(name in data$BamNames){
	bam.names =c(bam.names,name)
}
#chromosomes to look into
chroms = paste0("chr",c(1:22,'X','Y'))
#windowsize & organism
windowSize =1000
organism=Hsapiens

#find windows across all chromosomes
totalwindows <- GRanges()
for(i in 1:length(chroms))
{
  wstart <- seq(1,seqlengths(organism)[[chroms[i]]] - windowSize, windowSize)
  windows <- GRanges(chroms[i], IRanges(wstart, wstart+(windowSize-1)))
  totalwindows <-append(totalwindows,windows)

}

#create GRanges object for aligned reads
aln <- lapply(bamfiles, readGAlignments)
# generate counts
counts <- lapply(aln, GenomicRanges::countOverlaps, query=totalwindows)
names(counts) <- bam.names
sapply(counts, length)

#convert to data.frame
df <- data.frame(do.call("cbind",counts))
#remove rows if any element is zero
df.filt <- df[ !rowSums(df[,colnames(df)[(1:ncol(df))]]==0)>=1, ]

#plots
#find correlation
corr.matrix <- cor(do.call("cbind",df.filt))
corr.matrix
write.table(corr.matrix, file="correlation.tsv",sep="\t")

colfunc <- colorRampPalette(c("grey98",brewer.pal(6,"Blues")))
#heatmap
png(filename = "heatmap.png",width = 480, height = 480)
heatmap.2(corr.matrix, col=colfunc(10),margins=c(15,15),symm=F,
          trace="none",distfun=function (x) as.dist(1-x),
          hclust=function(x) hclust(x,method="ward.D"))
dev.off()
png(filename = paste0(bam.names[1],"-correlation.png"),width = 480, height = 480)
plot(df.filt[,1],df.filt[,2], col=rgb(0,100,0,50,maxColorValue=255),
         pch=16, log="xy", main=bam.names[1] , xlab = bam.names[1], ylab= bam.names[2],
          xaxt ='n',yaxt = 'n')
dev.off()
