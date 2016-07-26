#RScript to find the correlation of replicates using non overlapping windows
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
library(GenomicAlignments)
#input bam files
bamfiles <- c('68H.bam','68HR.bam',
             '68N3.bam','68N3R.bam',
             '68N8.bam','68N8R.bam')
#bam file names
bam.names <- c('68H','68HR',
	    '68N3','68N3R',
            '68N8','68N8R')
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
library(gplots)
library(RColorBrewer)
#find correlation
corr.matrix <- cor(do.call("cbind",df.filt))
corr.matrix
write.table(corr.matrix, file="correlation.tsv",sep="\t")

colfunc <- colorRampPalette(c("grey98",brewer.pal(6,"Blues")))
#heatmap
png(filename = "heatmap.png",width = 480, height = 480)
heatmap.2(corr.matrix, col=colfunc(10),margins=c(15,15),symm=F,
          trace="none",distfun=function (x) as.dist(1-x),hclust=function(x) hclust(x,method="ward.D"))
dev.off()
png(filename = "68H-correlation.png",width = 480, height = 480)
plot(df.filt[,1],df.filt[,2], col=rgb(0,100,0,50,maxColorValue=255), 
         pch=16, log="xy", main="68H" , xlab = "68H.R1", ylab= "68H.R2",
          xaxt ='n',yaxt = 'n')
dev.off()

png(filename = "68N3-correlation.png",width = 480, height = 480)
plot(df.filt[,3],df.filt[,4], col=rgb(0,100,0,50,maxColorValue=255), 
         pch=16, log="xy", main="68N3" ,xlab = "68N3.R1", ylab= "68N3.R2",
          xaxt ='n',yaxt = 'n')
dev.off()

png(filename = "68N8-correlation.png",width = 480, height = 480)
plot(df.filt[,5],df.filt[,6], col=rgb(0,100,0,50,maxColorValue=255), 
         pch=16, log="xy", main="68N8" ,xlab = "68N8.R1", ylab= "68N8.R2",
          xaxt ='n',yaxt = 'n')
dev.off()
#plot(log10(df.filt[,1]),log10(df.filt[,2]), col=rgb(0,100,0,50,maxColorValue=255), 
#         pch=16, main="68H" ,xlab = "log(no. of cuts) 68H",
#          xaxt ='n',yaxt = 'n')#labels=FALSE)
#dev.off()
