#  !/usr/bin/env Rscript
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
#  Rscript to identify the differential accessible site between two atac-seq samples.
#  This script is based on a Bioconductor workflow created by Aaron Lun and Gordon Smyth.
#  The code was adapted for our ATAC-seq data.
#  For more details, check https://www.bioconductor.org/help/workflows/chipseqDB
#  Paper link: http://nar.oxfordjournals.org/content/42/11/e95

#library(session)
#restore.session("RSession.Rda")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0)
{
  stop("Usage: Rscript differential_atac.R <XMLconfig.file> ", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input XML file not found!\
Usage: Rscript differential_atac.R <XML file>", call.=FALSE)
}
suppressPackageStartupMessages({
library(rtracklayer)
library(csaw)
library(edgeR)
library(XML)
})

#parse from XML
data = xmlToList(xmlParse(args[1]))
label1 = data[["File1"]][["Name"]]
bam1.R1 = data [["File1"]][["Bam.R1"]]
bam1.R2 = data [["File1"]][["Bam.R2"]]
label2 = data [["File2"]][["Name"]]
bam2.R1 = data [["File2"]][["Bam.R1"]]
bam2.R2 = data [["File2"]][["Bam.R2"]]
blacklistfile = data[["Blacklistfile"]]

#check whether blacklist file exists
if(!file.exists(blacklistfile)){
stop("Blacklistfile not found!!")
}

#create vector
bam.files = c(bam1.R1, bam1.R2, bam2.R1, bam2.R2)
celltype=c(label1,label1,label2,label2)
#dataframe
dd=data.frame(BAM=bam.files,CellType=celltype)

#outputheading
heading=paste(c(label1,'-vs-',label2),collapse='')

sink(paste(heading,"_output",sep=""))
sink(paste(heading,"_output",sep=""))
getwd()
list.files(path='.')
print("Dataframe")
print(dd)

#### blacklist region
data <- read.table(blacklistfile,header=F) 
colnames(data) <- c('chr','start','end','id','score','strand')
blacklist_granges <- with(data, GRanges(chr, IRanges(start+1, end),strand, score, id=id))
head(blacklist_granges)  
#saveRDS(file = "hg19-blacklist.rds", blacklist_granges)  

#Correlate reads
param <- readParam(discard=blacklist_granges) 
x <- correlateReads(bam.files,param=param)
frag.len <- which.max(x) - 1
print("Frag. Length")
print(frag.len)
#plot number of reads
plotno=1
png(paste(c(heading,'-',plotno,'.distance_between_cuts.png'),collapse=''))
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")
dev.off()
plotno=plotno+1

#Computing insert size histogram.
png(paste(c(heading,'-',plotno,'.fragmentsize.png'),collapse=''))
out <- getPESizes(bam.files[1])
frag.sizes <- out$sizes[out$sizes <= 1000]
hist(frag.sizes, breaks = 50, xlab = "Fragment sizes (bp)", ylab = "Frequency", main = "")
dev.off()
plotno=plotno+1

win.data <- windowCounts(bam.files, width=150L, ext=frag.len,param=param)
print("win.data ")
print(win.data)
print("wind.data$totals")
print(win.data$totals)
bins <- windowCounts(bam.files, bin=TRUE, width=2000L,param=param)
filter.stat <- filterWindows(win.data, background=bins, type="global")
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
print("summary")
summary(keep)


#pdf(file = "adjusted-bin-log-cpm.pdf")
png(paste(c(heading,'-',plotno,'.Adjusted bin log-CPM.png'),collapse=''))
hist(filter.stat$back.abundances, xlab="Adjusted bin log-CPM", breaks=100, main="", 
		   xlim=c(min(filter.stat$back.abundances), 0))
global.bg <- filter.stat$abundances - filter.stat$filter
abline(v=global.bg[1], col="red", lwd=2)
abline(v=global.bg[1]+log2(min.fc), col="blue", lwd=2)
legend("topright", lwd=2, col=c('red', 'blue'), legend=c("Background", "Threshold"))
dev.off()
plotno=plotno+1

##filtering using CPM cutofff
abundances <- aveLogCPM(asDGEList(win.data))
summary(abundances)
keep1 <- abundances > aveLogCPM(5, lib.size=mean(win.data$totals))
summary(keep)
#filtered.data <- win.data[keep,]

filtered.data <- win.data[keep & keep1,]
gc()

row.sums <- rowSums(assay(filtered.data))
print("min row sums")
print(min(row.sums))
png(paste(c(heading,'-',plotno,'.histogram of row sum.png'),collapse=''))
hist(row.sums[row.sums < 100], breaks = 30)
dev.off()
plotno=plotno+1

#Normalization
win.ab <- filter.stat$abundances[keep]
adjc <- log2(assay(filtered.data)+0.5)
logfc <- adjc[,1] - adjc[,4]
png(paste(c(heading,'-',plotno,'.log-fold-changes-before-norm.png'),collapse=''))
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),xlab="Average abundance", ylab="Log-fold change")
plotno=plotno+1
dev.off()

offsets <- normOffsets(filtered.data, type="loess")
print("offsets")
print(head(offsets))

#pdf(file = "log-fold-changes-after-normalization.pdf")
png(paste(c(heading,'-',plotno,'.log-fold-changes-after-norm.png'),collapse=''))
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
			 xlab="Average abundance", ylab="Log-fold change")
dev.off()
plotno=plotno+1

celltype <- factor(celltype)
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
print("design")
print(design)


y <- asDGEList(filtered.data)
y$offset <- offsets
y <- estimateDisp(y, design)
print("summary")
summary(y$trended.dispersion)

png(paste(c(heading,'-',plotno,'.BCV.png'),collapse=''))
plotBCV(y)
dev.off()
plotno=plotno+1

#fitting the model
fit <- glmQLFit(y, design, robust=TRUE)
#pdf(file = "qldisp.pdf")
png(paste(c(heading,'-',plotno,'.QLDispersion.png'),collapse=''))
plotQLDisp(fit)
dev.off()
plotno=plotno+1

print("summary(fit$dif.prior)")
print(summary(fit$df.prior))

#MDS Plot
#pdf(file = "mds.pdf")
png(paste(c(heading,'-',plotno,'.MDS.png'),collapse=''))
plotMDS(norm.adjc, labels=bam.files,
		  col=c("red", "blue")[as.integer(celltype)])
dev.off()
plotno=plotno+1

print("design")
print(design)

contrast <- makeContrasts(paste(c(label1,'-',label2),collapse=''), levels=design)
res <- glmQLFTest(fit, contrast=contrast)
print("head F and P value table")
print(head(res$table))

merged <- mergeWindows(rowRanges(filtered.data), tol=100, max.width=5000)
tabcom <- combineTests(merged$id, res$table)
print("head P value and FDR")
print(head(tabcom))

#Multiple test correction
is.sig <- tabcom$FDR <= 0.05
print("Summary")
print(summary(is.sig))

print("Total-Up-Down-both Table")
has.up <- tabcom$logFC.up > 0
has.down <- tabcom$logFC.down > 0
dafr=data.frame(Total.DB=sum(is.sig), DB.up=sum(is.sig & has.up & !has.down),
	 DB.down=sum(is.sig & !has.up & has.down), DB.both=sum(is.sig & has.up & has.down))
print(dafr)

out.ranges <- merged$region

simplified_up <- out.ranges[is.sig & has.up & !has.down]
simplified_up$score <- -10*log10(tabcom$PValue[is.sig & has.up & !has.down])

simplified_down <- out.ranges[is.sig & !has.up & has.down]
simplified_down$score <- -10*log10(tabcom$PValue[is.sig & !has.up & has.down])

simplified_both <- out.ranges[is.sig & has.up & has.down]
simplified_both$score <- -10*log10(tabcom$PValue[is.sig & has.up & has.down])

export(con=paste(c(heading,"-up-fdr.bed"),collapse=''),object=simplified_up)
export(con=paste(c(heading,"-down-fdr.bed"),collapse=''), object=simplified_down)
export(con=paste(c(heading,"-both-fdr.bed"),collapse=''), object=simplified_both)

sink()
#save.session("RSession.Rda")
