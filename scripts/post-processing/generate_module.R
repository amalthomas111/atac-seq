#!/usr/bin/env Rscript
#Rscript for performing module analysis
#input file countdata(no of reads falling into region of interest)
args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)==0)
{
  stop("At least one argument must be supplied (input file)", call.=FALSE)
} 
if (!file.exists(args[1]))
{
  stop("Input file not found!\
Usage: Rscript find_module.R <count-data-file>", call.=FALSE)
}
#print(args[1])
suppressPackageStartupMessages({
library(WGCNA)
library(DESeq2)
library(ggplot2)
})
allowWGCNAThreads(nThreads = 4)
#countdata = scan(args[1])
atacData <- read.table(args[1])#,header=TRUE)
head(atacData)
#create designmatrix
celltype <- colnames(atacData)
design.matrix <- data.frame( celltype = gsub('R','',celltype)) 
rownames(design.matrix) <- colnames(atacData)
design.matrix

dds <- DESeqDataSetFromMatrix(countData = atacData,
                              colData = design.matrix,
                              design = ~ celltype)
#variable stablization
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
#plot pca
pdf(file="pca.pdf")
plotPCA(vsd, intgroup="celltype") + theme_bw() 
datExpr = as.data.frame(t(assay(vsd)))
dev.off()

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK  

pdf(file="hclust-plot.pdf") 
sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 
dev.off()
#wgcna power
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "signed", minModuleSize = 5,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       networkType = "signed",
                       verbose = 3) 
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]] 

samples <- rownames(datExpr)

module.plot <- function(which.module) {
  sizeGrWindow(8,7)
  ME <- MEs[, paste("ME",which.module, sep="")]
  
  pdf(paste("Module-", which.module, ".pdf", sep=""))
  par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[,moduleLabels==which.module ]) ),
          nrgcols=30,rlabels=F, clabels=samples,
          main="", cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab=paste("Module-", which.module, sep=""))
  dev.off()
}

module.ids <- as.character(0:(ncol(MEs) - 1))
lapply(module.ids, module.plot) 

moduleLabels <- ordered(moduleLabels, levels = c(0:35))
genes.and.colors <- data.frame(gene=colnames(datExpr), color=moduleColors, label=moduleLabels)
genes.and.colors <- genes.and.colors[with(genes.and.colors, order(label)), ]
write.table(genes.and.colors, file = "gene-color-label.txt", quote = FALSE, sep = '\t', row.names = FALSE) 

datExpr.reord <- datExpr[, match(genes.and.colors$gene, colnames(datExpr))]
#datExpr <- datExpr[, which(colnames(datExpr) %in% diff.expr.genes)]
write.table(t(datExpr.reord), file="diff-genes-expression.txt", quote=FALSE, sep = '\t')  
pdf(file = "diff-expr-heatmap.pdf")
plotMat(t(scale(datExpr.reord)),
        nrgcols=30,rlabels=rep("*", 1, length(genes.and.colors$color)),
        clabels=samples,
        rcols=genes.and.colors$color,
        main="", cex.main=2)
dev.off() 