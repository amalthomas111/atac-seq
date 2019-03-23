#!/usr/bin/env Rscript
#Scirpt to plot insert size
#A.T
args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)==0)
 {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}
#print(args[1])
mydata = scan(args[1])
#summary(mydata)
outputname = gsub(".txt","",basename(args[1]))
mainDir=getwd()
subDir="quality"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

#histogram
png(filename=paste0("quality/hist_",outputname,".png"))
hist(mydata, main="Insert Size Distribution", xlab="Insert size (bp)",
     ylab="No. of reads", breaks=2000,
     xlim=c(0,600))
abline(v=38)
dev.off()

#density plot
png(filename=paste0("quality/hist_",outputname,"_density.png"))
plot(density(mydata,from=0,to=quantile(x=mydata,probs=0.99)),
     xlab="Insert Size (bp)",main="Insert Size Distribution")
dev.off()
