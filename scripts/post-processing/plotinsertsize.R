#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
if (length(args)==0)
 {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 
#print(args[1])
mydata = scan(args[1])
#summary(mydata)

#histogram
png(filename=paste0("hist_",args[1],".png"))
hist(mydata, main="Insert Size Distribution", xlab="Insert size (bp)", 
     ylab="No. of reads", breaks=2000,
     xlim=c(0,600))
abline(v=38)
dev.off()

#density plot
png(filename=paste0("hist_",args[1],"_density.png"))
plot(density(mydata,from=0,to=quantile(x=mydata,probs=0.99)),xlab="Insert Size (bp)",main="Insert Size Distribution")
dev.off()
