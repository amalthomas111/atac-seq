#!/usr/bin/env Rscript
#  Author Amal Thomas
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
#  This RScript plot the insert size distribution and density
#  Input file is insertsize file which can be created from mapped bam file
#  by the command:
#  samtools view -f66 <inputbam> |cut -f 9|sed 's/^-//' > insertsize.txt

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0)
 {
  stop("Usage: Rscript plotinsertsize.R <insertsize file>\
 To create insert file:\
 samtools view -f66 <inputbam> |cut -f 9|sed 's/^-//' > insertsize.txt", call.=FALSE)
}
if (!file.exists(args[1]))
{
  stop("Input file not found!\
Usage: Rscript plotinsertsize.R <insertsize file>", call.=FALSE)
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
