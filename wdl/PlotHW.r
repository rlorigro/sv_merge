#!/usr/bin/env Rscript
#
# Just a copy of the following script, with minor fixes:
# <https://github.com/broadinstitute/gatk-sv/blob/bd85b0a87b8c751759a92aa97ef7ba375cbbf8db/src/sv-pipeline/scripts/vcf_qc/plot_sv_vcf_distribs.R#L1018>.
#
# It runs without problems in the following cloud environment: 
# Legacy R / Bioconductor (R 4.2.2, Bioconductor 3.16, Python 3.7.12) 
#
# Rscript PlotHW.r input output
# @param args 0=input GT count matrix; must have header $AA,AB,BB$, where
# AA=homref, AB=het, BB=homvar; 1=output PNG.
#
args=commandArgs(trailingOnly=TRUE)
install.packages(pkgs='HardyWeinberg', repos="http://cran.us.r-project.org", verbose=TRUE, Ncpus=4)
library(HardyWeinberg)

png(args[2],res=300,height=1800,width=1800)

lab.cex=1

# Prep HW matrix
HW.mat <- as.matrix(read.table(args[1], sep=",", header=TRUE))

# Gather HW p-values & colors
HW.p <- HWChisqStats(X=HW.mat, x.linked=F, pvalues=T)
HW.cols <- rep("#4DAC26", times=length(HW.p))
HW.cols[which(HW.p<0.05)] <- "#81F850"
HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"

# Generate HW plot frame
par(mar=c(1,3.5,3,0.5),bty="n")
plot(x=1.15*c(-1/sqrt(3),1/sqrt(3)),y=c(-0.15,1.15),type="n",xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
segments(x0=c(-1/sqrt(3),0,1/sqrt(3)), x1=c(0,1/sqrt(3),-1/sqrt(3)), y0=c(0,1,0),y1=c(1,0,0))
HWTernaryPlot(X=HW.mat,newframe=F, vbounds=F,mafbounds=F, region=1,vertexlab=NA, alpha=0.05, curvecols=c("#4DAC26","#81F850",NA,NA),pch=NA)

# Add axes
text(x=c(-1/sqrt(3),1/sqrt(3)),y=0,labels=c("Ref.","Hom."),pos=1,font=4)
text(x=0,y=1,labels="Het.",pos=3,font=4)

# Finish HW plot
HWTernaryPlot(X=HW.mat,newframe=F, vbounds=F,mafbounds=F, region=1,vertexlab=NA, alpha=0.05/nrow(HW.mat), curvecols=c("#4DAC26","#AC26A1",NA,NA), pch=21,cex=0.5,signifcolour=F,markercol=HW.cols, markerbgcol=adjustcolor(HW.cols,alpha=0.25))

# Add legend. Treating NaN p-values as in HW. Calls with NaN have GT counts
# $1027,0,0$, so they might have been filtered out before this script in
# gatk-sv.
nan_flags = is.nan(HW.p)
n.pass <- length(which(HW.p>=0.05)) + length(which(nan_flags))
n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HW.mat)))
n.bonf <- length(which(HW.p<0.05/nrow(HW.mat)))
legend("topright",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=2, legend=c(paste("SV in H-W equilibrium\n(n=", prettyNum(n.pass,big.mark=","),"; ", round(100*(n.pass/nrow(HW.mat)),2),"%)\n",sep=""), paste("SV not in H-W equilibrium\n(Nominal; n=", prettyNum(n.nom,big.mark=","),"; ", round(100*(n.nom/nrow(HW.mat)),2),"%)\n",sep=""), paste("SV not in H-W equilibrium\n(Bonferroni; n=", prettyNum(n.bonf,big.mark=","),"; ", round(100*(n.bonf/nrow(HW.mat)),2),"%)\n",sep="")), bty="n",bg=NA,cex=0.7)

# Add number of SV on plot
axis(3,at=mean(par("usr")[1:2]),line=-0.9,tick=F,cex.axis=0.8, labels=paste("n=",prettyNum(nrow(HW.mat),big.mark=","),sep=""))

dev.off()