rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
load("2016-01-19/FSTR01.import.data.RData")
#plot(density(as.numeric(GC.mean$CC)),ylim=c(0,0.035))
#lines(density(as.numeric(GC.mean$SC)),col=2)
#lines(density(as.numeric(GC.mean$LC)),col=4)

# ANOVA test
for(i in 1:nrow(GC.mean)){
  GC.mean$pval[i] <- summary(aov(FS.mean[,i]~Gut.Comp))[[1]][1,5]
}

# FDR adj
# p.adjust(p, method = p.adjust.methods, n = length(p))
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
GC.mean$FDR <- p.adjust(GC.mean$pval,method="fdr")

# Significant threshold: FDR < 0.01
FS.sig.1 <- GC.mean[which(GC.mean$FDR<0.01),]
# Significant threshold: FDR < 0.05
FS.sig.5 <- GC.mean[which(GC.mean$FDR<0.05),]

# Save RData
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save(GC.mean,FS.sig.1,FS.sig.5,file=paste0(Sys.Date(),"/FSTR02.SIG.RData"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save(GC.mean,FS.sig.1,FS.sig.5,file=paste0(Sys.Date(),"/FSTR02.SIG.RData"))
}
  

# venn diagram?
# load("2016-01-19/FSTR02.SIG.RData")
library(VennDiagram)
venn.diagram(x=list(Cecum=lib.cc,LargeIntestine=lib.lc,SmallIntestine=lib.sc),filename="ts_venn.tiff",
             fill=c("cornflowerblue","green","yellow"),alpha = 0.50,
             fontfamily = "serif", fontface = "bold",cat.cex = 1.5,cat.fontfamily = "serif")