rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
load("2016-01-19/FSTR02.SIG.RData")
library(gage)
kegg.gs <- kegg.gsets(species = "ko", id.type = "kegg") # "ko":reference dataset
kegg.met <- kegg.gs$kg.sets[kegg.gs$met.idx] # metabolic pathway set
kegg.gs <- kegg.gs[[1]] # all pathways, without subset indice

# Metabolic GSEA for geneset with FDR <0.01
# HYPERGEOMETRIC TEST
# pr = 1 - phyper(q,m,n,k)
# N: total no. of  genes = m+n
# q: no. of genes "hit" within a certain gene set 
# m: no. of DE genes or comply with a certain standard genes
# n: N-m
# k: no. of genes within a certain gene set
FS.sig.1[,1] <- as.numeric(as.character(FS.sig.1[,1]))
FS.sig.1[,2] <- as.numeric(as.character(FS.sig.1[,2]))
FS.sig.1[,3] <- as.numeric(as.character(FS.sig.1[,3]))
FS.CC.N.1 <- rownames(FS.sig.1)[which(FS.sig.1[,1]>FS.sig.1[,2] & FS.sig.1[,1]>FS.sig.1[,3] & FS.sig.1[,1]>=1)]
FS.LC.N.1 <- rownames(FS.sig.1)[which(FS.sig.1[,2]>FS.sig.1[,1] & FS.sig.1[,2]>FS.sig.1[,3] & FS.sig.1[,2]>=1)]
FS.SC.N.1 <- rownames(FS.sig.1)[which(FS.sig.1[,3]>FS.sig.1[,2] & FS.sig.1[,3]>FS.sig.1[,1] & FS.sig.1[,3]>=1)]
FS.GSEA.1 <- data.frame("Metabolic.GS"=names(kegg.met),
                        "CC.q"=0,"CC.m"=length(FS.CC.N.1),"CC.n"=nrow(FS.sig.1)-length(FS.CC.N.1),"CC.k"=0,"CC.pr"=NA,
                        "LC.q"=0,"LC.m"=length(FS.LC.N.1),"LC.n"=nrow(FS.sig.1)-length(FS.LC.N.1),"LC.k"=0,"LC.pr"=NA,
                        "SC.q"=0,"SC.m"=length(FS.SC.N.1),"SC.n"=nrow(FS.sig.1)-length(FS.SC.N.1),"SC.k"=0,"SC.pr"=NA)
for(i in 1:nrow(FS.GSEA.1)){
  FS.GSEA.1$CC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.CC.N.1)))
  FS.GSEA.1$CC.k[i] <- length(kegg.met[[i]])
  FS.GSEA.1$LC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.LC.N.1)))
  FS.GSEA.1$LC.k[i] <- length(kegg.met[[i]])
  FS.GSEA.1$SC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.SC.N.1)))
  FS.GSEA.1$SC.k[i] <- length(kegg.met[[i]])
}
FS.GSEA.1$CC.pr <- 1 - phyper(FS.GSEA.1$CC.q,FS.GSEA.1$CC.m,FS.GSEA.1$CC.n,FS.GSEA.1$CC.k)
FS.GSEA.1$LC.pr <- 1 - phyper(FS.GSEA.1$LC.q,FS.GSEA.1$LC.m,FS.GSEA.1$LC.n,FS.GSEA.1$LC.k)
FS.GSEA.1$SC.pr <- 1 - phyper(FS.GSEA.1$SC.q,FS.GSEA.1$SC.m,FS.GSEA.1$SC.n,FS.GSEA.1$SC.k)


# Metabolic GSEA for geneset with FDR <0.05
# HYPERGEOMETRIC TEST
# pr = 1 - phyper(q,m,n,k)
# N: total no. of  genes = m+n
# q: no. of genes "hit" within a certain gene set 
# m: no. of DE genes or comply with a certain standard genes
# n: N-m
# k: no. of genes within a certain gene set
FS.sig.5[,1] <- as.numeric(as.character(FS.sig.5[,1]))
FS.sig.5[,2] <- as.numeric(as.character(FS.sig.5[,2]))
FS.sig.5[,3] <- as.numeric(as.character(FS.sig.5[,3]))
FS.CC.N.5 <- rownames(FS.sig.5)[which(FS.sig.5[,1]>FS.sig.5[,2] & FS.sig.5[,1]>FS.sig.5[,3] & FS.sig.5[,1]>=1)]
FS.LC.N.5 <- rownames(FS.sig.5)[which(FS.sig.5[,2]>FS.sig.5[,1] & FS.sig.5[,2]>FS.sig.5[,3] & FS.sig.5[,2]>=1)]
FS.SC.N.5 <- rownames(FS.sig.5)[which(FS.sig.5[,3]>FS.sig.5[,2] & FS.sig.5[,3]>FS.sig.5[,1] & FS.sig.5[,3]>=1)]
FS.GSEA.5 <- data.frame("Metabolic.GS"=names(kegg.met),
                        "CC.q"=0,"CC.m"=length(FS.CC.N.5),"CC.n"=nrow(FS.sig.5)-length(FS.CC.N.5),"CC.k"=0,"CC.pr"=NA,
                        "LC.q"=0,"LC.m"=length(FS.LC.N.5),"LC.n"=nrow(FS.sig.5)-length(FS.LC.N.5),"LC.k"=0,"LC.pr"=NA,
                        "SC.q"=0,"SC.m"=length(FS.SC.N.5),"SC.n"=nrow(FS.sig.5)-length(FS.SC.N.5),"SC.k"=0,"SC.pr"=NA)
for(i in 1:nrow(FS.GSEA.5)){
  FS.GSEA.5$CC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.CC.N.5)))
  FS.GSEA.5$CC.k[i] <- length(kegg.met[[i]])
  FS.GSEA.5$LC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.LC.N.5)))
  FS.GSEA.5$LC.k[i] <- length(kegg.met[[i]])
  FS.GSEA.5$SC.q[i] <- sum(!is.na(match(kegg.met[[i]],FS.SC.N.5)))
  FS.GSEA.5$SC.k[i] <- length(kegg.met[[i]])
}
FS.GSEA.5$CC.pr <- 1 - phyper(FS.GSEA.5$CC.q,FS.GSEA.5$CC.m,FS.GSEA.5$CC.n,FS.GSEA.5$CC.k)
FS.GSEA.5$LC.pr <- 1 - phyper(FS.GSEA.5$LC.q,FS.GSEA.5$LC.m,FS.GSEA.5$LC.n,FS.GSEA.5$LC.k)
FS.GSEA.5$SC.pr <- 1 - phyper(FS.GSEA.5$SC.q,FS.GSEA.5$SC.m,FS.GSEA.5$SC.n,FS.GSEA.5$SC.k)
