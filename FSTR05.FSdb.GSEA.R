rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
load("2016-01-22/FSdb.GS.RData") # load FS GeneSet Database
load("2016-01-22/FSdb.MET.RData") # load FS MetabolicSet Database
load("2016-01-19/FSTR02.SIG.RData") # load differential expression lists w/o TukeyHSD test

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
FS.GSEA.1 <- data.frame("Metabolic.GS"=names(met.db.cc),
                        "CC.q"=0,"CC.m"=length(FS.CC.N.1),"CC.n"=nrow(FS.sig.1)-length(FS.CC.N.1),"CC.k"=0,"CC.pr"=NA,
                        "LC.q"=0,"LC.m"=length(FS.LC.N.1),"LC.n"=nrow(FS.sig.1)-length(FS.LC.N.1),"LC.k"=0,"LC.pr"=NA,
                        "SC.q"=0,"SC.m"=length(FS.SC.N.1),"SC.n"=nrow(FS.sig.1)-length(FS.SC.N.1),"SC.k"=0,"SC.pr"=NA)
for(i in 1:nrow(FS.GSEA.1)){
  FS.GSEA.1$CC.q[i] <- sum(!is.na(match(met.db.cc[[i]],FS.CC.N.1)))
  FS.GSEA.1$CC.k[i] <- length(met.db.cc[[i]])
  FS.GSEA.1$LC.q[i] <- sum(!is.na(match(met.db.lc[[i]],FS.LC.N.1)))
  FS.GSEA.1$LC.k[i] <- length(met.db.lc[[i]])
  FS.GSEA.1$SC.q[i] <- sum(!is.na(match(met.db.sc[[i]],FS.SC.N.1)))
  FS.GSEA.1$SC.k[i] <- length(met.db.sc[[i]])
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
FS.GSEA.5 <- data.frame("Metabolic.GS"=names(met.db.cc),
                        "CC.q"=0,"CC.m"=length(FS.CC.N.5),"CC.n"=nrow(FS.sig.5)-length(FS.CC.N.5),"CC.k"=0,"CC.pr"=NA,
                        "LC.q"=0,"LC.m"=length(FS.LC.N.5),"LC.n"=nrow(FS.sig.5)-length(FS.LC.N.5),"LC.k"=0,"LC.pr"=NA,
                        "SC.q"=0,"SC.m"=length(FS.SC.N.5),"SC.n"=nrow(FS.sig.5)-length(FS.SC.N.5),"SC.k"=0,"SC.pr"=NA)
for(i in 1:nrow(FS.GSEA.5)){
  FS.GSEA.5$CC.q[i] <- sum(!is.na(match(met.db.cc[[i]],FS.CC.N.5)))
  FS.GSEA.5$CC.k[i] <- length(met.db.cc[[i]])
  FS.GSEA.5$LC.q[i] <- sum(!is.na(match(met.db.lc[[i]],FS.LC.N.5)))
  FS.GSEA.5$LC.k[i] <- length(met.db.lc[[i]])
  FS.GSEA.5$SC.q[i] <- sum(!is.na(match(met.db.sc[[i]],FS.SC.N.5)))
  FS.GSEA.5$SC.k[i] <- length(met.db.sc[[i]])
}
FS.GSEA.5$CC.pr <- 1 - phyper(FS.GSEA.5$CC.q,FS.GSEA.5$CC.m,FS.GSEA.5$CC.n,FS.GSEA.5$CC.k)
FS.GSEA.5$LC.pr <- 1 - phyper(FS.GSEA.5$LC.q,FS.GSEA.5$LC.m,FS.GSEA.5$LC.n,FS.GSEA.5$LC.k)
FS.GSEA.5$SC.pr <- 1 - phyper(FS.GSEA.5$SC.q,FS.GSEA.5$SC.m,FS.GSEA.5$SC.n,FS.GSEA.5$SC.k)

##### FISHER'S EXACT TEST #####
# Metabolic GSEA for geneset with FDR <0.01
# FISHER'S EXACT TEST
# fisher.test(matrix(c(q,M,K,N),2,2),alternative='greater')$p.val
# N: total no. of  genes = m+n
# q: no. of genes "hit" within a certain gene set 
# m: no. of DE genes or comply with a certain standard genes
# n: N-m
# k: no. of genes within a certain gene set
# parameters for fisher'test: M = m-q, K = k-q
FS.FISHER.1 <- data.frame("Metabolic.GS"=names(met.db.cc),
                          "CC.q"=0,"CC.M"=0,"CC.K"=0,"CC.N"=nrow(FS.sig.1),"CC.pr"=NA,
                          "LC.q"=0,"LC.M"=0,"LC.K"=0,"LC.N"=nrow(FS.sig.1),"LC.pr"=NA,
                          "SC.q"=0,"SC.M"=0,"SC.K"=0,"SC.N"=nrow(FS.sig.1),"SC.pr"=NA)
for(i in 1:nrow(FS.FISHER.1)){
  FS.FISHER.1$CC.q[i] <- sum(!is.na(match(met.db.cc[[i]],FS.CC.N.1)))
  FS.FISHER.1$CC.M[i] <- length(FS.CC.N.1)-FS.FISHER.1$CC.q[i]
  FS.FISHER.1$CC.K[i] <- length(met.db.cc[[i]])-FS.FISHER.1$CC.q[i]
  FS.FISHER.1$LC.q[i] <- sum(!is.na(match(met.db.lc[[i]],FS.LC.N.1)))
  FS.FISHER.1$LC.M[i] <- length(FS.LC.N.1)-FS.FISHER.1$LC.q[i]
  FS.FISHER.1$LC.K[i] <- length(met.db.lc[[i]])-FS.FISHER.1$LC.q[i]
  FS.FISHER.1$SC.q[i] <- sum(!is.na(match(met.db.sc[[i]],FS.SC.N.1)))
  FS.FISHER.1$SC.M[i] <- length(FS.SC.N.1)-FS.FISHER.1$SC.q[i]
  FS.FISHER.1$SC.K[i] <- length(met.db.sc[[i]])-FS.FISHER.1$SC.q[i]
  # FISHER'S EXACT TEST
  FS.FISHER.1$CC.pr[i] <- fisher.test(matrix(c(FS.FISHER.1$CC.q[i],FS.FISHER.1$CC.M[i],FS.FISHER.1$CC.K[i],FS.FISHER.1$CC.N[i]),2,2),
                                      alternative='greater')$p.val
  FS.FISHER.1$LC.pr[i] <- fisher.test(matrix(c(FS.FISHER.1$LC.q[i],FS.FISHER.1$LC.M[i],FS.FISHER.1$LC.K[i],FS.FISHER.1$LC.N[i]),2,2),
                                      alternative='greater')$p.val
  FS.FISHER.1$SC.pr[i] <- fisher.test(matrix(c(FS.FISHER.1$SC.q[i],FS.FISHER.1$SC.M[i],FS.FISHER.1$SC.K[i],FS.FISHER.1$SC.N[i]),2,2),
                                      alternative='greater')$p.val
}

# Metabolic GSEA for geneset with FDR <0.05
# FISHER'S EXACT TEST
# fisher.test(matrix(c(q,M,K,N),2,2),alternative='greater')$p.val
# N: total no. of  genes = m+n
# q: no. of genes "hit" within a certain gene set 
# m: no. of DE genes or comply with a certain standard genes
# n: N-m
# k: no. of genes within a certain gene set
# parameters for fisher'test: M = m-q, K = k-q
FS.FISHER.5 <- data.frame("Metabolic.GS"=names(met.db.cc),
                          "CC.q"=0,"CC.M"=0,"CC.K"=0,"CC.N"=nrow(FS.sig.5),"CC.pr"=NA,
                          "LC.q"=0,"LC.M"=0,"LC.K"=0,"LC.N"=nrow(FS.sig.5),"LC.pr"=NA,
                          "SC.q"=0,"SC.M"=0,"SC.K"=0,"SC.N"=nrow(FS.sig.5),"SC.pr"=NA)
for(i in 1:nrow(FS.FISHER.5)){
  FS.FISHER.5$CC.q[i] <- sum(!is.na(match(met.db.cc[[i]],FS.CC.N.5)))
  FS.FISHER.5$CC.M[i] <- length(FS.CC.N.5)-FS.FISHER.5$CC.q[i]
  FS.FISHER.5$CC.K[i] <- length(met.db.cc[[i]])-FS.FISHER.5$CC.q[i]
  FS.FISHER.5$LC.q[i] <- sum(!is.na(match(met.db.lc[[i]],FS.LC.N.5)))
  FS.FISHER.5$LC.M[i] <- length(FS.LC.N.5)-FS.FISHER.5$LC.q[i]
  FS.FISHER.5$LC.K[i] <- length(met.db.lc[[i]])-FS.FISHER.5$LC.q[i]
  FS.FISHER.5$SC.q[i] <- sum(!is.na(match(met.db.sc[[i]],FS.SC.N.5)))
  FS.FISHER.5$SC.M[i] <- length(FS.SC.N.5)-FS.FISHER.5$SC.q[i]
  FS.FISHER.5$SC.K[i] <- length(met.db.sc[[i]])-FS.FISHER.5$SC.q[i]
  # FISHER'S EXACT TEST
  FS.FISHER.5$CC.pr[i] <- fisher.test(matrix(c(FS.FISHER.5$CC.q[i],FS.FISHER.5$CC.M[i],FS.FISHER.5$CC.K[i],FS.FISHER.5$CC.N[i]),2,2),
                                      alternative='greater')$p.val
  FS.FISHER.5$LC.pr[i] <- fisher.test(matrix(c(FS.FISHER.5$LC.q[i],FS.FISHER.5$LC.M[i],FS.FISHER.5$LC.K[i],FS.FISHER.5$LC.N[i]),2,2),
                                      alternative='greater')$p.val
  FS.FISHER.5$SC.pr[i] <- fisher.test(matrix(c(FS.FISHER.5$SC.q[i],FS.FISHER.5$SC.M[i],FS.FISHER.5$SC.K[i],FS.FISHER.5$SC.N[i]),2,2),
                                      alternative='greater')$p.val
}

##### DISCUSSION #####
# PROBLEM:: (1-1) too large DE gene datasets [To Solve!] (1-2) post-hoc multiple test
# PROBLEM:: (2) gene lists in metabolic sets are small [Needn't solve!]


# Save Metabolic GSEA tables
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  write.csv(FS.GSEA.1,file=paste0(Sys.Date(),"/FS.GSEA.1.csv"))
  write.csv(FS.GSEA.5,file=paste0(Sys.Date(),"/FS.GSEA.5.csv"))
  write.csv(FS.FISHER.1,file=paste0(Sys.Date(),"/FS.FISHER.1.csv"))
  write.csv(FS.FISHER.5,file=paste0(Sys.Date(),"/FS.FISHER.5.csv"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  write.csv(FS.GSEA.1,file=paste0(Sys.Date(),"/FS.GSEA.1.csv"))
  write.csv(FS.GSEA.5,file=paste0(Sys.Date(),"/FS.GSEA.5.csv"))
  write.csv(FS.FISHER.1,file=paste0(Sys.Date(),"/FS.FISHER.1.csv"))
  write.csv(FS.FISHER.5,file=paste0(Sys.Date(),"/FS.FISHER.5.csv"))
}