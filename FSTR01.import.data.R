rm(list=ls())
setwd("~/FS_transcriptome/kaas/KO")
FSfile <- list.files()
file.names <- gsub(".txt","",FSfile)
FS.gut <- lapply(FSfile,read.table)


for(i in 1:length(FS.gut)){
  FS.gut[[i]] <- as.data.frame(table(FS.gut[[i]][,2]))
  colnames(FS.gut[[i]])[2] <- file.names[i]
}

FS.table <- merge.data.frame(FS.gut[[1]],FS.gut[[2]],all=TRUE)
for(j in 3:length(FS.gut)){
  FS.table <- merge.data.frame(FS.table,FS.gut[[j]],all=TRUE)
}
rownames(FS.table) <- FS.table[,1]
#colnames(FS.table)[1] <- "Kegg.ID"

# mean by gut compartments
FS.mean <- t(FS.table[,-1])
FS.mean[which(is.na(FS.mean))]=0
Gut.Comp <- substr(file.names,5,6)
GC.mean <- aggregate(FS.mean,list(Gut.Comp),mean)
GC.sd <- aggregate(FS.mean,list(Gut.Comp),sd)


GC.mean <- t(GC.mean)
colnames(GC.mean) <- GC.mean[1,]
GC.mean <- as.data.frame(GC.mean[-1,])

# Save RData
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save.image(paste0("~/FS_transcriptome/GSEA/",Sys.Date(),"/FSTR01.import.data.RData"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save.image(paste0("~/FS_transcriptome/GSEA/",Sys.Date(),"/FSTR01.import.data.RData"))
}

save(FS.mean,Gut.Comp,file="2016-01-19/FS.mean.RData")

