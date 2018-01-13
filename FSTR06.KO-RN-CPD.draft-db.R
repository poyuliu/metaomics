rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
library(KEGGgraph)
library(gage)
kegg.gs <- kegg.gsets(species = "ko", id.type = "kegg") # "ko":reference dataset
kegg.met <- kegg.gs$kg.sets[kegg.gs$met.idx] # metabolic pathway set
kegg.cpd <- kegg.met # empty metabolite/compound libraries
ko.met <- substr(names(kegg.met),1,7)

for(x in 1:length(kegg.met)){
  kgml <- parseKGML(paste0("http://rest.kegg.jp/get/",ko.met[x],"/kgml"))
#  kgml <- parseKGML(paste0("http://www.kegg.jp/kegg-bin/download?entry=",ko.met[x],"&format=kgml"))
  
  reaction <- getReactions(kgml)
  node <-nodes(kgml)

  # REACTION LIST
  if(length(reaction)==0){
    rx <- data.frame(rn=NA,type=NA,substrate.1=NA,substrate.2=NA,substrate.3=NA,substrate.4=NA,substrate.5=NA,product.1=NA,product.2=NA,product.3=NA,product.4=NA,product.5=NA)
  } else if(length(reaction)!=0){
    rx <- data.frame(rn=NA,type=NA,substrate.1=NA,substrate.2=NA,substrate.3=NA,substrate.4=NA,substrate.5=NA,product.1=NA,product.2=NA,product.3=NA,product.4=NA,product.5=NA)
    for(i in 1:length(reaction)){
      if(length(reaction[[i]]@substrateName)==5 && length(reaction[[i]]@productName)==5){
        rx[i,] <- c(reaction[[i]]@name,reaction[[i]]@type,reaction[[i]]@substrateName,reaction[[i]]@productName)
      } else if(length(reaction[[i]]@substrateName)==5 && length(reaction[[i]]@productName)!=5){
        rx[i,] <- c(reaction[[i]]@name,reaction[[i]]@type,reaction[[i]]@substrateName,reaction[[i]]@productName,rep(NA,5-length(reaction[[i]]@productName)))
      } else if(length(reaction[[i]]@substrateName)!=5 && length(reaction[[i]]@productName)==5){
        rx[i,] <- c(reaction[[i]]@name,reaction[[i]]@type,reaction[[i]]@substrateName,rep(NA,5-length(reaction[[i]]@substrateName)),reaction[[i]]@productName)
      } else
        rx[i,] <- c(reaction[[i]]@name,reaction[[i]]@type,reaction[[i]]@substrateName,rep(NA,5-length(reaction[[i]]@substrateName)),reaction[[i]]@productName,rep(NA,5-length(reaction[[i]]@productName)))
    }
    rx[,1] <- gsub("rn:","",rx[,1])
    rx[,3] <- gsub("cpd:","",rx[,3]);rx[,3] <- gsub("gl:","",rx[,3])
    rx[,4] <- gsub("cpd:","",rx[,4]);rx[,4] <- gsub("gl:","",rx[,4])
    rx[,5] <- gsub("cpd:","",rx[,5]);rx[,5] <- gsub("gl:","",rx[,5])
    rx[,6] <- gsub("cpd:","",rx[,6]);rx[,6] <- gsub("gl:","",rx[,6])
    rx[,7] <- gsub("cpd:","",rx[,7]);rx[,7] <- gsub("gl:","",rx[,7])
    rx[,8] <- gsub("cpd:","",rx[,8]);rx[,8] <- gsub("gl:","",rx[,8])
    rx[,9] <- gsub("cpd:","",rx[,9]);rx[,9] <- gsub("gl:","",rx[,9])
    rx[,10] <- gsub("cpd:","",rx[,10]);rx[,10] <- gsub("gl:","",rx[,10])
    rx[,11] <- gsub("cpd:","",rx[,11]);rx[,11] <- gsub("gl:","",rx[,11])
    rx[,12] <- gsub("cpd:","",rx[,12]);rx[,12] <- gsub("gl:","",rx[,12])
  }
  
  # KO+REACTION LIST
  if(length(reaction)==0){
    ko <- data.frame(idx=NA,KO=NA,rn=NA)
    for(k in 1:length(node)){
      if(k==1){
        ko <- data.frame(idx=k,KO=node[[k]]@name,rn=node[[k]]@reaction)
      } else
        ko <- rbind(ko,data.frame(idx=k,KO=node[[k]]@name,rn=node[[k]]@reaction))
    }
    cpd <- ko[which(substr(ko[,2],1,3)=="cpd"|substr(ko[,2],1,3)=="gl:"),2]; cpd <- gsub("cpd:","",cpd); cpd <- gsub("gl:","",cpd) # get compound IDs in this kgml
    ko <- ko[which(substr(ko[,2],1,3)=="ko:"),]
    ko[,2] <- gsub("ko:","",ko[,2])
    
  } else if(length(reaction)!=0){
    ko <- data.frame(idx=NA,KO=NA,rn=NA)
    for(k in 1:length(node)){
      if(k==1){
        ko <- data.frame(idx=k,KO=node[[k]]@name,rn=node[[k]]@reaction)
      } else
        ko <- rbind(ko,data.frame(idx=k,KO=node[[k]]@name,rn=node[[k]]@reaction))
    }
    cpd <- ko[which(substr(ko[,2],1,3)=="cpd"|substr(ko[,2],1,3)=="gl:"),2]; cpd <- gsub("cpd:","",cpd); cpd <- gsub("gl:","",cpd) # get compound IDs in this kgml
    ko <- ko[which(!is.na(ko[,3])),]
    ko[,2] <- gsub("ko:","",ko[,2])
    ko[,3] <- gsub("rn:","",ko[,3])
  }
  
  # merge ko+rx
  korx <- merge(ko,rx,by="rn",all=TRUE)
  korx <- korx[which(!duplicated(korx)),]
  #korx[,1] <- substr(korx[,1],1,6) # remove 2nd rx of dual reactions
  rownames(korx) <- NULL
  
  # Maunal curation IDs for reaction diretions
  p <- which(duplicated(korx[,1:3]))
  
  # renew kegg.met/kegg.cpd
  kegg.met[[x]] <- korx
  kegg.cpd[[x]] <- cpd
}

# save KO-RN-CPD.draft.db
# IS NOT "FS"-KO-RN-CPD.db
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save(kegg.met,kegg.cpd,file=paste0(Sys.Date(),"/KO-RN-CPD.draft.db.RData"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save(kegg.met,kegg.cpd,file=paste0(Sys.Date(),"/KO-RN-CPD.draft.db.RData"))
}