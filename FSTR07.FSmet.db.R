rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
load("2016-01-22/FSdb.MET.RData")
load("2016-01-25/KO-RN-CPD.draft.db.RData")

cpd.db.cc <- met.db.cc
cpd.db.lc <- met.db.lc
cpd.db.sc <- met.db.sc

# match FS transcriptome data to KO-RN-CPD.draft.db and generate FS metabolite dbs
for(i in 1:152){
  cpd.cc <- c(as.matrix(kegg.met[[i]][which(kegg.met[[i]]$KO %in% met.db.cc[[i]]),5:14]))
  cpd.db.cc[[i]] <- unique(cpd.cc[!is.na(cpd.cc)])
  cpd.lc <- c(as.matrix(kegg.met[[i]][which(kegg.met[[i]]$KO %in% met.db.lc[[i]]),5:14]))
  cpd.db.lc[[i]] <- unique(cpd.lc[!is.na(cpd.lc)])
  cpd.sc <- c(as.matrix(kegg.met[[i]][which(kegg.met[[i]]$KO %in% met.db.sc[[i]]),5:14]))
  cpd.db.sc[[i]] <- unique(cpd.sc[!is.na(cpd.sc)])
}

# output cpd list without pathway category (txt format)
cc.db <- unique(unlist(cpd.db.cc))
lc.db <- unique(unlist(cpd.db.lc))
sc.db <- unique(unlist(cpd.db.sc))

# remove "dr"
cc.db[which(substr(cc.db,1,3)=="dr:")] <- substr(cc.db[which(substr(cc.db,1,3)=="dr:")],11,16)
cc.db <- unique(cc.db)
lc.db[which(substr(lc.db,1,3)=="dr:")] <- substr(lc.db[which(substr(lc.db,1,3)=="dr:")],11,16)
lc.db <- unique(lc.db)
sc.db[which(substr(sc.db,1,3)=="dr:")] <- substr(sc.db[which(substr(sc.db,1,3)=="dr:")],11,16)
sc.db <- unique(sc.db)


# save FSmet.db.db
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save(cpd.db.cc,cpd.db.lc,cpd.db.sc,file=paste0(Sys.Date(),"/FSmet.db.RData"))
  write(cc.db,file=paste0(Sys.Date(),"/cc.db.txt"))
  write(lc.db,file=paste0(Sys.Date(),"/lc.db.txt"))
  write(sc.db,file=paste0(Sys.Date(),"/sc.db.txt"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save(cpd.db.cc,cpd.db.lc,cpd.db.sc,file=paste0(Sys.Date(),"/FSmet.db.RData"))
  write(cc.db,file=paste0(Sys.Date(),"/cc.db.txt"))
  write(lc.db,file=paste0(Sys.Date(),"/lc.db.txt"))
  write(sc.db,file=paste0(Sys.Date(),"/sc.db.txt"))
}



##################################################
# Note:
# DID NOT REMOVE COUPLED SUBSTRATE OR PRODUCT REACTIONS. ex. "ko00830 Retinol metabolism"
# DID NOT REMOVE "dr:" ids of DRUG METABOLISMs. ex. "ko00983 Drug metabolism - other enzymes"
##################################################
# check compound list with coupled IDs
# char.cc <- list()
# char.lc <- list()
# char.sc <- list()
# for(i in 1:152){
#   char.cc[[i]] <- nchar(cpd.db.cc[[i]])
#   char.lc[[i]] <- nchar(cpd.db.lc[[i]])
#   char.sc[[i]] <- nchar(cpd.db.sc[[i]])
# }
# print(char.cc)
# print(char.lc)
# print(char.sc)