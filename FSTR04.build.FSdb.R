rm(list=ls())
setwd("~/FS_transcriptome/GSEA/")
load("2016-01-19/FSTR02.SIG.RData")
library(gage)
# Load KEGG KO sets
kegg.gs <- kegg.gsets(species = "ko", id.type = "kegg") # "ko":reference dataset
kegg.met <- kegg.gs$kg.sets[kegg.gs$met.idx] # metabolic pathway set
kegg.gs <- kegg.gs[[1]] # all pathways, without subset indice

#####20160905#####
## output K IDs by pathways
c1 <- substr(khier[which(khier[,1]=="Metabolism"),3],1,5)
c2 <- substr(khier[which(khier[,1]=="Genetic Information Processing"),3],1,5)
c3 <- substr(khier[which(khier[,1]=="Environmental Information Processing"),3],1,5)
c4 <- substr(khier[which(khier[,1]=="Cellular Processes"),3],1,5)

kegg.gs <- kegg.gsets(species = "ko", id.type = "kegg") # "ko":reference dataset
kegg.met <- kegg.gs$kg.sets[substr(names(kegg.gs$kg.sets),3,7) %in% c1]# metabolic pathway set
kegg.gip <- kegg.gs$kg.sets[substr(names(kegg.gs$kg.sets),3,7) %in% c2] # Genetic Information Processing
kegg.eip <- kegg.gs$kg.sets[substr(names(kegg.gs$kg.sets),3,7) %in% c3] # Environmental Information Processing
kegg.cep <- kegg.gs$kg.sets[substr(names(kegg.gs$kg.sets),3,7) %in% c4] # Cellular Processes
kegg.gs <- kegg.gs[[1]] # all pathways, without subset indice

save(kegg.gs,kegg.met,kegg.gip,kegg.eip,kegg.cep,file="~/FS_transcriptome/GSEA/kegg_path_20160905.RData")
#####20160905#####


# Gut compartmental K ID libraries
# QUALITATIVE!!!
lib.cc <- rownames(GC.mean)[which(GC.mean$CC!=0)]
lib.lc <- rownames(GC.mean)[which(GC.mean$LC!=0)]
lib.sc <- rownames(GC.mean)[which(GC.mean$SC!=0)]
# QUANTITATIVE >>> threshold??

# venn diagram
# library(VennDiagram)
# venn.diagram(x=list(Cecum=lib.cc,LargeIntestine=lib.lc,SmallIntestine=lib.sc),filename="ts_venn.tiff",
#              fill=c("cornflowerblue","green","yellow"),alpha = 0.50,imagetype = "tiff",
#              fontfamily = "serif", fontface = "bold",cat.cex = 1.5,cat.fontfamily = "serif")



# DBs with all pathway sets (signalling, metabolic, disease pathways)
# Cecum
db.cc <- kegg.gs
for(i in 1:length(db.cc)){
  db.cc[[i]] <- db.cc[[i]][which(db.cc[[i]] %in% lib.cc)]
}
# Large intestine
db.lc <- kegg.gs
for(i in 1:length(db.lc)){
  db.lc[[i]] <- db.lc[[i]][which(db.lc[[i]] %in% lib.lc)]
}
# Small intestine
db.sc <- kegg.gs
for(i in 1:length(db.sc)){
  db.sc[[i]] <- db.sc[[i]][which(db.sc[[i]] %in% lib.sc)]
}
# Save FSdb.GS
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save(db.cc,db.lc,db.sc,file=paste0(Sys.Date(),"/FSdb.GS.RData"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save(db.cc,db.lc,db.sc,file=paste0(Sys.Date(),"/FSdb.GS.RData"))
}

# DBs with Metabolic pathway set
# Cecum
met.db.cc <- kegg.met
for(i in 1:length(met.db.cc)){
  met.db.cc[[i]] <- met.db.cc[[i]][which(met.db.cc[[i]] %in% lib.cc)]
}
# Large intestine
met.db.lc <- kegg.met
for(i in 1:length(met.db.lc)){
  met.db.lc[[i]] <- met.db.lc[[i]][which(met.db.lc[[i]] %in% lib.lc)]
}
# Small intestine
met.db.sc <- kegg.met
for(i in 1:length(met.db.sc)){
  met.db.sc[[i]] <- met.db.sc[[i]][which(met.db.sc[[i]] %in% lib.sc)]
}

# Save FSdb.MET
if(paste0("./",Sys.Date()) %in% list.dirs()){
  print("Saving folder does exist!")
  save(met.db.cc,met.db.lc,met.db.sc,file=paste0(Sys.Date(),"/FSdb.MET.RData"))
} else{
  print("Saving folder does not exist!")
  system(paste("mkdir",Sys.Date()))
  save(met.db.cc,met.db.lc,met.db.sc,file=paste0(Sys.Date(),"/FSdb.MET.RData"))
}