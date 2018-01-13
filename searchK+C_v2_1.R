rm(list=ls())
load("./KO-RN-CPD.draft.db.RData") # load KO+C database
# function
CKO.search <- function(KCID){
  path.idx <- grep(KCID,kegg.path)
  path.list <- kegg.path[path.idx]
  path.ko <- list()
  if(substr(KCID,1,1)=="C"){
  	for(i in 1:length(path.list)){
  	path.ko[[i]] <- path.list[[i]][which(path.list[[i]]$substrate.1==KCID|path.list[[i]]$substrate.2==KCID|path.list[[i]]$substrate.3==KCID|
  	path.list[[i]]$substrate.4==KCID|path.list[[i]]$substrate.5==KCID|path.list[[i]]$product.1==KCID|
        path.list[[i]]$product.2==KCID|path.list[[i]]$product.3==KCID|path.list[[i]]$product.4==KCID|
        path.list[[i]]$product.5==KCID),"KO"]
    path.ko[[i]] <- unique(path.ko[[i]])
    names(path.ko)[[i]] <- names(path.list)[[i]]
  	} 
  } else if(substr(KCID,1,1)=="K"){
  		for(i in 1:length(path.list)){
  			tempK <- path.list[[i]][which(path.list[[i]]$KO==KCID),5:14]
  			tempK <- tempK[!is.na(tempK)]
  			if(class(path.list[[i]])=="data.frame"){
  				path.ko[[i]] <- tempK
  			} else path.ko[[i]] <- NA
  			names(path.ko)[[i]] <- names(path.list)[[i]]
  		}
  	}
  path.ko
}

## search KO by C
#####EXAMPLE#####
# > CKO.search("C00230")
# $`ko00362 Benzoate degradation`
# [1] "K00481" "K19065" "K00448" "K00449" "K04100" "K04101"
# 
# $`ko00400 Phenylalanine, tyrosine and tryptophan biosynthesis`
# [1] "K09483" "K15652"
# 
# $`ko00623 Toluene degradation`
# [1] "K20218"
# 
# $`ko00624 Polycyclic aromatic hydrocarbon degradation`
# [1] "K19065" "K00448" "K00449" "K04100" "K04101" "K18076" "K18256" "K04102"
# 
# $`ko00627 Aminobenzoate degradation`
# [1] "K04101" "K04100" "K03862" "K03863" "K15066"
