KEGG <- function(mz,mode="pos",add="H",tolerance="mDa",Q=10){
  #H <- 1.007825032
  H <- 1.0078
  mz.n <- round(mz,4)+H
  mz.p <- round(mz,4)-H
  if(tolerance=="mDa"){
    if(mode=="pos"){
      kcpd.list.p <- list()
      for(i in 1:length(mz.p)){
        x <- try(read.table(paste0("http://rest.kegg.jp/find/compound/",round(mz.p,1)[i],"/exact_mass")),silent=T)
        if (!inherits(x, 'try-error')){
          x <- x
          colnames(x) <- c("CPD","Exact_mass")
          x$dppm <- abs(((x$Exact_mass-mz.p[i])/mz.p[i])*10^6)
          x$delta <- abs(x$Exact_mass-mz.p[i])*10^3
          x<- cbind(QID=i,mz=mz[i],x)
          kcpd.list.p[[i]] <- x
          rm(x)
        } else x <- NA
      }
      
      kegg.pos <- do.call(rbind,kcpd.list.p)
      kegg.pos <- kegg.pos[!is.na(kegg.pos[,1]),]
      kegg.pos[which(kegg.pos$delta<=Q),]
    } else
      if(mode=="neg"){
        kcpd.list.p <- list()
        for(i in 1:length(mz.n)){
          x <- try(read.table(paste0("http://rest.kegg.jp/find/compound/",round(mz.n,1)[i],"/exact_mass")),silent=T)
          if (!inherits(x, 'try-error')){
            x <- x
            colnames(x) <- c("CPD","Exact_mass")
            x$dppm <- abs(((x$Exact_mass-mz.n[i])/mz.n[i])*10^6)
            x$delta <- abs(x$Exact_mass-mz.n[i])*10^3
            x<- cbind(QID=i,mz=mz[i],x)
            kcpd.list.p[[i]] <- x
            rm(x)
          } else x <- NA
        }
        
        kegg.neg <- do.call(rbind,kcpd.list.p)
        kegg.neg <- kegg.neg[!is.na(kegg.neg[,1]),]
        kegg.neg[which(kegg.neg$delta<=Q),]
      } else NULL
  } else
    if(tolerance=="ppm"){
      if(mode=="pos"){
        kcpd.list.p <- list()
        for(i in 1:length(mz.p)){
          x <- try(read.table(paste0("http://rest.kegg.jp/find/compound/",round(mz.p,1)[i],"/exact_mass")),silent=T)
          if (!inherits(x, 'try-error')){
            x <- x
            colnames(x) <- c("CPD","Exact_mass")
            x$dppm <- abs(((x$Exact_mass-mz.p[i])/mz.p[i])*10^6)
            x$delta <- abs(x$Exact_mass-mz.p[i])*10^3
            x<- cbind(QID=i,mz=mz[i],x)
            kcpd.list.p[[i]] <- x
            rm(x)
          } else x <- NA
        }
        
        kegg.pos <- do.call(rbind,kcpd.list.p)
        kegg.pos <- kegg.pos[!is.na(kegg.pos[,1]),]
        kegg.pos[which(kegg.pos$dppm<=Q),]
      } else
        if(mode=="neg"){
          kcpd.list.p <- list()
          for(i in 1:length(mz.n)){
            x <- try(read.table(paste0("http://rest.kegg.jp/find/compound/",round(mz.n,1)[i],"/exact_mass")),silent=T)
            if (!inherits(x, 'try-error')){
              x <- x
              colnames(x) <- c("CPD","Exact_mass")
              x$dppm <- abs(((x$Exact_mass-mz.n[i])/mz.n[i])*10^6)
              x$delta <- abs(x$Exact_mass-mz.n[i])*(10^3)
              x<- cbind(QID=i,mz=mz[i],x)
              kcpd.list.p[[i]] <- x
              rm(x)
            } else x <- NA
          }
          
          kegg.neg <- do.call(rbind,kcpd.list.p)
          kegg.neg <- kegg.neg[!is.na(kegg.neg[,1]),]
          kegg.neg[which(kegg.neg$dppm<=Q),]
        } else NULL
    } else NULL
}
