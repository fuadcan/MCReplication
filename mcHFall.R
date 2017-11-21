# importing datas
# Tm=75;n=30;clsize=5; alphaInt <- c(.01,.1); frho=.6; egeFact=1; ind=925
setwd("C:/Users/fuadcan/Documents/CADF_TablesN")
inDirHF  <- "C:/Users/fuadcan/Documents/Datas/noConsTrial/singleClub/"
outDirHF <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Outputs_noConsTrial/EGE_1/"

library("stringr")
mcHFall <- function(Tm,n,clsize,frho,egeFact){
  dirname  <- inDirHF
  
  
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",frho,"-e",egeFact,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz<-length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]][1:n]
  
  
  write.table(c(Tm,n),"C:/Users/fuadcan/Documents/CADF_TablesN/dims.csv",row.names = FALSE,col.names = FALSE)
  
  
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
  
    
    write.table(z,file = "C:/Users/fuadcan/Documents/CADF_TablesN/deneme1.csv",row.names = FALSE,col.names = FALSE)
    ############################# Hobijn & Franses ###########################
    tempHF01 <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main3.gss'), intern=TRUE,wait=TRUE)
    tempHF05 <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main05.gss'), intern=TRUE,wait=TRUE)
    tempHF1 <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main1.gss'), intern=TRUE,wait=TRUE)
    
    tempHF01 <- tempHF01[1:(which(tempHF01=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    tempHF05 <- tempHF05[1:(which(tempHF05=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    tempHF1 <- tempHF1[1:(which(tempHF1=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    
    ######################################################################
    
    tempHF<-tempHF01
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    ############################## Report ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(maxlAbs,list))
    maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(maxlRel,list))
    repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
    
    gmmlHF <- t(matrix(c(gammas,rep("",2*n)),n,3))
    gmmlHF[2,][maxlAbs] <- "c";  gmmlHF[3,][maxlRel] <- "c"
    

    
    repHF01<-repHF; gmmlHF01<-gmmlHF
    
    #######################################################
    
    tempHF<-tempHF05
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    ############################## Report ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(maxlAbs,list))
    maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(maxlRel,list))
    repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
    
    gmmlHF <- t(matrix(c(gammas,rep("",2*n)),n,3))
    gmmlHF[2,][maxlAbs] <- "c";  gmmlHF[3,][maxlRel] <- "c"
    
    repHF05<-repHF; gmmlHF05<-gmmlHF
    
    tempHF<-tempHF1
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    ############################## Report ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(setdiff(maxlAbs,list),0))
    maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(setdiff(maxlRel,list),0))
    repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
    
    gmmlHF <- t(matrix(c(gammas,rep("",2*n)),n,3))
    gmmlHF[2,][maxlAbs] <- "c";  gmmlHF[3,][maxlRel] <- "c"
    
    
    repHF1<-repHF; gmmlHF1<-gmmlHF
    return(list(repHF01,gmmlHF01,repHF05,gmmlHF05,repHF1,gmmlHF1))
    ############################## END REPORT ##############################
  }
  
  
  
  ############################# CONSOLIDATION ######################
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating\n")

  reportHF01 <- lapply(1:lenz, function(x) t(consConcs[[x]][[1]]))
  reportHF01 <- t(matrix(unlist(reportHF01),4,lenz))
  
  gmmlsHF01  <- lapply(1:lenz, function(x) t(consConcs[[x]][[2]]))
  gmmlsHF01  <- t(matrix(unlist(gmmlsHF01),n,3*lenz))
  
  reportHF05 <- lapply(1:lenz, function(x) t(consConcs[[x]][[3]]))
  reportHF05 <- t(matrix(unlist(reportHF05),4,lenz))
  
  gmmlsHF05  <- lapply(1:lenz, function(x) t(consConcs[[x]][[4]]))
  gmmlsHF05 <- t(matrix(unlist(gmmlsHF05),n,3*lenz))

  reportHF1 <- lapply(1:lenz, function(x) t(consConcs[[x]][[5]]))
  reportHF1 <- t(matrix(unlist(reportHF1),4,lenz))
  
  gmmlsHF1  <- lapply(1:lenz, function(x) t(consConcs[[x]][[6]]))
  gmmlsHF1 <- t(matrix(unlist(gmmlsHF1),n,3*lenz))
  
  
  savedir <- outDirHF
  outdir  <- paste0(savedir,"Results_",n,"-",clsize,"-e",egeFact,"_HFall/")
  filedir <- paste0(Tm,"-",frho,"-e",egeFact,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  save(reportHF01,file = paste0(outdir,"reportHF-",filedir))
  save(reportHF05,file = paste0(outdir,"reportHF05-",filedir))
  save(reportHF1,file = paste0(outdir,"reportHF1-",filedir))
  
  save(gmmlsHF01,file = paste0(outdir,"gmmlHF-",filedir))
  save(gmmlsHF05,file = paste0(outdir,"gmmlHF05-",filedir))
  save(gmmlsHF1,file = paste0(outdir,"gmmlHF1-",filedir))
  
#   return(reportHF)
}
