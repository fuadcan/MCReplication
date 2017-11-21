# setwd("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes")
setwd("C:/Users/fuadcan/Documents/CADF_TablesN")
# Tm=100;n=10;k=3; frho=.6; egeFact=1; ind=51
library("stringr")

indirHFall  <- "C:/Users/fuadcan/Documents/Datas/thirdTrial/multiClub/"
outdirHFall <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Outputs_thirdTrial/EGE_1_plus/"

mcHFallplus <- function(Tm,n,k,frho,egeFact){
  
  
  
  dirname  <- indirHFall
  #   dirname  <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Datas/"
  
  filename <- paste0("zZ_",Tm,"-",n,"-",k,"-",frho,"-e",egeFact,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz<-length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]]

  fsts<- sapply(list, function(x) x[1])  
  list<- list[order(fsts)]
  
  
  write.table(c(Tm,n),"dims.csv",row.names = FALSE,col.names = FALSE)
  
  cat("initializing parameters\n")
  
  
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]

    
    
    write.table(z,file = "deneme1.csv",row.names = FALSE,col.names = FALSE)
    
    cat("Analyzing\n")
    
    ############################# Hobijn & Franses ###########################
    tempHF <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main05.gss'), intern=TRUE,wait=TRUE)
    tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
    absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    absList05<- absList; relList05<- relList

    tempHF <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main1.gss'), intern=TRUE,wait=TRUE)
    tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
    absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    absList1<- absList; relList1<- relList

    
    tempHF <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main3.gss'), intern=TRUE,wait=TRUE)
    tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
    aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
    aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
    absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
    
    rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
    rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
    rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
    relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
    
    fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
    absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
    
    if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
    if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
    
    absList01<- absList; relList01<- relList
    
    ############################## Evaluation ##############################
    cat("Evaluation\n")
    listCode <- unlist(lapply(list, function(x) c(x,0)))
    
    sccHFabs01<- suppressWarnings(mean(unlist(lapply(absList01, function(x) c(sort(x),0)))==listCode)==1)
    sccHFrel01<- suppressWarnings(mean(unlist(lapply(relList01, function(x) c(sort(x),0)))==listCode)==1)
    
    repHF01<- matrix(c(sccHFabs01,sccHFabs01==1,sccHFrel01,sccHFrel01==1),1)
    gmmlHF01 <- t(matrix(c(gammas,rep("",2*n)),n,3))
    for(i in 1:length(absList01)){gmmlHF01[2,][absList01[[i]]] <- paste0("c",i)} 
    for(i in 1:length(relList01)){gmmlHF01[3,][relList01[[i]]] <- paste0("c",i)}
    
    sccHFabs05<- suppressWarnings(mean(unlist(lapply(absList05, function(x) c(sort(x),0)))==listCode)==1)
    sccHFrel05<- suppressWarnings(mean(unlist(lapply(relList05, function(x) c(sort(x),0)))==listCode)==1)
    
    repHF05<- matrix(c(sccHFabs05,sccHFabs05==1,sccHFrel05,sccHFrel05==1),1)
    gmmlHF05 <- t(matrix(c(gammas,rep("",2*n)),n,3))
    for(i in 1:length(absList05)){gmmlHF05[2,][absList05[[i]]] <- paste0("c",i)} 
    for(i in 1:length(relList05)){gmmlHF05[3,][relList05[[i]]] <- paste0("c",i)}
    

    sccHFabs1<- suppressWarnings(mean(unlist(lapply(absList1, function(x) c(sort(x),0)))==listCode)==1)
    sccHFrel1<- suppressWarnings(mean(unlist(lapply(relList1, function(x) c(sort(x),0)))==listCode)==1)
    
    repHF1<- matrix(c(sccHFabs1,sccHFabs1==1,sccHFrel1,sccHFrel1==1),1)
    gmmlHF1 <- t(matrix(c(gammas,rep("",2*n)),n,3))
    for(i in 1:length(absList1)){gmmlHF1[2,][absList1[[i]]] <- paste0("c",i)} 
    for(i in 1:length(relList1)){gmmlHF1[3,][relList1[[i]]] <- paste0("c",i)}
    
    
    return(list(repHF05,gmmlHF05,repHF1,gmmlHF1,repHF01,gmmlHF01))
    ############################## END REPORT ##############################
  }
  
  ############################# CONSOLIDATION ######################
  cat("Consolidating\n")
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating\n")
  
  reportHF05 <- lapply(1:lenz, function(x) t(consConcs[[x]][[1]]))
  reportHF05 <- t(matrix(unlist(reportHF05),4,lenz))
  
  gmmlsHF05  <- lapply(1:lenz, function(x) t(consConcs[[x]][[2]]))
  gmmlsHF05 <- t(matrix(unlist(gmmlsHF05),n,3*lenz))
  
  reportHF1 <- lapply(1:lenz, function(x) t(consConcs[[x]][[3]]))
  reportHF1 <- t(matrix(unlist(reportHF1),4,lenz))
  
  gmmlsHF1  <- lapply(1:lenz, function(x) t(consConcs[[x]][[4]]))
  gmmlsHF1 <- t(matrix(unlist(gmmlsHF1),n,3*lenz))
  
  reportHF01 <- lapply(1:lenz, function(x) t(consConcs[[x]][[5]]))
  reportHF01 <- t(matrix(unlist(reportHF01),4,lenz))
  
  gmmlsHF01  <- lapply(1:lenz, function(x) t(consConcs[[x]][[6]]))
  gmmlsHF01 <- t(matrix(unlist(gmmlsHF01),n,3*lenz))
  
#   return(list(reportHF05,reportHF1))
  
  savedir <- outdirHFall
  outdir  <- paste0(savedir,"Results_",n,"-",k,"-e",egeFact,"_HFall/")
  filedir <- paste0(Tm,"-",frho,"-e",egeFact,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  save(reportHF05,file = paste0(outdir,"reportHF05-",filedir))
  save(reportHF1,file = paste0(outdir,"reportHF1-",filedir))
  save(gmmlsHF05,file = paste0(outdir,"gmmlHF05-",filedir))
  save(gmmlsHF1,file = paste0(outdir,"gmmlHF1-",filedir))

  save(reportHF01,file = paste0(outdir,"reportHF-",filedir))
  save(gmmlsHF01,file = paste0(outdir,"gmmlHF-",filedir))

  cat("Finished\n")
  
}