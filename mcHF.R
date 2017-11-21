# Tm=75;n=10;clsize=5; frho=.6; ind=1
#########################


library("stringr")
mcHF <- function(Tm,n,clsize,frho,noCons){

  inDirHF  <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  outDirHF <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  dirname  <- inDirHF
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  
  cat("Initializing variables\n")
  
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",frho,nocStr,".rda")
  zz       <- get(load(paste0(dirname,filename)))
  lenz<-length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]][1:n]
  
  #########################
  write.table(c(Tm,n),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  
  repForDat <- function(ind){
    
    z      <- zz[[1]][[ind]]
    
    
  write.table(z,file = "hfcodes/datt.csv",row.names = FALSE,col.names = FALSE)
  ########################################################

  cat("Analyzing")
  
  RtoGauss <- function(crit){
    tempHF <- shell(paste0("C:/gauss10/tgauss -b ",'RUN hfcodes\\main',crit,'.gss'), intern=TRUE,wait=TRUE)
    tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
    
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
    
    ############################## Pre-evaluation ##############################
    
    indAbs    <- sapply(1:length(absList), function(x) length(intersect(absList[[x]],list)))
    maxIndAbs <- sample(which(max(indAbs)==indAbs),1)
    maxlAbs   <- absList[[maxIndAbs]]
    indRel    <- sapply(1:length(relList), function(x) length(intersect(relList[[x]],list)))
    maxIndRel <- sample(which(max(indRel)==indRel),1)
    maxlRel   <- relList[[maxIndRel]]
    
    maxHFabs <- length(intersect(maxlAbs,list)); excHFabs<- length(setdiff(maxlAbs,list))
    maxHFrel <- length(intersect(maxlRel,list)); excHFrel<- length(setdiff(maxlRel,list))
    repHF <- matrix(c(maxHFabs,excHFabs,maxHFrel,excHFrel),1,4)
    
    gmmlHF <- t(matrix(rep("",2*n),n))
    gmmlHF[1,][maxlAbs] <- "c";  gmmlHF[2,][maxlRel] <- "c"
    return(list(repHF,gmmlHF))
    
  }
    
    temp <- lapply(c("01","05","1"), RtoGauss) 
    repHFs  <- do.call(rbind,lapply(temp, function(t) t[[1]]))
    gmmlHFs <- do.call(rbind,lapply(temp, function(t) t[[2]]))
    rownames(gmmlHFs) <- sapply(c("01","05","1"), function(c) paste0(c("abs","rel"),c))
    rownames(repHFs)  <- c("crit01","crit05","crit1")
  return(list(repHFs,gmmlHFs))
    
    
    ############################## END REPORT ##############################
  }
  
  
  
  ############################# CONSOLIDATION ######################
  
  consConcs <- lapply(1:lenz,function(x){cat(paste0(frho,"-",x,"\n")); repForDat(x)})
  cat("Consolidating\n")
  
  reportHF  <- do.call(rbind,lapply(consConcs, function(c) c[[1]]))
  gmmlsHF   <- do.call(rbind,lapply(consConcs, function(c) c[[2]]))
  gmmlsHF   <- rbind(gammas,gmmlsHF)
  
  
  savedir <- outDirHF
  outdir  <- paste0(savedir,"Results_",n,"-",clsize,nocStr,"_HF/")
  filedir <- paste0(Tm,"-",frho,".rda")
  
  
  dir.create(outdir,recursive = TRUE)
  cat("Saving\n")
  
  save(reportHF,file = paste0(outdir,"reportHF-",filedir))
  save(gmmlsHF,file = paste0(outdir,"gmmlHF-",filedir))
  
}
