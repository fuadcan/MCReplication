# Tm<-50;n<-10;k<-2;frho<-.2;ind=1; noCons=T; nopois <- T

source("panelURPassFailMCA.R")
source("mclapplyhack.R")
source("clusterAGK.R")
library("igraph")
library("urca")


options(cores=8)

mcAplus <- function(Tm,n,k,frho,noCons,nopois){
  
  inDirM  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  outDirM <- if(noCons){"Output/noCons/multiClub/"} else {"Output/withCons/multiClub/"}
  
  cat("initializing parameters\n")
  dataDir <- inDirM
  nocStr  <- if(noCons){"-noConsMlt"} else {"-withConsMlt"}
  poistr  <- if(nopois){NULL} else {"pois"}
  
  zDir <- paste0(dataDir,"zZ_",Tm,"-",n,"-",k,"-",frho,nocStr,poistr,".rda")
  zz   <- get(load(zDir))
  
  lenz   <- length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]]
  
  list <- list[sapply(list, length)!=1]
  fsts <- sapply(list, function(x) x[1])  
  list <- list[order(fsts)]
  
  
  repForDatAGK <- function(ind){
    
    z      <- zz[[1]][[ind]]
    
    # Construction of panels containing all pairwise differences
    
    pPanel<-matrix(,nrow(z),(ncol(z)*(ncol(z)-1)/2))
    colnm<-matrix(,1,(ncol(z)*(ncol(z)-1)/2))
    k=0
    
    for(i in 1:(ncol(z)-1)){
      
      for(j in (i+1):ncol(z)){
        k=k+1
        pPanel[,k]<-z[,i]-z[,j] 
        colnm[,k]<-paste(i,j,sep="-")
        
      }
      
    }
    colnames(pPanel) <- colnm
    
    # Analyze and clustering 
    listOfpmats <- panelURPassFailMCA(pPanel,noCons)

    clists <- lapply(listOfpmats, clusterAGK)
    # Pre-evaluation 
    
    listCode <- unlist(lapply(list, function(x) c(x,0)))
    sccAGK   <- sapply(clists, function(clist)suppressWarnings(mean(unlist(lapply(clist, function(x) c(x,0)))==listCode)==1))
    clists
    
    gmmlAGK  <- sapply(1:length(clists), function(i){temp <- rep("",n);
      for(j in 1:length(clists[[i]])){temp[clists[[i]][[j]]] <- paste0("c",j)}; return(temp)})
    
    colnames(gmmlAGK) <- paste0("crt-",c(.01,.05,.1))
    return(list(sccAGK,gmmlAGK))
    
  } 
  
  cat("Analyzing\n") 
  # Application for ntrial times 
  conc <- mclapply.hack(1:lenz, function(x){cat(paste0(Tm,"-",frho,"_",x,"\n")); repForDatAGK(x)})
  
  
  # Consolidation of outputs
  report <- sapply(conc, function(c) c[[1]])
  
  gmmls  <- do.call(rbind,lapply(conc, function(c) t(c[[2]])))
  gmmls  <- rbind(gammas,gmmls)
  
  # Save
  dirName <- paste0(outDirM,"Results_",n,"-",k,"-A",nocStr)
  dir.create(dirName,recursive = TRUE)
  
  reportname <- paste0(dirName,"/reportA-",Tm,"-",frho,poistr,".rda")
  gmmlname   <- paste0(dirName,"/gmmlA-",Tm,"-",frho,poistr,".rda")
  
  save(report, file = reportname)
  save(gmmls,  file = gmmlname)
  cat("Finished\n\n") 
  
  
}

