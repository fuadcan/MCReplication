# Tm<-50;n<-10;clsize<-5;frho=.2;ind=1

source("panelURPassFail.R")
source("mclapplyhack.R")
library("igraph")
library("urca")

options(cores=8)

mcAGK <- function(Tm,n,clsize,frho,noCons){

# Inputs:
  
  #   Tm     <- Length of time series
  #   n      <- number of country
  #   clsize <- number of club members
  #   frho   <- AR term of the coefficient of nonstationary factor's error term
  #   noCons <- TRUE if logGDP series doesn't contain intercept
  
# Outputs:
  
  # Report containing values to evaluate success
  # Gmml: Matrix containing gammas with detected convergent and nonconvergent country lists
  
  inDirAGK  <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  outDirAGK <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  cat("initializing variables\n")
  
  dataDir <- inDirAGK
  nocStr  <- if(noCons){"-noCons"} else {"-withCons"}
  zDir    <- paste0(dataDir,"zZ_",Tm,"-",n,"-",clsize,"-",frho,nocStr,".rda")
  
  zz     <- get(load(zDir))
  lenz   <- length(zz[[1]])
  list   <- zz[[2]]
  gammas <- zz[[3]]
  
  
  repForDatAGK <- function(ind){
    z      <- zz[[1]][[ind]]
    
    # Construction of panels containing all pairwise differences
    pPanel <- matrix(,nrow(z),(ncol(z)*(ncol(z)-1)/2))
    
    colnm  <- matrix(,1,(ncol(z)*(ncol(z)-1)/2))
    k=0
    
    for(i in 1:(ncol(z)-1)){
      
      for(j in (i+1):ncol(z)){
        k=k+1
        pPanel[,k]<-z[,i]-z[,j] 
        colnm[,k]<-paste(i,j,sep="-")
        
      }
      
    }
    colnames(pPanel) <- colnm
    
    cat("Analyzing\n")
    # Analyze and clustering 
    
    listOfpmats    <- panelURPassFail(pPanel,noCons)
      
    lenOfigrs <- 3
    
    # Pre-evaluation 
    listOfigraphs <- lapply(1:lenOfigrs, function(x) graph.adjacency(listOfpmats[,,x],"undirected",diag = F))
    lclqs         <- lapply(1:lenOfigrs ,function(x) sample(largest.cliques(listOfigraphs[[x]]),1))
    
    # This part may slowen the script but left for generalization...
    matrixialize <- function(lclq){
      len  <- length(lclq)
      lon  <- length(lclq[[1]])
      matrixed <- matrix(unlist(lclq),lon,len)
      matrixed <- sapply(1:len, function(x) matrixed[,x] <- sort(matrixed[,x]))
      matrixed <- t(matrixed)
      
      return(matrixed)
      
    }
    
    lgrps1 <- lapply(1:lenOfigrs, function(x) matrixialize(lclqs[[x]]) )
    maxInds <- sapply(1:length(lgrps1), function(y){ 
      maxInts <- sapply(1:dim(lgrps1[[y]])[1],function(x) length(intersect(lgrps1[[y]][x,],list)))
      indM    <- sample(which(maxInts==max(maxInts)),1); return(indM)})
    
    maxSuccess <- function(y,gr){max(sapply(1:dim(gr[[y]])[1],function(x) length(intersect(gr[[y]][x,],list))))}
    excOfMax   <- function(y,gr){maxInts <- sapply(1:dim(gr[[y]])[1],function(x) length(intersect(gr[[y]][x,],list)))
                                 indM    <- sample(which(maxInts==max(maxInts)),1)
                                 excofmax<- length(setdiff(gr[[y]][maxInds[y],],list))
                                 return(excofmax)}
    
    
    maxSccRate1 <- sapply(1: length(lgrps1),function(x) maxSuccess(x,lgrps1))
    excOfMax1   <- sapply(1: length(lgrps1),function(x) excOfMax(x,lgrps1))
      
    
    reportAGK      <- matrix(c(maxSccRate1,excOfMax1),,2)  
    
    gmml <- matrix("",3,n)
    for(i in 1:nrow(gmml)){gmml[i,][lgrps1[[i]][maxInds[i],]]<-"c"}
    
    
    return(list(reportAGK,gmml))
    
  }
  
  cat("Analyzing\n")     
  # Application for ntrial times 
  conc <- mclapply.hack(1:lenz, function(x){cat(paste0(Tm,"-",frho,"_",x,"\n")); repForDatAGK(x)})
#   conc <- lapply(1:lenz, function(x){cat(paste0(Tm,"-",frho,"_",x,"\n")); repForDatAGK(x)})
  
  # Consolidation of outputs
  lenOfigrs <- 3
  report <- lapply(1:lenz, function(x) t(conc[[x]][[1]]))
  report <- t(matrix(unlist(report),2,lenz*lenOfigrs))
  
  gmmls  <- lapply(1:lenz, function(x) t(conc[[x]][[2]]))
  gmmls  <-  rbind(gammas,t(matrix(unlist(gmmls),n)))
    
  clnamesR  <- rep(c("adf-0.01","adf-0.05","adf-0.10"),lenz)
  clnamesGL <- rep(c("adf-0.01","adf-0.05","adf-0.10"),lenz)
  rownames(report)<-clnamesR; rownames(gmmls)<- c("gmm",clnamesGL)
  
  dirName <- paste0(outDirAGK,"Results_",n,"-",clsize,"-AGK",nocStr)
  dir.create(dirName,recursive = T)
  
  
  reportname <- paste0(dirName,"/reportAGK-",Tm,"-",frho,".rda")
  gmmlname   <- paste0(dirName,"/gmmlAGK-",Tm,"-",frho,".rda")
  
  # Save
  save(report,file=reportname)
  save(gmmls,file=gmmlname)
 
  cat("Finished\n\n")
}