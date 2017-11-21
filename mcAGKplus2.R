# setwd("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes")
# Tm<-50;n<-10;k<-2;frho<-.6;egeFact<-1;ind=1;

source("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes/panelURPassFailMC2.R")
source("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes/mclapplyhack.R")
source("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes/clusterAGK.R")
library("igraph")
library("fUnitRoots")

inDirM  <- "C:/Users/fuadcan/Documents/Datas/noConsTrial/multiClub/"
outDirM <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Outputs_noConsTrial/"
options(cores=8)

mcAGKplus2 <- function(Tm,n,k,frho,egeFact){
cat("initializing parameters\n")
dataDir<- inDirM
zDir <- paste0(dataDir,"zZ_",Tm,"-",n,"-",k,"-",frho,"-e",egeFact,".rda")
zz <- get(load(zDir))

lenz<-length(zz[[1]])
list   <- zz[[2]]
gammas <- zz[[3]]

fsts<- sapply(list, function(x) x[1])  
list<- list[order(fsts)]


repForDatAGK <- function(ind){
  
  
  z      <- zz[[1]][[ind]]
  
  
  pPanel<-matrix(,nrow(z),(ncol(z)*(ncol(z)-1)/2))
  
  colnm<-matrix(,1,(ncol(z)*(ncol(z)-1)/2))
  k=0
  
  for(i in 1:(ncol(z)-1)){
    
    for(j in (i+1):ncol(z)){
      k=k+1
      pPanel[,k]<-z[,i]-z[,j]#!!!!!!!!!!!!?NEML?!!!!!!!!!!!! 
      colnm[,k]<-paste(i,j,sep="-")
      
    }
    
  }
  colnames(pPanel)<-colnm
  
  cat("Constructing Adjacency Matrices\n") 
  listOfpmats <- panelURPassFailMC2(pPanel)
  cat("Constructed\n") 
  
  lp <- dim(listOfpmats[[1]])[3]
  
  listOftvals <- listOfpmats[[2]]
  
  cat("Clustering\n") 
  clists <- lapply(1:lp, function(x) clusterAGK(listOfpmats[[1]][,,x]))
  
  cat("Evaluating\n") 
  
  listCode <- unlist(lapply(list, function(x) c(x,0)))
  sccAGK<- suppressWarnings(sapply(clists, function(cl) mean(unlist(lapply(cl, function(x) c(x,0)))==listCode)==1))
  
  repAGK  <- matrix(c(sccAGK,sccAGK==1),,2)
  gmmlAGK <- t(matrix(c(gammas,rep("",dim(repAGK)[1]*n)),n))
  for(i in 1:dim(repAGK)[1]){
    for(j in 1:length(clists[[i]])){gmmlAGK[1+i,][clists[[i]][[j]]] <- paste0("c",j)}} 
  
  return(list(repAGK,gmmlAGK,listOfpmats,listOftvals))
  
}

cat("Consolidating\n") 

conc <- mclapply.hack(1:lenz, function(x){cat(paste0(Tm,"-",frho,"_",x,"\n")); repForDatAGK(x)})


lenOfigrs <- 9
report <- lapply(1:lenz, function(x) t(conc[[x]][[1]]))
report <- t(matrix(unlist(report),2,lenz*lenOfigrs))

gmmls <- lapply(1:lenz, function(x) t(conc[[x]][[2]]))
gmmls <-  t(matrix(unlist(gmmls),n,lenz*(lenOfigrs+1)))

pmats <- simplify2array(lapply(1:lenz, function(x) conc[[x]][[3]]))
tvals <- simplify2array(lapply(1:lenz, function(x) conc[[x]][[4]]))

alphaRepVec <- rep(c(.01,.05,.1),3); typeRepVec  <- sort(rep(c("adf","gls","kpss"),3))
clnamesR  <- paste0(typeRepVec,alphaRepVec); clnamesR <- rep(clnamesR,lenz)
clnamesGL <- paste0(typeRepVec,alphaRepVec); clnamesGL<- rep(c("gmm",clnamesGL),lenz)
rownames(report)<-clnamesR; rownames(gmmls)<-clnamesGL
colnames(report)<- c("Rate","Success?")

dirName <- paste0(outDirM,"EGE_",egeFact,"_plus/Results_",n,"-",k,"-e",egeFact,"AGK")
dir.create(dirName,recursive = TRUE)

cat("Saving\n") 
reportname <- paste0(dirName,"/reportAGK-",Tm,"-",frho,"-e",egeFact,".rda")
gmmlname  <- paste0(dirName,"/gmmlAGK-",Tm,"-",frho,"-e",egeFact,".rda")
pmatname  <- paste0(dirName,"/pmatAGK-",Tm,"-",frho,"-e",egeFact,".rda")
tvalname  <- paste0(dirName,"/tvalAGK-",Tm,"-",frho,"-e",egeFact,".rda")

save(report,file=reportname)
save(gmmls,file=gmmlname)
save(pmats,file=pmatname)
save(tvals,file=tvalname)
cat("Finished\n\n") 


}








