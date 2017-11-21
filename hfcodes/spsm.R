# setwd("C:/Users/fuadcan/Dropbox/FuatHarun")
library("parallel")
library("igraph")

options(core=8)

tstat <- function(y,model){
  y<-matrix(y,length(y),1)
  dy<-diff(y) # a vertical vector
  ny<-length(dy)
  shy<- matrix(y[-length(y),1],length(y)-1,1)
  
  t     <- if(model=="constant"){rep(1,ny)} else {t(matrix(c(rep(1,ny),1:ny),ny,2))}
  idM   <- matrix(0,ny,ny);diag(idM)<-1
  Mt    <- idM-t%*%solve(t(t)%*%t)%*%t(t) 
  varT  <- (t(dy)%*%Mt%*%dy)/ny
  tstat <- (t(dy)%*%Mt%*%shy)/(sqrt(varT*t(shy)%*%Mt%*%shy))
  
  return(tstat)
}

imp_crit <- read.table("imp_crit.csv",sep = ";",as.is = rep(T,5))
imp_crit <- imp_crit[-1,]
impcrit  <- apply(imp_crit[,2:5],2,as.numeric)


imp <- function(zp,model){
  Tler  <- as.numeric(imp_crit[,1][2:15])
  dzp   <- dim(zp); Tm<- dzp[1]; N<- dzp[2]
  int   <- sum(Tm>Tler)
  wght  <- (Tm-Tler[int])/(Tler[int+1]-Tler[int])
  wcrit <- (1-wght)*impcrit[int,]+wght*impcrit[int+1,]
  ts    <- apply(zp,2,function(y) tstat(y,model))
  tbar  <- mean(ts)
  zstat <- (sqrt(dzp[2])*(tbar-wcrit[3]))/(sqrt(wcrit[4])) 
  
  return(list(zstat,ts))
}

spsm<- function(z,model,alpha){
  Tler  <- as.numeric(imp_crit[,1][1:14])
  dz    <- dim(z)
  Tm    <- dz[1]; N<- dz[2]
  int   <- sum(Tm>Tler)
  wght  <- (Tm-Tler[int])/(Tler[int+1]-Tler[int])
  wcrit <- (1-wght)*impcrit[int,]+wght*impcrit[int+1,]
  
  crit   <- qnorm(alpha,0,1)
  indVec <- 1:N
  
  
  zt    <- z 
  ts    <- simplify2array(mclapply(1:N, function(y) tstat(z[,y],model)))
  tbar  <- mean(ts)
  zstat <- (sqrt(dzp[2])*(tbar-wcrit[3]))/(sqrt(wcrit[4])) 
  
  
  if(zstat>=crit){nstatvec <- indVec}
  lind=length(ts)
  cat(paste0(lind,"\n"))
  while((zstat < crit) & (lind>1)){
    tbar  <- mean(ts)
    zstat <- (sqrt(dzp[2])*(tbar-wcrit[3]))/(sqrt(wcrit[4])) 

    lind    <- lind-1
    cat(paste(lind,"\n",zstat,"\n"))
    
   if((zstat>=crit)|(lind==1)){nstatvec <- indVec} else {
    
    ind     <- which(ts==min(ts))
    indVec  <- indVec[-ind]
    zpt     <- zpt[,-ind]
    ts      <- ts[-ind]

   }
    
    
  }
    
  return(nstatvec)
  
}

  
  

pairsPanel<-function(z){
  ppan<-matrix(,nrow(z),(ncol(z)*(ncol(z)-1)/2))
  
  colnm<-matrix(,1,(ncol(z)*(ncol(z)-1)/2))
  k=0
  
  for(i in 1:(ncol(z)-1)){
    
    for(j in (i+1):ncol(z)){
      k=k+1
      ppan[,k]<-z[,i]-z[,j]
      colnm[,k]<-paste(i,j,sep="-")
      
    }
    
  }
  colnames(ppan)=colnm
  return(ppan)
}

