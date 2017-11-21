# setwd("C:/Users/fuadcan/Dropbox/FuatHarun/CADF_TablesN")
# Tm=100;n=20;clsize=5; alphaInt <- c(.01,.1); d=1; egeFact=1; ind=1
mcCIPS <- function(Tm,n,clsize,d,alphaInt,egeFact){

  dirname  <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Datas/"
  filename <- paste0("zZ_",Tm,"-",n,"-",clsize,"-",d,"-e",egeFact,"_NO_DENS.rda")
  zz       <- get(load(paste0(dirname,filename)))

  lenz<-length(zz)
#   lenz<-10
  
    
  cat("initializing parameters\n")
  if(!is.numeric(alphaInt)){if(alphaInt=="all"){alphaInt<-c(.01,.1)}}
  
  alphaVec  <- c(0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1)
  alphaVec  <- alphaVec[alphaVec<=max(alphaInt) & alphaVec>=min(alphaInt)]
  alphaVec1 <- intersect(alphaVec,c(.01,.05,.1))
  
  critVals     <- read.table("critvals.csv",sep = ";")
  cadfcritvals <- list(); for(i in 1:3){cadfcritvals[[i]] <- critVals[10*(i-1)+2:9,2:9]}
  cipscritvals <- list(); for(i in 1:3){cipscritvals[[i]] <- critVals[10*(i-1)+2:9,12:19]}
  tmvec        <- as.numeric(as.numeric_version(critVals[2:9,1]))
  
  critvalsips  <- read.table("critvalsTbar.csv",sep=";")
  ipscritvals  <- list(); for(i in 1:3) ipscritvals[[i]] <- critvalsips[9*(i-1)+2:9,2:12]
  Tips         <- as.numeric(as.numeric_version(critvalsips[2:9,1]))
  Nips         <- as.numeric(as.numeric_version(critvalsips[1,2:12]))
  
  repForDatCIPS <- function(ind){
  
        
    zz     <- zz[[ind]]
    
    z      <- zz[[1]]
    z      <- z+10000
    list   <- zz[[2]]
    gammas <- zz[[3]]
    
  write.table(z,file = "deneme1.csv",row.names = FALSE,col.names = FALSE)
  
tabl <- shell(paste("C:/gauss10/tgauss.exe -b ",'RUN C:\\Users\\fuadcan\\Dropbox\\FuatHarun\\CADF_TablesN\\cadftbls.txt'), intern=TRUE,wait=TRUE)


rowN     <- which(tabl=="CADF_i(p) Statistics: Case 2 , WITH INTERCEPT")+3
tempCADF <- as.numeric(unlist(strsplit(tabl[(rowN):(rowN+n)],"  ")))
cadf     <- t(matrix(tempCADF[is.na(tempCADF)==FALSE],5,n))


rowNi <- which(tabl=="ADF_i(p) Statistics: Case 2 , WITH INTERCEPT")+3
tempADF <- as.numeric(unlist(strsplit(readLines("nanay.txt")[(rowNi):(rowNi+n)],"  ")))
adf     <- t(matrix(tempADF[is.na(tempADF)==FALSE],5,n))


tempPS <- shell(paste("C:/gauss10/tgauss -b",'RUN C:\\Users\\fuadcan\\Dropbox\\FuatHarun\\ps\\exm.pgm'), intern=TRUE,wait=TRUE)
indPS  <- seq(5,length(tempPS),13); indPS <- indPS[-tail(1:length(indPS),2)]
tempPS <- strsplit(tempPS[indPS]," ")
tempPS <- lapply(1:length(tempPS), function(x){tmp<-as.numeric(tempPS[[x]]); tmp <- tmp[is.na(tmp)==FALSE];return(tmp)})



spsm<- function(alpha){
  alphaInd <- sum(alpha>=alphaVec1)
  cipscritvals[alphaInd][[1]]<-t(cipscritvals[alphaInd][[1]])
  intT     <- sum(Tm>tmvec)
  intN     <- sum(n>tmvec)
  if(intT==length(tmvec)){wght<-1;intT<-intT-1} else {wght <- (Tm-tmvec[intT])/(tmvec[intT+1]-tmvec[intT])}
  wcrit  <- (1-wght)*cipscritvals[alphaInd][[1]][intT,]+wght*cipscritvals[alphaInd][[1]][intT+1,]
  if(intN==length(tmvec)){wghtN<-1;intN<-intN-1} else {wghtN <- (n-tmvec[intN])/(tmvec[intN+1]-tmvec[intN])}
  wcritN <- (1-wghtN)*wcrit[intN]+wghtN*wcrit[intN+1]
  
  
  indVec <- 1:n
    
  
  ts    <- cadf[,5]
  tbar  <- mean(ts)
  
  
  if(tbar>=wcritN){nstatvec <- indVec}
  lind=length(ts)
  cat(paste0(lind,"\n"))
  while(((tbar < wcritN) & (lind>1))==TRUE){
    tbar  <- mean(ts)
    
    lind    <- lind-1
    cat(paste(lind,"\n",tbar,"\n"))
  
      
      ind     <- sample(which(ts==min(ts)),1)
      indVec  <- indVec[-ind]
      ts      <- ts[-ind]
      
    
   nstatvec <- indVec 
    
  }
  
  return(nstatvec)
  
}


ips<- function(alpha){
  alphaInd <- sum(alpha>=alphaVec1)
  intT     <- sum(Tm>Tips)
  intN     <- sum(n>Nips)
  if(intT==length(Tips)){wght<-1;intT<-intT-1} else {wght <- (Tm-Tips[intT])/(Tips[intT+1]-Tips[intT])}
  wcrit  <- (1-wght)*ipscritvals[alphaInd][[1]][intT,]+wght*ipscritvals[alphaInd][[1]][intT+1,]
  if(intN==length(Nips)){wghtN<-1;intN<-intN-1} else {wghtN <- (n-Nips[intN])/(Nips[intN+1]-Nips[intN])}
  wcritN <- (1-wghtN)*wcrit[intN]+wghtN*wcrit[intN+1]
  
  
  indVec <- 1:n
  
  
  ts    <- adf[,5]
  tbar  <- mean(ts)
  
  
  if(tbar>=wcritN){nstatvec <- indVec}
  lind=length(ts)
  cat(paste0(lind,"\n"))
  while(((tbar < wcritN) & (lind>1))==TRUE){
    tbar  <- mean(ts)
    
    lind    <- lind-1
    cat(paste(lind,"\n",tbar,"\n"))
    
    
    ind     <- sample(which(ts==min(ts)),1)
    indVec  <- indVec[-ind]
    ts      <- ts[-ind]
    
    
    nstatvec <- indVec 
    
  }
  
  return(nstatvec)
  
}


############################### SPSM #####################################
# 
# conc  <- lapply(alphaVec1, function(x) spsm(x))
# conci <- lapply(alphaVec1, function(x) ips(x))
# 
# sccss      <-  function(cnc){length(intersect(cnc,list))}
# excss      <-  function(cnc){length(setdiff(cnc,list))}
# gammaOfExc <- function(cnc){gmexc <- gammas[setdiff(cnc,list)]; return(gmexc)}
# lofExc     <- function(cnc){lofexc<- rep("",n);lofexc[setdiff(cnc,list)] <- "exc";
#                         lofexc[is.na(lofexc)] <- "";   return(lofexc)}
# 
# 
# sccvec   <- sapply(1:length(alphaVec1), function(x) sccss(conc[[x]]))
# excvec   <- sapply(1:length(alphaVec1), function(x) excss(conc[[x]]))
# # gexcmatr <- t(sapply(1:length(alphaVec1), function(x) gammaOfExc(conc[[x]])))
# lexcvec  <- t(sapply(1:length(alphaVec1), function(x) lofExc(conc[[x]])))
# 
# sccveci   <- sapply(1:length(alphaVec1), function(x) sccss(conci[[x]]))
# excveci   <- sapply(1:length(alphaVec1), function(x) excss(conci[[x]]))
# # gexcmatri <- t(sapply(1:length(alphaVec1), function(x) gammaOfExc(conci[[x]])))
# lexcveci  <- t(sapply(1:length(alphaVec1), function(x) lofExc(conci[[x]])))



# cat("reporting\n")
# # Das report!
# report <- matrix(,length(alphaVec1),3)
# report[,1] <- alphaVec1
# report[,2] <- sccvec
# report[,3] <- excvec
# 
# gmmnlist <- matrix(0,length(alphaVec1)+1,n+1)
# gmmnlist[2:(1+length(alphaVec1)),1] <- alphaVec1
# gmmnlist[1,2:(n+1)] <- gammas
# gmmnlist[2:(1+length(alphaVec1)),2:(n+1)] <- lexcvec
# 
# reporti <- matrix(,length(alphaVec1),3)
# reporti[,1] <- alphaVec1
# reporti[,2] <- sccveci
# reporti[,3] <- excveci
# 
# gmmnlisti <- matrix(0,length(alphaVec1)+1,n+1)
# gmmnlisti[2:(1+length(alphaVec1)),1] <- alphaVec1
# gmmnlisti[1,2:(n+1)] <- gammas
# gmmnlisti[2:(1+length(alphaVec1)),2:(n+1)] <- lexcveci
# 
# colnames(report) <- c("alpha","Success", "Exceeders")
# 
# cat("finished\n\n")
# #   colnames(dasReportante) <- c("Data Type", "Test Type", "alpha", "Avg Success", "Max Success", "Exceeders of Max")
# ############################### PS #####################################
# 

maxInds  <- sapply(1:length(tempPS), function(x) length(intersect(tempPS[[x]],list)))
maxInd   <- sample(which(max(maxInds)==maxInds),1)
maxI     <- intersect(tempPS[[maxInd]],list)
maxsccps <- maxInds[maxInd]
excvecps <- setdiff(tempPS[[maxInd]],list)
excps    <- length(excvecps)
gexcps   <- gammas[excvecps]
lexcps   <- rep("",n); lexcps[maxI] <- "c"


reportPS <- matrix(,1,2)
reportPS[,1]<- maxsccps
reportPS[,2]<- excps

gmmnlistPS <- matrix(,2,n)
gmmnlistPS[1,] <- gammas
gmmnlistPS[2,] <- lexcps


return(list(reportPS,gmmnlistPS))
}


consConcs <- lapply(1:lenz,function(x){cat(paste0(x,"\n")); repForDatCIPS(x)})
lconcon <- length(consConcs)

# repspsmTMP <- sapply(1:lconcon, function(x) t(consConcs[[x]][[1]]))
# repspsmTMP <- t(matrix(repspsmTMP,3,3*lconcon))
# 
# gmmlspsmTMP <- sapply(1:lconcon, function(x) t(consConcs[[x]][[2]]))
# gmmlspsmTMP <- t(matrix(gmmlspsmTMP,n+1,4*lconcon))

reppsTMP <- t(sapply(1:lconcon, function(x) consConcs[[x]][[1]]))

gmmlpsTMP  <- lapply(1:lconcon, function(x) t(consConcs[[x]][[2]]))
gmmlpsTMP  <- t(matrix(unlist(gmmlpsTMP),n,2*lconcon))

# repipsTMP <- sapply(1:lconcon, function(x) t(consConcs[[x]][[5]]))
# repipsTMP <- t(matrix(repipsTMP,3,3*lconcon))
# 
# gmmlipsTMP <- sapply(1:lconcon, function(x) t(consConcs[[x]][[6]]))
# gmmlipsTMP <- t(matrix(gmmlipsTMP,n+1,4*lconcon))
# 
# 
# repspsm     <- matrix(,dim(repspsmTMP)[1],dim(repspsmTMP)[2]+1) 
# repspsm[,1] <- sort(rep(1:lconcon,length(alphaVec1)))
# repspsm[,2:(dim(repspsm)[2])] <- repspsmTMP
# 
# gmmlspsm  <- matrix(,dim(gmmlspsmTMP)[1],dim(gmmlspsmTMP)[2]+1) 
# gmmlspsm[,1]  <- sort(rep(1:lconcon,length(alphaVec1)+1))
# gmmlspsm[,2:(dim(gmmlspsm)[2])] <- gmmlspsmTMP
# 
# repips     <- matrix(,dim(repipsTMP)[1],dim(repipsTMP)[2]+1) 
# repips[,1] <- sort(rep(1:lconcon,length(alphaVec1)))
# repips[,2:(dim(repips)[2])] <- repipsTMP
# 
# gmmlips  <- matrix(,dim(gmmlipsTMP)[1],dim(gmmlipsTMP)[2]+1) 
# gmmlips[,1]  <- sort(rep(1:lconcon,length(alphaVec1)+1))
# gmmlips[,2:(dim(gmmlips)[2])] <- gmmlipsTMP


repps     <- matrix(,dim(reppsTMP)[1],dim(reppsTMP)[2]+1) 
repps[,1] <- 1:lconcon
repps[,2:(dim(repps)[2])] <- reppsTMP


gmmlps      <- matrix(,dim(gmmlpsTMP)[1],dim(gmmlpsTMP)[2]+1) 
gmmlps[,1]  <- sort(rep(1:lconcon,2))
gmmlps[,2:(dim(gmmlps)[2])] <- gmmlpsTMP

colnames(repps) <- c("trial","Success","# Exceeders")
# colnames(repspsm) <- c("trial","alpha","Success","# Exceeders")

savedir <- "C:/Users/fuadcan/Dropbox/Ege-Proje/maxClique/hf/Outputs/EGE_1/"
outdir  <- paste0(savedir,"Results_",n,"-",clsize,"-e",egeFact,"_PS/")
filedir <- paste0(Tm,"-",d,"-e",egeFact,".rda")


dir.create(outdir,recursive = TRUE)

# save(repspsm,file = paste0(outdir,"reportSPSM-",filedir))
# save(gmmlspsm,file = paste0(outdir,"gmmlSPSM-",filedir))
save(repps,file = paste0(outdir,"reportPS-",filedir))
save(gmmlps,file = paste0(outdir,"gmmlPS-",filedir))
# save(repips,file = paste0(outdir,"reportIPS-",filedir))
# save(gmmlips,file = paste0(outdir,"gmmlIPS-",filedir))

return(list(repps,gmmlps))
}