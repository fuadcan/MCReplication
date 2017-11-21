# n<-10;clsize=5;Tm=100;frho=.2;noCOns=F
library("rugarch")
hitRat <- function(maxScc,exc,clsize,n){
  UU <- maxScc; DU <- clsize-maxScc ; UD <- exc; DD <- n-clsize-exc
  hitR   <- (UU+DD)/(UU+DD+UD+DU)
  H      <- (UU)/(UU+DU)
  F      <- (UD)/(UD+DD)
  KS     <- H-F
  return(list(hitR,H,F,KS))
  
}

PTestHF<- function(dat){

  gmms    <- dat[1,]; dat<-dat[-1,]
  indlstA <- seq(1,nrow(dat),2); lstA <- dat[indlstA,]
  indlstR <- seq(2,nrow(dat),2); lstR <- dat[indlstR,]
  
  
  actua  <- (gmms==1)*1-(gmms!=1)*1
  forecA <- matrix(1,nrow(lstA),ncol(lstA)); forecA[lstA!="c"]<- -1
  forecR <- matrix(1,nrow(lstR),ncol(lstR)); forecR[lstR!="c"]<- -1   
  
  actua  <- rep(actua,nrow(dat)/2)
  dacA   <- DACTest(c(t(forecA)),actua,"PT",.05)$Stat

  dacR   <- DACTest(c(t(forecR)),actua,"PT",.05)$Stat
  dac    <- c(dacA,dacR)
    return(dac)
}

PTestAGK<- function(dat){
  
  gmms    <- dat[1,]
  dat     <- dat[-1,]
  

  indAdf  <- seq(1,nrow(dat),3)
  indErs  <- seq(2,nrow(dat),3)
  indKpss <- seq(3,nrow(dat),3)
  
  lists <- list()
  lists <- list(dat[indAdf,],dat[indErs,],dat[indKpss,])
  
  actua <- (gmms==1)*1-(gmms!=1)*1
  forec <- list()
  forec <- lapply(1:length(lists), function(i) forec[[i]] <- matrix(1,nrow(lists[[1]]),ncol(lists[[1]])))
  for(i in 1:length(lists)){forec[[i]][which(lists[[i]]=="")]<- -1}

  actuaVec <- rep(actua,nrow(dat)/length(lists))
  forecVecs <- lapply(1:3, function(i) c(t(forec[[i]])))
  
  
  ptRess <- sapply(1:3, function(x) DACTest(forecVecs[[x]],actuaVec,"PT",.05)$Stat)
  
  return(ptRess)
  
}


anlysofoutptAGK<- function(Tm,n,clsize,noCons){
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  dirName  <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  calcRep<-list()
  for(frho in c(2,6)){calcRep[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,"-AGK",nocStr,"/reportAGK-",Tm,"-",frho/10,".rda")))}
  
  gmml<-list()
  for(frho in c(2,6)){gmml[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,"-AGK",nocStr,"/gmmlAGK-",Tm,"-",frho/10,".rda")))}
  
  

  trials<-100
  TTm <- trials*3
  nT  <- 3
  
  anlysOutp <- function(calcRep,clsize,frho,gmml){
    
    ind <- frho*10
    
    
    indadf.01  <- seq(1,TTm,nT)
    indadf.05  <- seq(2,TTm,nT)
    indadf.1   <- seq(3,TTm,nT)
    
    adf.01  <- apply(calcRep[[ind]][indadf.01,],2,sum)
    adf.05  <- apply(calcRep[[ind]][indadf.05,],2,sum)
    adf.1   <- apply(calcRep[[ind]][indadf.1,],2,sum)
    
    adf.01prf  <- sum(calcRep[[ind]][indadf.01,1]==clsize & calcRep[[ind]][indadf.01,2]==0)
    adf.05prf  <- sum(calcRep[[ind]][indadf.05,1]==clsize & calcRep[[ind]][indadf.05,2]==0)
    adf.1prf   <- sum(calcRep[[ind]][indadf.1,1]==clsize & calcRep[[ind]][indadf.1,2]==0)
    
    prfs  <- c(adf.01prf,adf.05prf,adf.1prf)
    repss <- cbind(adf.01,adf.05,adf.1)
    
    totclsize <- clsize*trials; totn <- n*trials
 
    hitReps <- sapply(1:3, function(x) hitRat(repss[1,x],repss[2,x],totclsize,totn))
    rownames(hitReps) <- c("HitR", "H", "F", "KS")
    
    anlysPT <- PTestAGK(gmml[[ind]])
    rep     <- rbind(hitReps,anlysPT,prfs)
    
    return(rep)
  }
  
  tempRep <- lapply(c(2,6)/10, function(frho) anlysOutp(calcRep,clsize,frho,gmml))
  reps    <- lapply(1:6, function(x) rbind(tempRep[[1]][x,],tempRep[[2]][x,]))  
  
  
  
  return(reps)
  
}

anlysofoutptHF<- function(Tm,n,clsize,noCons){
  nocStr   <- if(noCons){"-noCons"} else {"-withCons"}
  dirName  <- if(noCons){"Output/noCons/singleClub/"} else {"Output/withCons/singleClub/"}
  
  calcRep<-list()
  for(frho in c(2,6)){calcRep[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_HF/reportHF-",Tm,"-",frho/10,".rda")))}
  
  gmml<-list()
  for(frho in c(2,6)){gmml[[frho]]<- get(load(paste0(dirName,"Results_",n,"-",clsize,nocStr,"_hf/gmmlHF-",Tm,"-",frho/10,".rda")))}
  
  trials <- 100
  TTm    <- trials*3
  nT     <- 3
  
  anlysOutp <- function(calcRep,clsize,frho,gmml){
    ind <- frho*10
    
    calcRepA <- calcRep[[ind]][,c(1,2)]  
    calcRepR <- calcRep[[ind]][,c(3,4)]
    
    HFA.01   <- calcRepA[rownames(calcRepA) == "crit01",]; HFR.01   <- calcRepR[rownames(calcRepR) == "crit01",]
    HFA.05   <- calcRepA[rownames(calcRepA) == "crit05",]; HFR.05   <- calcRepR[rownames(calcRepR) == "crit05",]
    HFA.1    <- calcRepA[rownames(calcRepA) == "crit1", ]; HFR.1    <- calcRepR[rownames(calcRepR) == "crit1", ]
    
    perfA    <- sapply(list(HFA.01,HFA.05,HFA.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    perfR    <- sapply(list(HFR.01,HFR.05,HFR.1), function(d) sum(d[,1]==clsize & d[,2]==0 ))
    
    hitRatA  <- sapply(list(HFA.01,HFA.05,HFA.1), function(d) apply(d,2,sum))
    hitRatR  <- sapply(list(HFR.01,HFR.05,HFR.1), function(d) apply(d,2,sum))
    
    
    repA <- apply(hitRatA, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*trials,n*trials)))
    repR <- apply(hitRatR, 2, function(hr) simplify2array(hitRat(hr[1],hr[2],clsize*trials,n*trials)))
    
    gmmls   <- lapply(c("01","05","1"), function(crit) gmml[[ind]][rownames(gmml[[ind]]) %in% c("gammas",paste0(c("abs","rel"), crit)), ])
    ptestAR <- sapply(gmmls, PTestHF)
    
    rownames(repA) <- rownames(repR) <- c("HitRat","H","F","KS")
    
    repA <- rbind(repA,ptestAR[1,],perfA)
    repR <- rbind(repR,ptestAR[2,],perfR)
    
    anlys <- cbind(repA,repR)
    colnames(anlys) <- unlist(lapply(c("A","R"), function(l) paste0(l,c(".01",".05",".1"))))
    
    return(anlys)
  }
  
  tempRepp <- lapply((c(2,6))/10, function(d) anlysOutp(calcRep,clsize,d,gmml))
  
  reps     <- lapply(1:6, function(x) rbind(tempRepp[[1]][x,],tempRepp[[2]][x,]))
  
  return(reps)
}

# TmVec <- c(50,75,100,200); nVec <- c(10,20,30,40); clsize <- c(3,5,7,10)

overallAnlys <- function(TmVec,n,clsize,noCons) {
  nlyAGK   <- lapply(TmVec, function(Tm) anlysofoutptAGK(Tm,n,clsize,noCons))
  nlyAGK   <- lapply(1:length(nlyAGK[[1]]), function(x) do.call(rbind,lapply(nlyAGK, function(n) n[[x]])))
  nlyHF    <- lapply(TmVec, function(Tm) anlysofoutptHF(Tm,n,clsize,noCons))
  nlyHF    <- lapply(1:length(nlyHF[[1]]), function(x)  do.call(rbind,lapply(nlyHF, function(n) n[[x]])))
  
  nocstr <- if(noCons){"A"} else {"R"}
  ovrNly   <- lapply(1:length(nlyAGK), function(x) cbind(nlyAGK[[x]],nlyHF[[x]][,grepl(nocstr,colnames(nlyHF[[x]]))]))
  ovrNly   <- lapply(1:length(nlyAGK), function(x) ovrNly[[x]][order(rep(c(.2,.6),length(TmVec))),])
  for(i in 1:length(ovrNly)){colnames(ovrNly[[i]])[1:3] <- c("adf.01","adf.05","adf.1")}
  
  return(ovrNly)
  
}


overall <- function(TmVec,nVec,clsize,noCons) {
  
  tempList <- lapply(nVec,  function(n) overallAnlys(TmVec,n,clsize,noCons))
  tempList <- lapply(1:length(tempList[[1]]), function(y) lapply(1:length(nVec), function(x) tempList[[x]][[y]]))
  tempList <- lapply(tempList, function(tlist) do.call(rbind,tlist))
  
  return(tempList)
}


overallRep<-function(){
  rep3T  <- overall(c(50,75,100,200),c(10,20,30,40), 3,T)
  rep5T  <- overall(c(50,75,100,200),c(10,20,30,40), 5,T)
  rep7T  <- overall(c(50,75,100,200),c(20,30,40)   , 7,T)
  rep10T <- overall(c(50,75,100,200),c(20,30,40)   , 10,T)
  rep3F  <- overall(c(50,75,100,200),c(10,20,30,40), 3,F)
  rep5F  <- overall(c(50,75,100,200),c(10,20,30,40), 5,F)
  rep7F  <- overall(c(50,75,100,200),c(20,30,40)   , 7,F)
  rep10F <- overall(c(50,75,100,200),c(20,30,40)   , 10,F)
  
  noConst   <- lapply(1:length(rep3T), function(x) rbind(rep3T[[x]],rep5T[[x]],rep7T[[x]],rep10T[[x]]))
  withConst <- lapply(1:length(rep3F), function(x) rbind(rep3F[[x]],rep5F[[x]],rep7F[[x]],rep10F[[x]]))
  
  dtType3  <- sapply(c(10,20,30,40), function(n) 
    sapply(c(.2,.6), function(frho) sapply(c(50,75,100,200), function(Tm) paste(3,n,frho,Tm,sep = "-"))))
  dtType5  <- sapply(c(10,20,30,40), function(n) 
    sapply(c(.2,.6), function(frho) sapply(c(50,75,100,200), function(Tm) paste(5,n,frho,Tm,sep = "-"))))
  dtType7  <- sapply(c(20,30,40), function(n) 
    sapply(c(.2,.6), function(frho) sapply(c(50,75,100,200), function(Tm) paste(7,n,frho,Tm,sep = "-"))))
  dtType10 <- sapply(c(20,30,40), function(n) 
    sapply(c(.2,.6), function(frho) sapply(c(50,75,100,200), function(Tm) paste(10,n,frho,Tm,sep = "-"))))
  dtType   <-c(c(dtType3),c(dtType5),c(dtType7),c(dtType10))
  
  for(i in 1:6){rownames(noConst[[i]]) <- dtType}; noConst   <- noConst[2:6]
  for(i in 1:6){rownames(withConst[[i]])<-dtType} ;withConst <- withConst[2:6]
  
  for(i in 1:5){colnames(noConst[[i]]) <- paste0(colnames(noConst[[i]]),"_",c("H","F","KS","PT","Perf")[i])}
  for(i in 1:5){colnames(withConst[[i]]) <- paste0(colnames(withConst[[i]]),"_",c("H","F","KS","PT","Perf")[i])}
  
  
  
  KSrepNoConst <- do.call(cbind,noConst[1:3]); KSrepWithConst <- do.call(cbind,withConst[1:3])
  PTSnoConst   <- cbind(noConst[[4]],rep(NA,nrow(noConst[[4]])),noConst[[5]]); 
  PTSwithConst <- cbind(withConst[[4]],rep(NA,nrow(withConst[[4]])),withConst[[5]])
  
  dir.create("Results",recursive = TRUE)
  write.csv(KSrepNoConst, file="Results/single_kupiersScoreNoCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(PTSnoConst, file="Results/single_PTandPerfSuccNoCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(KSrepWithConst, file="Results/single_kupiersScoreWithCons.csv",sep = ";",row.names = T, col.names = T)  
  write.csv(PTSwithConst, file="Results/single_PTandPerfSuccWithCons.csv",sep = ";",row.names = T, col.names = T)  
  
#   return(list(KSrepNoConst,PTSnoConst,KSrepWithConst,PTSwithConst))
}

