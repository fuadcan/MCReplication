IDList  <-get(load("IDlist.rda"))
varVect <- IDList[[9]]
rhoVect <- IDList[[11]]

IDListM <- get(load("IDlistM.rda"))

varVectM<- IDListM[[3]]
rhoVectM<- IDListM[[4]]


ind2   <- c(6,3)
ind3   <- c(4,4,2)
ind4   <- c(5,5,3,2)
ind5.1 <- c(4,4,4,3,2)
ind5.2 <- c(6,6,4,4,3)
ind6   <- c(8,7,5,4,3,3)
ind7   <- c(7,7,7,6,6,4,2)
ind8   <- c(6,6,4,4,3,3,3,2)

k2 <- unlist(lapply(1:length(ind2), function(i) rep(IDListM[[1]][i],ind2[i])))
k2 <- c(k2, tail(IDListM[[1]],16)[0:(10 - length(k2))])
k3 <- unlist(lapply(1:length(ind3), function(i) rep(IDListM[[1]][i],ind3[i])))
k3 <- c(k3, tail(IDListM[[1]],16)[0:(10 - length(k3))])

k4 <- unlist(lapply(1:length(ind4), function(i) rep(IDListM[[1]][i],ind4[i])))
k4 <- c(k4, tail(IDListM[[1]],16)[0:(20 - length(k4))])
k5.1 <- unlist(lapply(1:length(ind5.1), function(i) rep(IDListM[[1]][i],ind5.1[i])))
k5.1 <- c(k5.1, tail(IDListM[[1]],16)[0:(20 - length(k5.1))])

k5.2 <- unlist(lapply(1:length(ind5.2), function(i) rep(IDListM[[1]][i],ind5.2[i])))
k5.2 <- c(k5.2, tail(IDListM[[1]],16)[0:(30 - length(k5.2))])
k6 <- unlist(lapply(1:length(ind6), function(i) rep(IDListM[[1]][i],ind6[i])))
k6 <- c(k6, tail(IDListM[[1]],16)[0:(30 - length(k6))])

k7 <- unlist(lapply(1:length(ind7), function(i) rep(IDListM[[1]][i],ind7[i])))
k7 <- c(k7, tail(IDListM[[1]],16)[0:(40 - length(k7))])
k8 <- unlist(lapply(1:length(ind8), function(i) rep(IDListM[[1]][i],ind8[i])))
k8 <- c(k8, tail(IDListM[[1]],16)[0:(40 - length(k8))])

# n=20; Tm=100; trial=2; clubSize=7; d=1;k=3;frho=.2

dataGen <- function(Tm,n,clubSize,frho,trial,noCons){
  
  # Script for generating data involving single club 
  
  # Inputs:
  
  #   Tm is the time step
  #   n is the # of the countries included in data
  #   clubSize is the # of the club members
  #   frho is the AR coefficient of the nonstationary factor's error terms
  #   trial is the quantity of the data that will be produced
  #   noCons is TRUE if logGDP series is not wanted to contain intercept
  
  # Outputs:
  # z <- Matrix of n countries with the amount of clubsize converging
  # ind <- indice vector states which countries are converging
  # gammas <- vector of country factor loadings 
  
  
  # Initial Parameters
  Tm <- Tm + 1000
  
  
  if(clubSize==10){gammaVect <- IDList[[7]]; ind <- IDList[[8]]} else 
    {gammaVect <- IDList[[clubSize-2]];ind<-IDList[[clubSize-1]]}
  
  gammaVect<- gammaVect[1:n]
  consVect <- IDList[[10]][1:n]
  
  quantt <- function(X){
    
    #  Generation of epsilon
    
    epsMat <- sapply(1:n, function(x) arima.sim(Tm, model = list(ar=rhoVect[x]), sd = sqrt((1-rhoVect[x]^2)*varVect[x])))
    
    
    #  Generation of the factor
    v_t <- arima.sim(n=Tm-1,model=list(ar=frho),sd=sqrt(1-frho^2))@.Data
    f_t <- 5 + c(0,cumsum(v_t))
    
    # Generation of y_it series
    
    y<-sapply(1:n,function(x) gammaVect[x]*f_t+epsMat[,x])
    
    if(noCons){ z <- y[1001:Tm,] } else {z <- t(consVect + t(y[1001:Tm,]))}
    
    
    
    return(z)
  }
  
  dats<- lapply(1:trial, function(x){ cat(paste0(x,"\n")); quantt(x)} )
  
  outFile <- list(dats,ind,gammaVect)
  outDir <- if(noCons){"Data/noCons/singleClub/"} else {"Data/withCons/singleClub/"}
  dir.create(outDir,recursive = TRUE)
  nocStr<- if(noCons){"-noCons.rda"} else {"-withCons.rda"}
  fileName <-  paste0(outDir,"zZ_",Tm-1000,"-",n,"-",clubSize,"-",frho,nocStr)
  
  save(outFile, file = fileName)
  
}

# Tm = 50;n = 30;k = 5;frho = 0.2;trial = 100;noCons = T

dataGenplus<-function(Tm,n,k,frho,trial,noCons){
  
  # Script for generating data involving multiple clubs
  
  # Inputs:
  
  #   Tm is the time step
  #   n is the # of the countries included in data
  #   k is the number of clubs 
  #   frho is the AR coefficient of the nonstationary factor's error term
  #   trial is the quantity of the data that will be produced
  #   noCons is TRUE if logGDP series is not wanted to contain intercept
  
  # Outputs:
  #   z <- Matrix of n countries form k convergence clubs
  #   ind <- indice vector states which countries are converging
  #   gammas <- vector of country factor loadings 
  
  Tm <- Tm+1000
  
  
  if(k==5 & n==20){indc <- ind5.1; gammaVect <- k5.1} else if(k==5 & n==30){indc<- ind5.2; gammaVect <- k5.2} else {indc <- get(paste0("ind",k)); gammaVect <- get(paste0("k",k))} 
  
  consVect <- c(IDListM[[2]][1:sum(indc)], tail(IDListM[[2]],16)[0:(n - sum(indc))])
  
  
  # Shuffled gammas
  gammac <- gammaVect[which(!duplicated(gammaVect))]
  
  shfs <- sample(1:n)
  
  gammaVect <- gammaVect[shfs]; consVect<-consVect[shfs]
  
  ind <- lapply(gammac, function(g) which(gammaVect==g))  
  
  
  quantt <- function(x){
    # Generation of epsilon
    
    epsMat <- sapply(1:n, function(x) arima.sim(Tm, model = list(ar=rhoVectM[x]), sd = sqrt((1-rhoVectM[x]^2)*varVectM[x])))
    
    epsMat<- epsMat[,shfs]
    
    #  Generation of the factor
    
    v_t <- arima.sim(n=Tm-1,model=list(ar=frho),sd=sqrt(1-frho^2))@.Data
    f_t <- 5 + c(0,cumsum(v_t))
    
    
    # Generation of y_it series
    if(noCons){y <- sapply(1:n, function(x) gammaVect[x]*f_t+epsMat[,x])} else
    {y<-sapply(1:n, function(x) consVect[x] + gammaVect[x]*f_t+epsMat[,x])}
    
    z <- y[1001:Tm,]
    
    return(z)
  }
  dat <- lapply(1:trial, function(x){cat(paste0(x,"\n")); quantt(x)})
  
  
  outFile <- list(dat,ind,gammaVect)
  
  outDir   <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  dir.create(outDir,recursive = TRUE)
  nocStr   <- if(noCons){"-noConsMltpois.rda"} else {"-withConsMltpois.rda"}
  fileName <-  paste0(outDir,"zZ_",Tm-1000,"-",n,"-",k,"-",frho,nocStr)
  
  save(outFile, file = fileName)
  
}