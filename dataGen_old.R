IDList<-get(load("IDlist.rda"))
varVect<- IDList[[5]]
rhoVect<- IDList[[7]]

IDListM <- get(load("IDlistMulti.rda"))

varVectM<- IDListM[[3]]
rhoVectM<- IDListM[[4]]


k2   <- c(1:4, 8:11, 53:54)
k3   <- c(1:3, 8:10, 15:16, 53:54)
k4   <- c(1:5, 8:12, 15:18, 22:25, 53:54)
k5.1 <- c(1:4, 8:11, 15:18, 22:24, 29:31, 53:54)
k5.2 <- c(1:6, 8:13, 15:20, 22:26, 29:33, 53:54)
k6   <- c(1:5, 8:12, 15:19, 22:26, 29:32, 36:39, 53:54)
k7   <- c(1:6, 8:13, 15:20, 22:27, 29:34, 36:39, 43:46, 53:54)
k8   <- c(1:5, 8:12, 15:19, 22:26, 29:33, 36:40, 43:46, 49:52, 53:54)


# n=10; Tm=100; trial=2; clubSize=5; d=1;k=3;frho=.2


dataGenplus_old<-function(Tm,n,k,frho,trial,noCons){
  
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
  
  Tm <- Tm + 1000

  if(k==5 & n==20){indc<- k5.1} else if(k==5 & n==30){indc<- k5.2} else {indc<- get(paste0("k",k))} 
  gammaVect<- IDListM[[1]][indc]
  consVect <- IDListM[[2]][indc]
  

# Shuffled gammas
  gammac <- gammaVect[which(!duplicated(gammaVect[1:(n-2)]))]
  
  shfs <- sample(1:n)
  
  gammaVect <- gammaVect[shfs]; consVect<-consVect[shfs]
  
  ind <- lapply(1:length(gammac), function(x) which(gammaVect==gammac[x]))  


  quantt <- function(x){
  # Generation of epsilon
    
  epsMat <- sapply(indc, function(x) arima.sim(Tm, model = list(ar=rhoVectM[x]), sd = sqrt((1-rhoVectM[x]^2)*varVectM[x])))
    
  epsMat<- epsMat[,shfs]
  
  #  Generation of the factor
  
  v_t <- arima.sim(n=Tm-1,model=list(ar=frho),sd=sqrt(1-frho^2))@.Data
  f_t <- 5 + c(0,cumsum(v_t))
  
  
  # Generation of y_it series
  if(noCons){y <- sapply(1:n, function(x) gammaVect[x]*f_t+epsMat[,x])} else
    {y<-sapply(1:n, function(x) consVect[x] + gammaVect[x]*f_t+epsMat[,x])}
  
  z<-y[1001:Tm,]
  
  return(z)
  }
  dat <- lapply(1:trial, function(x){cat(paste0(x,"\n")); quantt(x)})


  outFile <- list(dat,ind,gammaVect)

  outDir  <- if(noCons){"Data/noCons/multiClub/"} else {"Data/withCons/multiClub/"}
  dir.create(outDir,recursive = TRUE)
  nocStr<- if(noCons){"-noConsMlt.rda"} else {"-withConsMlt.rda"}
  fileName <-  paste0(outDir,"zZ_",Tm-1000,"-",n,"-",k,"-",frho,nocStr)

  save(outFile, file = fileName)

}