library("igraph")
library("urca")


panelURPassFailMCA<-function(z,noCons){
  #   Finds the pairs for given type of test (adf, kpss or ADF-GLS) and constructs the igraph
  #   
  #   Inputs: 
  #     z->Panel dataset
  #   
  #   
  #   Outputs:
  #     pairMat-> Matrix with binaries identifies the pairs
  #   
  
  Tmm <- nrow(z)
  
  if(noCons){tpA<-"none"; tpK<-NULL; tpG<-NULL} else {tpA<-"drift"; tpK<-"mu"; tpG<-"constant"} 
  testRsltA <- apply(z,2,function(x) (ur.df(x,type = tpA)@teststat)[1])
  cvalsA    <- (ur.df(rnorm(Tmm),type = tpA)@cval)[1,]
  
  passFailA <- sapply(cvalsA, function(c) c > testRsltA)

  ncolZ <- ceiling(sqrt(ncol(z)*2))
  
  
  
  pmGen <- function(pf){
    
    pairMat<-matrix(,ncolZ,ncolZ)
    pairMat[lower.tri(pairMat)] <- pf
    t(pairMat) -> pairMat
    pairMat[lower.tri(pairMat)] <- pf
    diag(pairMat) <- 0  
    
    return(pairMat)
  }
  
  pairMat <- lapply(1:ncol(passFailA), function(x) pmGen(passFailA[,x]))
  
  return(pairMat) 
  
}