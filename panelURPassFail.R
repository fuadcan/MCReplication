library("igraph")
library("urca")


panelURPassFail<-function(z,noCons){
#   Finds the pairs for given type of test (adf, kpss or ADF-GLS) and constructs the igraph
#   
#   Inputs: 
#     z->Panel dataset
#     noCons -> TRUE if lonGDP series doesn't contain intercept
#   
#   Outputs:
#     pairMat-> Matrix with binaries identifies the pairs
#     gr-> the graph constructed by igraph
#   

Tmm <- nrow(z)
if(noCons){tpA<-"none"; tpK<-NULL; tpG<-NULL} else {tpA<-"drift"; tpK<-"mu"; tpG<-"constant"} 
testRsltA <- apply(z,2,function(x) (ur.df(x,type = tpA)@teststat)[1])
cvalsA    <- (ur.df(rnorm(Tmm),type = tpA)@cval)[1,]
# 
# testRsltK <- apply(z,2,function(x) ur.kpss(x,type=tpK)@teststat)
# cvalsK    <- ur.kpss(rnorm(Tmm),type=tpK)@cval[-3]
# 
# testRsltE <- apply(z,2,function(x) ur.ers(x,model=tpG)@teststat)
# cvalsE    <- ur.ers(rnorm(Tmm),model=tpG)@cval[1,]

passFailA <- sapply(cvalsA, function(c) c > testRsltA)
# passFailK <- sapply(cvalsK, function(c) c > testRsltK)
# passFailE <- sapply(cvalsE, function(c) c > testRsltE)

ncolZ <- ceiling(sqrt(ncol(z)*2))

pmGen <- function(pf){
  
  pairMat<-matrix(,ncolZ,ncolZ)
  pairMat[lower.tri(pairMat)] <- pf
  t(pairMat) -> pairMat
  pairMat[lower.tri(pairMat)] <- pf
  diag(pairMat) <- 0  
  
  return(pairMat)
}

pmsA <- lapply(1:ncol(passFailA), function(x) pmGen(passFailA[,x]))
# pmsE <- lapply(1:ncol(passFailE), function(x) pmGen(passFailE[,x]))
# pmsK <- lapply(1:ncol(passFailK), function(x) pmGen(passFailK[,x]))

# pms  <- simplify2array(c(pmsA,pmsE,pmsK))
pms  <- simplify2array(pmsA)


return(pms) 

}