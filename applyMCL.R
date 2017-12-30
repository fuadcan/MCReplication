library("urca")
library("igraph")

# yearOrRegion <- "Penn"
# yearOrRegion <- 1930
applyMCL <- function(yearOrRegion){
# Reading data  
  if(is.numeric(yearOrRegion)){year      <- yearOrRegion
  filename  <- paste0("Application/madisonFrom-",year,".csv")
  z         <- read.table(filename,header = T,sep = ";")
  z         <- data.matrix(z)} else if(yearOrRegion=="Maddison"){
    
    filename  <- "Application/madisonFromMaxNA2-1950.csv"
    z         <- read.table(filename,header = T,sep = " ")
    z         <- data.matrix(z)
  } else if(yearOrRegion=="Penn"){
    z         <- read.csv("Application/penn.csv")
    z         <- data.matrix(z)} else if(yearOrRegion == "Europe"){
      z <- data.matrix(read.table("Application/mds_Europe-1950.csv",header = T,sep = ";"))} else {
        fname1 <- paste0("Application/mds_G7-1950.csv"); fname2 <- paste0("Application/mds_Europe-1950.csv"); fname3 <- paste0("Application/mds_S&P-1950.csv")
  z_g7         <- read.table(fname1,header = T,sep = ";"); z_g7 <- data.matrix(z_g7)
  z_eu         <- read.table(fname2,header = T,sep = ";"); z_eu <- data.matrix(z_eu)
  z_sp         <- read.table(fname3,header = T,sep = ";"); z_sp <- data.matrix(z_sp)
  
  
  if(yearOrRegion=="Europe+G7"|yearOrRegion=="G7+Europe")
  {z<- matrix(c(z_g7,z_eu),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_eu))} else
    if(yearOrRegion=="Europe+S&P"|yearOrRegion=="S&P+Europe")
    {z<- matrix(c(z_eu,z_sp),dim(z_eu)[1],); colnames(z)<- c(colnames(z_eu),colnames(z_sp))} else
      if(yearOrRegion=="G7+S&P"|yearOrRegion=="S&P+G7")
      {z<- matrix(c(z_g7,z_sp),dim(z_g7)[1],); colnames(z)<- c(colnames(z_g7),colnames(z_sp))} else {stop("Unknown Country List")}
  if(sum(duplicated(t(z)))!=0){z<- z[,-which(duplicated(t(z)))]}
    }
  
  z <- log(z)
  n <- ncol(z)
  cmbn <- combn(n,2)
  
  # Panel of Pairs
  pPanel <- sapply(1:ncol(cmbn), function(i) z[,cmbn[1,i]] - z[,cmbn[2,i]])
  cnames <- colnames(z)
  colnames(pPanel) <- sapply(1:ncol(cmbn), function(i) paste0(cnames[cmbn[1,i]],"-",cnames[cmbn[2,i]]))
  
  # Function for generating adjacency matrices
  panelURPassFail<-function(panel,noCons){
    
    if(noCons){tpA<-"none"} else {tpA<-"drift"} 
    
    # Calculating converging pairs
    testRsltA <- apply(panel, 2, function(s) {x <- s[!is.na(s)]; if(length(x) == 0){return(10000)} else {(ur.df(x,type = tpA)@teststat)[1]}})
    cvalsA    <- apply(panel, 2, function(x) {x <- x[!is.na(x)]; Tmm <- length(x); if(Tmm == 0){res <- -10000} else {res <- (ur.df(rnorm(Tmm),type = tpA)@cval)[1,2]}})
    
    passFailA <- cvalsA > testRsltA
    
    ncolZ <- ceiling(sqrt(ncol(panel)*2))
    
    # Converting to adjacency matrix
    pmGen<-function(pf){
      pairMat<-matrix(,ncolZ,ncolZ)
      pairMat[lower.tri(pairMat)] <- pf
      t(pairMat) -> pairMat
      pairMat[lower.tri(pairMat)] <- pf
      diag(pairMat)=0  
      
      return(pairMat)
    }
    
    pm <- pmGen(passFailA)
    
    
    return(pm) 
    
  }
  
  # Adjacency matrix of the 
  pm   <- panelURPassFail(pPanel,F)
  
  # Graphs
  gr   <- graph.adjacency(pm,"undirected",diag = F)
  
  # Maximal Clubs
  clubs <- lapply(maximal.cliques(gr), function(c) cnames[c])
  ord   <- order(sapply(clubs,length),decreasing = T)
  clubs <- matrix(sapply(clubs[ord], function(c) paste0(c,collapse = " - ")),,1)
  write.csv(clubs, file = paste0("Results/",yearOrRegion,"MCL.csv"))
  
  # Maximal Club Counts per sizes
  lcl    <- largest.cliques(gr)
  nlcl   <- if(class(lcl)=="list"){length(lcl[[1]])} else {length(lcl)}
  counts <- sapply(2:nlcl, function(n)  sum(sapply(maximal.cliques(gr), length) == n))
  
  # Membership Counts
  clubs     <- maximal.cliques(gr)
  clubs     <- clubs[sapply(clubs, length) > 1]
  overlaps  <- sapply(1:length(cnames), function(c) sum(sapply(clubs, function(cl) sum(c ==cl))))
  overlaps  <- cbind(cnames, overlaps)
  overlaps  <- overlaps[order(overlaps[,1]),]
  
  return(list(counts,overlaps))
}

ress1 <- lapply(c("Maddison","Penn","Europe","Europe+G7","Europe+S&P","G7+S&P"), applyMCL)
ress2 <- lapply(c(1930,1940), applyMCL)
ress <- c(ress2,ress1) # Results for all datasets

countsAll <- lapply(ress, function(r) r[[1]])
countsAll <- t(sapply(countsAll, function(cnts) c(cnts,rep(NA,max(sapply(countsAll,length))-length(cnts)))))
rownames(countsAll) <- c(1930,1940,"Maddison","Penn","Europe","Europe+G7","Europe+S&P","G7+S&P")
colnames(countsAll) <- paste0('# ',2:(ncol(countsAll)+1))

# Writing the reports
write.csv(countsAll,"Results/countsAllMCL.csv")
overlapsAll <- lapply(ress, function(r) r[[2]]) 
save(overlapsAll,file = "Results/overlapsAllMCL.rda")
