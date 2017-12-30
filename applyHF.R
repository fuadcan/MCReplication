library("urca")
library("igraph")

# yearOrRegion <- "Europe"
applyHF <- function(yearOrRegion){
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
  cnames <- colnames(z)
  
  # Saving data to disc for Gauss to read
  write.table(c(nrow(z),ncol(z)),"hfcodes/dims.csv",row.names = FALSE,col.names = FALSE)
  write.table(z,file = "hfcodes/datt.csv",row.names = FALSE,col.names = FALSE)
  
    cat("Analyzing\n")
    
      ########################################################
      tempHF <- shell(paste0("C:/gauss6.0/tgauss -b ",'RUN hfcodes\\main05.gss'), intern=TRUE,wait=TRUE)
      tempHF <- tempHF[1:(which(grepl("GAUSS",tempHF))[1]-2)]
      
      # Processing outputs for absolute and relative convergence 
      aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
      aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
      aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
      absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]
      
      rCrude<- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
      rCrude<-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
      rCrude<-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
      relList<- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]
      
      fstsABS <- sapply(absList, function(x) x[1]); fstsREL <- sapply(relList, function(x) x[1])  
      absList  <- absList[order(fstsABS)]; relList  <- relList[order(fstsREL)]
      
      if(length(absList)==0){absList=list(c(0,0,0),c(0,0,0))}
      if(length(relList)==0){relList=list(c(0,0,0),c(0,0,0))}
      
      
      ############################## Evaluation ##############################
      
      
      gmmlHF <- t(matrix(rep("",2*n),n))
      for(i in 1:length(absList)){gmmlHF[1,][absList[[i]]] <- paste0("c",i)} 
      for(i in 1:length(relList)){gmmlHF[2,][relList[[i]]] <- paste0("c",i)}
      
      
      ############################## END REPORT ##############################
    
    temp    <- cbind(cnames,gmmlHF[2,])
    clubs <- temp[!duplicated(temp[,2]),2]
    clubs <- clubs[clubs!=""]
    
    clubs <- lapply(clubs, function(cl) temp[temp[,2]==cl,1])
    clubs <- clubs[order(sapply(clubs,length),decreasing = T)]
    save(clubs, file = paste0("Results/",yearOrRegion,"_HF.rda")) # Saving clubs to disc  
    counts <- table(temp[,2])
    sizes <- sapply(2:max(counts), function(s) sum(counts==s)) # Count of clubs per sizes
    return(sizes)
  }
  
  # Report for all
  ress1 <- lapply(c("Maddison","Penn","Europe","Europe+G7","Europe+S&P","G7+S&P"), applyHF)
  ress2 <- lapply(c(1930,1940), applyHF)
  ress  <- c(ress2,ress1)
  countsAll <- t(sapply(ress, function(cnts) c(cnts,rep(NA,max(sapply(ress,length))-length(cnts)))))

  rownames(countsAll) <- c(1930,1940,"Maddison","Penn","Europe","Europe+G7","Europe+S&P","G7S&P")
  colnames(countsAll) <- paste0('# ',2:(ncol(countsAll)+1))

  write.csv(countsAll,"Results/countsAllHF.csv")
  