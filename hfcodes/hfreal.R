# setwd("C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes")
# setwd("C:/Users/fuadcan/Documents/CADF_TablesN")
# importing datas
# Tm=50;n=20;clsize=5; alphaInt <- c(.01,.1); d=1; egeFact=1; ind=1
write.table(matrix(paste0("c",1:200),1,200), file = "C:/Users/fuadcan/Documents/CADF_TablesN/cnames.txt",row.names = F,col.names = F)

library("stringr")
hfreal<-function(year){
  dirname <-"C:/Users/fuadcan/Dropbox/FuatHarun/MonteCarlo/R codes/"
filename  <- paste0(dirname,"madisonFrom-",year,"-2.csv")
# filename  <- paste0(dirname,"PWT.TXT")
z        <- read.table(filename,header = T,sep = ";")
# z        <- read.table(filename,header = T,sep = ",")
# z<-z[,-1]
n<-ncol(z)

#   shuffling z

# z<-log(z)
z <-z[,sample(n)]
shfnames <- colnames(z)

write.table(dim(z),"C:/Users/fuadcan/Documents/CADF_TablesN/dims.csv",row.names = FALSE,col.names = FALSE)
write.table(z,file = "C:/Users/fuadcan/Documents/CADF_TablesN/deneme1.csv",row.names = FALSE,col.names = FALSE)


tempHF <- shell(paste("C:/gauss10/tgauss -b ",'RUN C:\\Users\\fuadcan\\Documents\\CADF_TablesN\\main2.gss'), intern=TRUE,wait=TRUE)
tempHF <- tempHF[1:(which(tempHF=="GAUSS 10.0.3 (Dec 22 2009, 1345) 32-bit")-2)]
aCrude<- strsplit(tempHF[1:((which(tempHF=="brkpnt"))-1)]," ")
aCrude<-lapply(1:length(aCrude), function(x) aCrude[[x]] <- aCrude[[x]][aCrude[[x]]!=""])
aCrude<-lapply(1:length(aCrude), function(x) as.numeric(str_replace_all(aCrude[[x]],"c","")))
absList<- aCrude[sapply(1:length(aCrude), function(x) length(aCrude[[x]])!=1)]

rCrude  <- strsplit(tempHF[((which(tempHF=="brkpnt"))+1):length(tempHF)]," ")
rCrude  <-lapply(1:length(rCrude), function(x) rCrude[[x]] <- rCrude[[x]][rCrude[[x]]!=""])
rCrude  <-lapply(1:length(rCrude), function(x) as.numeric(str_replace_all(rCrude[[x]],"c","")))
relList <- rCrude[sapply(1:length(rCrude), function(x) length(rCrude[[x]])!=1)]

matrix(c(1:n,shfnames),n,2)

absNames<-list()
for(i in 1:length(absList)){absNames[[i]]<-shfnames[absList[[i]]]}

relNames<-list()
for(i in 1:length(relList)){relNames[[i]]<-shfnames[relList[[i]]]}

return(list(absNames))
}

