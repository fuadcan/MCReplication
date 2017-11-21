clusterAGK<- function(pmat){
#   pmat --> list of clubs
clist <- list()
k<-1
colnames(pmat)<- 1:ncol(pmat)
finitto<-0


  while(finitto==0){
  igr <- graph.adjacency(pmat,"undirected", diag = F)  
  lgr <- largest.cliques(igr) ; lgr<- (sample(lgr,1))[[1]]
  clist[[k]]<- sort(as.numeric(colnames(pmat)[lgr]))
  k<-k+1  
  pmat <- pmat[-lgr,-lgr] 
  if(length(clist[[k-1]])<2){finitto<-1;clist<-clist[1:(k-2)]} else if(prod(dim(pmat)==c(0,0))){finitto<-1}
  }

firsts <- sapply(clist, function(x) x[1])
clist<-clist[order(firsts)]

return(clist)

}