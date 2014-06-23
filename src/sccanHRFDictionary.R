print("########sccan HRF estimation########")
nl<-nrow(  featspace )
print(paste("NL rows:",nl))
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
longc<-0
nperm<-0
nv<-4; its<-20
mysparse<-c( 0, 0 )
mysparse<-c( -0.9, -0.9 )
myrob<-0
myestimatedhrfs<-matrix( rep(0,responselength*length(words)),ncol=responselength  )
ct<-1
for ( w in words ) {
  redlist<-c()
  redlist<-sort(c(redlist,grep(w, fspacenames[inds1] )))
  sentspace2<-cbind(  log( sentspace - min(sentspace) + 1 ) )
  ccamats1<-list( ( featspace[ redlist, ] ) , (sentspace2[redlist, ] ) )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, robust=myrob) #, nboot=50 )  # subaal
  mm<-matrix(fcca1$eig1[,1],nrow=responselength)
  myestimatedhrf<-apply(mm,FUN=median,MARGIN=1)
  myestimatedhrfs[ct,]<-myestimatedhrf
  ct<-ct+1
}
rownames(myestimatedhrfs)<-words
for ( ct in 1:nrow(myestimatedhrfs) )
  {
  plot( myestimatedhrfs[ct,] , type='l', main=rownames(myestimatedhrfs)[ct])
  Sys.sleep(1)
  }
