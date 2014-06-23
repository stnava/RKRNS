print("########sccan HRF estimation########")
ccafeatspace<-residuals(lm(featspace~ as.numeric( nchar[ eventsw > 0 ] ) ))
nl<-nrow(  featspace )
print(paste("NL rows:",nl))
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
longc<-0
nperm<-0
myrob<-0
myestimatedhrfs<-matrix( rep(0,responselength*length(words)),ncol=responselength  )
ct<-1
for ( w in words ) {
  redlist<-c()
  redlist<-sort(c(redlist,grep(w, fspacenames[inds1] )))
  sentspace2<-cbind(  log( sentspace - min(sentspace) + 1 ) )
  ccamats1<-list( ( ccafeatspace[ redlist, ] ) , (sentspace2[redlist, ] ) )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, robust=myrob) #, nboot=50 )  # subaal
  mm<-matrix(fcca1$eig1[,1],nrow=responselength)
  myestimatedhrf<-apply(mm,FUN=median,MARGIN=1)
  if ( mean(myestimatedhrf) < 0 ) myestimatedhrf=myestimatedhrf*(-1)
  myestimatedhrfs[ct,]<-myestimatedhrf
  ct<-ct+1
}
for ( ct in c(20,55,99,111,155, 180, 222, 240) )
  {
  plot( myestimatedhrfs[ct,] , type='l', main=rownames(myestimatedhrfs)[ct])
  Sys.sleep(0.25)
  }
rownames(myestimatedhrfs)<-words
print("#### now convolve the hrf and do a regression")
k1<-which( words == 'politician' ) 
k2<-which( words == 'red' )
myf1 <- function( x ) conv( x , myestimatedhrfs[k1,] )[1:length(x)]
myf2 <- function( x ) conv( x , myestimatedhrfs[k2,] )[1:length(x)]
dmatsblock1<-apply( dmats[30000:nrow(dmats),] , FUN=myf1, MARGIN=2 )
dmatsblock2<-apply( dmats[30000:nrow(dmats),] , FUN=myf2, MARGIN=2 )

bmat<-as.matrix(imat[30000:nrow(dmats),])
bmatprev<-bmat-ashift( bmat, v=c(12,0) )
bsig<-rowMeans(bmat)
mdl11<-lm( dmatsblock1[,k1]~bmat + bmatprev ) 
mdl12<-lm( dmatsblock1[,k2]~bmat + bmatprev ) 
mdl21<-lm( dmatsblock2[,k1]~bmat + bmatprev ) 
mdl22<-lm( dmatsblock2[,k2]~bmat + bmatprev ) 
print(summary(mdl11)$r.squared)
print(summary(mdl12)$r.squared)
print(summary(mdl21)$r.squared)
print(summary(mdl22)$r.squared)

if ( FALSE ) {
    pdf("myestimatedhrfscor.pdf",width=32,height=32)
pheatmap(cor(t(myestimatedhrfs)))
dev.off()
for ( ct in 1:nrow(myestimatedhrfs) )
  {
  plot( myestimatedhrfs[ct,] , type='l', main=rownames(myestimatedhrfs)[ct])
  Sys.sleep(0.1)
  }
}
