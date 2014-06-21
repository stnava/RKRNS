print("########basic validation########")
nl<-nrow(  featspace )
print(paste("NL rows:",nl))
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
meanfeat1<-apply(featspace[ inds1  , ], FUN=mean,MARGIN=2)
meanfeat2<-apply(featspace[ inds2  , ], FUN=mean,MARGIN=2)
docca<-T
if ( docca == TRUE ) {
  longc<-0
  nperm<-0
  nv<-4; its<-20
  mysparse<-c( 0, 0 )
  mysparse<-c( -0.9, -0.5 )
  myrob<-0
  redlist<-c()
  for ( w in c('terrorist','artist') ) redlist<-sort(c(redlist,grep(w, fspacenames )))
#  for ( w in sentences[c(25:55)] ) redlist<-sort(c(redlist,grep(w, fspacenames )))
  # c("child","woman") c('doctor','terrorist','artist')
  #  for ( w in c('child') ) redlist<-sort(c(redlist,grep(w, fspacenames )))
  redlist2<-c()
  for ( w in c('doctor','terrorist') ) redlist2<-sort(c(redlist2,grep(w, fspacenames )))
  l1<-1:(length(redlist)/2)
  l2<-((max(l1)+1):length(redlist))
  sentspace2<-cbind(  sentspace  )
  # multivariate correlation between global bold features and eigensentences
  ccamats1<-list( ( featspace[ redlist[l1], ] ) , (sentspace2[redlist[l1], ] ) )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, robust=myrob) #, nboot=50 )  # subaal
  pj1<-  featspace[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig1)
  pj2<- sentspace2[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig2)
  for ( k in 1:1 ) {
    plot( fcca1$projections[,k], fcca1$projections2[,k] )
    Sys.sleep(1)
    plot( pj2[,k], pj1[,k] )
    Sys.sleep(1)
    print( cor.test(  pj2[,k], pj1[,k] ) ) 
    mm<-matrix(fcca1$eig1[,k],nrow=responselength)
    myestimatedhrf<-rowMeans(mm)
    plot(myestimatedhrf,type='l')
    myestimatedbrainregions<-colMeans(mm)
    Sys.sleep(1)
#    for ( i in 1:24 ) { plot(mm[i,],type='l'); Sys.sleep(0.5) }
  }
}
