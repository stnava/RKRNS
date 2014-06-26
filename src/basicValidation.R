
print("########basic validation########")
ccafeatspace<-residuals(lm(featspace~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ] ))
# ccafeatspace<-featspace
nl<-nrow(  ccafeatspace )
print(paste("NL rows:",nl))
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
meanfeat1<-apply(ccafeatspace[ inds1  , ], FUN=mean,MARGIN=2)
meanfeat2<-apply(ccafeatspace[ inds2  , ], FUN=mean,MARGIN=2)
docca<-T
if ( docca == TRUE ) {
  longc<-0
  nperm<-0
  nv<-10; its<-100
  mysparse<-c( -0.005, -0.25 )
  myrob<-0
  # c('cross','lake')  ) 
  # c('lake','mountain','stone','beach','river')  ) c('politician')  ) 
  # c("child","woman") c('doctor','terrorist','artist')c('lake','mountain','woman') 
  redlist<-c()
  locwordlist<-c('bird','duck')
  locwordlist<-c('child','woman')   #
  locwordlist<-'.red.'
  locwordlist<-c('tree','bird','green','red') #
  locwordlist<-'coffee'
  locwordlist<-c('lake','mountain','stone','beach','river')
  locwordlist<-c('dime') # ,'green','red') #
  locwordlist<-'criminal'
  locwordlist<-c(  '.coffee.' )
  locwordlist<-c('politician','scientist')
  for ( w in locwordlist ) redlist<-sort(c(redlist,grep(w, fspacenames )))
  wct<-1
  wclasses<-rep(0,length(redlist))
  mywf <- function( x ) return(length(grep(w,x)))
  for ( w in locwordlist ) {
   nz<-wclasses==0
   wclass<-wct*as.numeric(unlist(lapply(fspacenames[redlist],FUN=mywf)))
   wclasses[nz]<-wclasses[nz]+wclass[nz]
   wct<-wct+1
  }
  l1<-1:(length(redlist)*1/2)
  l2<-((max(l1)+1):length(redlist))
  ###########################################################
  sentspace2<-cbind(  log( sentspace - min(sentspace) + 1 ) )
  # sentspace2<-sentspace 
  # multivariate correlation between global bold features and eigensentences
  ccamats1<-list( ( ccafeatspace[ redlist[l1], ] ) , (sentspace2[redlist[l1], ] ) )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, robust=0 ) #, nboot=50 )  # subaal
  pj1<-  ccafeatspace[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig1)
  pj2<- sentspace2[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig2)
  par(mfrow=c(2,1))
  brainregionlist<-list()
  for ( k in 1:nv ) {
#    plot( fcca1$projections[,k], fcca1$projections2[,k] )
#    Sys.sleep(1)
    plot( pj2[,k], pj1[,k] )
    Sys.sleep(1)
    print( cor.test(  pj2[,k], pj1[,k] ) ) 
    mm<-matrix(fcca1$eig1[,k],nrow=responselength)
    myestimatedhrf<-rowMeans(mm)
    plot(myestimatedhrf,type='l')
    mmmag<-sqrt( mm^2 )
    myestimatedbrainregionsval<-apply(mmmag,FUN=sum,MARGIN=2)
    selector<-myestimatedbrainregionsval >  mean(myestimatedbrainregionsval[myestimatedbrainregionsval>0])
    brainregionlist<-lappend(brainregionlist,colnames(imat)[ selector ])
#    for ( i in 1:24 ) { plot(mm[i,],type='l'); Sys.sleep(0.5) }
  }
}
decode<-T
if ( decode ) {
#  eanat1<-sparseDecom( (ccafeatspace[ redlist, ])   , sparseness=-0.9, nvecs=50, its=2, mycoption=1 )
#  decodemat<-as.matrix( eanat1$eig )
  decodemat<-as.matrix(fcca1$eig1)
  decodemat2<-as.matrix(fcca1$eig2)
  mydf <-data.frame( dx=sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l1], ]   %*% decodemat  )
  myudf<-data.frame( dx=sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l2], ]   %*% decodemat )
#  myrf<-RRF( dx ~ . , data=mydf )
  myrf<-svm( dx ~ . ,  mydf )
  pred<-predict( myrf, newdata=myudf )
#  myrf<-bgp(X=mydf[,2:ncol(mydf)],Z=mydf[,1],XX=myudf[,2:ncol(myudf)])
  
  mydata <- data.frame(Real=myudf$dx,Pred=pred,group=fspacenames[redlist[l2]])
  eigSz<-apply(sentspace2[ redlist[l2]  , ],FUN=max,MARGIN=1)*1.5
  chart_title<-locwordlist
  myqplot <- qplot(  x=Real, y=Pred, size=eigSz, colour=factor(group), data=mydata) +
                     scale_size(range=c(2, 5))+labs(title = chart_title)+
                     theme(text = element_text(size=5))
  ggsave("myqplot.pdf")
  # cbind(myudf$dx,pred) )
# print( sum( myudf$dx == pred )/length(pred) )
print( cor.test(myudf$dx,pred) )
}
# print(brainregionlist)
print("Final result of this script: we can predict the low-dimensional representation of the eigenwords as seen by bold.")
print(unique(fspacenames[redlist]))
