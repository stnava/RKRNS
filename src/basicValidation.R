#########################################
print("########basic validation########")
#########################################
ccafeatspace<-residuals(lm(featspace~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ] ))
# ccafeatspace<-featspace
nl<-nrow(  ccafeatspace )
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
meanfeat1<-apply(ccafeatspace[ inds1  , ], FUN=mean,MARGIN=2)
meanfeat2<-apply(ccafeatspace[ inds2  , ], FUN=mean,MARGIN=2)
docca<-T
if ( docca == TRUE ) {
  longc<-0
  nperm<-0
  nv<-5; its<-100
  mysparse<-c( -0.25, -0.25 )
  myrob<-0
  # c('cross','lake')  ) 
  # c('lake','mountain','stone','beach','river')  ) c('politician')  ) 
  # c("child","woman") c('doctor','terrorist','artist')c('lake','mountain','woman') 
  redlist<-c()
  locwordlist<-c('child','woman')   #
  locwordlist<-'.red.'
  locwordlist<-c('tree','bird','green','red') #
  locwordlist<-'coffee'
  locwordlist<-c('lake','mountain','stone','beach','river')
  locwordlist<-c('dime') # ,'green','red') #
  locwordlist<-'criminal'
  locwordlist<-c('hotel')
  locwordlist<-c('politician','scientist')
  locwordlist<-c('lake','mountain','stone','beach','river')
  locwordlist<-c('bird','duck')
  locwordlist<-c(  '.coffee.' )
  for ( w in locwordlist ) redlist<-sort(unique(c(redlist,grep(w, fspacenames ))))
  wct<-1
  wclasses<-rep(0,length(redlist))
  wclasslevs<-(unique(fspacenames[redlist]))
  print(paste("NL rows:",nl,"classes:",length(wclasslevs)))
  mywf <- function( x ) return(length(grep(w,x)))
  for ( w in fspacenames[redlist] ) {
      wclasses[wct]<-which( w == wclasslevs )
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
wclasses<-as.factor(wclasses)
decode<-T
if ( decode ) {
#  eanat1<-sparseDecom( (ccafeatspace[ redlist, ])   , sparseness=-0.9, nvecs=50, its=2, mycoption=1 )
#  decodemat<-as.matrix( eanat1$eig )
  decodemat<-as.matrix(fcca1$eig1)
  decodemat2<-as.matrix(fcca1$eig2)
  mydfc <-data.frame( dx=sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l1], ]   %*% decodemat  )
  myudfc<-data.frame( dx=sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l2], ]   %*% decodemat )
  myrf<-svm( dx ~ . , data=mydfc ,  ntree=5000 ) # , ntry=2000 , mtry=5 )
  pred<-predict( myrf, newdata=myudfc )
  mydata <- data.frame(Real=myudf$dx,Pred=pred,group=fspacenames[redlist[l2]])
  eigSz<-apply(sentspace2[ redlist[l2]  , ],FUN=max,MARGIN=1)*1.5
  chart_title<-locwordlist
  myqplot <- qplot(  x=Real, y=Pred, size=eigSz, colour=factor(group), data=mydata) +
                     scale_size(range=c(2, 5))+labs(title = chart_title)+
                     theme(text = element_text(size=5))
  ggsave("myqplot.pdf")

  
  mydf <-data.frame( dx=wclasses[l1],  # sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l1], ]   %*% decodemat  )
  myudf<-data.frame( dx=wclasses[l2], # sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l2], ]   %*% decodemat )
  myrf<-svm( dx ~ . ,  mydf )
#  myrf<-bgp(sX=mydf[,2:ncol(mydf)],Z=mydf[,1],XX=myudf[,2:ncol(myudf)])
  pred<-predict( myrf, newdata=myudf )
  print(paste("PredErr:",sum(wclasses[l2]==pred)/length(pred),1.0/length(wclasslevs)))
  
  # cbind(myudf$dx,pred) )
# print( sum( myudf$dx == pred )/length(pred) )
if ( typeof(pred) != "integer" ) print( cor.test(myudf$dx,pred) )
}
# print(brainregionlist)
print("Final result of this script: we can predict the low-dimensional representation of the eigenwords as seen by bold.")
print(unique(fspacenames[redlist]))


decode2<-FALSE # just a comparison to CCA 
if ( decode2 ) {
  svdred<-svd( featspace[ redlist[l1], ] , nv=200, nu=0 )
  mydf <-data.frame( dx=wclasses[l1],  fsp=featspace[ redlist[l1], ]  %*% svdred$v )
  myudf <-data.frame( dx=wclasses[l2], fsp=featspace[ redlist[l2], ]  %*% svdred$v )
  myrf<-svm( dx ~ . ,  mydf )
  pred2<-predict( myrf, newdata=myudf )
  print(paste("PredErr:",sum(wclasses[l2]==pred2)/length(pred2),1.0/length(wclasslevs)))
}
