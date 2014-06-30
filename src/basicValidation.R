#########################################
print("########basic validation########")
#########################################
#fglobsig<-apply(featspace,FUN=mean,MARGIN=1)
if ( !exists("ccafeatspace") ) ccafeatspace<-residuals(lm(featspace~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ]  ))
# ccafeatspace<-feat2$amplitudeTransform
# ccafeatspace<-residuals(lm(feat2$frequencyTransform~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ]  ))
# ccafeatspace<-featspace
nl<-nrow(  ccafeatspace )
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
meanfeat1<-apply(ccafeatspace[ inds1  , ], FUN=mean,MARGIN=2)
meanfeat2<-apply(ccafeatspace[ inds2  , ], FUN=mean,MARGIN=2)
docca<-T
sentenceids<-1:length(unique(sentences))
sentencedf<-data.frame( sentences=sentences, sentenceids=sentenceids )
wordids<-1:length(unique(words))
worddf<-data.frame( words=words, wordids=wordids )
if ( docca == TRUE ) {
  longc<-0
  nperm<-0
  nv<-50; its<-25
  nvsvm<-8
  mysparse<-c(  -1/(nv-1),  1/(nv-1) )
  myrob<-0
  # c('lake','mountain','stone','beach','river')  ) c('politician')  ) 
  # c("child","woman") c('doctor','terrorist','artist')c('lake','mountain','woman') 
  redlist<-c()
  locwordlist<-'.red.'
  locwordlist<-c('tree','bird','green','red') #
  locwordlist<-c('politician','scientist')
  locwordlist<-c('lake','mountain','stone','beach','river')
  locwordlist<-c('criminal','terrorist','doctor','artist')
  locwordlist<-c('child')   #
  locwordlist<-c('green','red') #
  locwordlist<-c('cross','lake') 
  locwordlist<-c("young","yellow")
  locwordlist<-c('dime' ,'green','red') #
  locwordlist<-c('bird','duck')
  locwordlist<-c('child','woman','man','boy','girl')   #
  locwordlist<-c('summer','winter')
  locwordlist<-c('eat','drink')
  locwordlist<-c('fish','bird','dog')
  locwordlist<-'coffee'
  locwordlist<-c('politician','scientist','judge','doctor','artist')
  locwordlist<-c('yellow','white','blue','black','green','red') #
  locwordlist<-c('lake','mountain','stone','beach','river','tree')
  locwordlist<-c('judge','criminal')
####################################  
  wordcounts<-rep(0,length(words))
  wct<-1 ; l1<-length(fspacenames)/2
  for ( w in words ) {
      myct<-length( unique( grep(w,unique(fspacenames[1:l1])) ) )
      wordcounts[wct]<-myct
      wct<-wct+1
  }
#  locwordlist<-words[ wordcounts > 5 ]
  print(locwordlist)
  for ( w in locwordlist ) redlist<-sort(unique(c(redlist,grep(w, fspacenames ))))
  l1<-length(redlist)/2
  l1<-1:l1
  l2<-((max(l1)+1):length(redlist))
  wct<-1
  wclasses<-rep(0,length(redlist))
  wclasslevs<-(unique(fspacenames[redlist]))
  print(paste("NL rows:",nl,"classes:",length(wclasslevs)))
  mywf <- function( x ) return(length(grep(w,x)))
  for ( w in fspacenames[redlist] ) {
      wclasses[wct]<-sentencedf$sentenceids[  sentencedf$sentences == w ]
      wct<-wct+1
  }
  wclassesf<-as.factor(wclasses)
  decode2<-T # just a comparison to CCA 
  pltsz<-8
  if ( decode2 )
    {
    svdred<-svd( ccafeatspace[ redlist[l1], ] , nv=nvsvm, nu=0 )
    mydf <-data.frame( dx=wclasses[l1],  fsp=ccafeatspace[ redlist[l1], ]  %*% svdred$v )
    myudf <-data.frame( dx=wclasses[l2], fsp=ccafeatspace[ redlist[l2], ]  %*% svdred$v )
    myrf<-svm(dx ~ . ,mydf, kernel='linear',type='C-classification',probability=TRUE)
                                        #  myrf<-RRF( dx ~ . ,  mydf ,  ntree=5000 )
    pred2<-predict( myrf, newdata=myudf )
    svmerr<-sum(wclassesf[l2]==pred2)/length(pred2)
    randerr<-1.0/length(wclasslevs)
    svmresult<-paste("SVM-PredErr:",svmerr*100,"%, vs random",randerr*100,"%")
    print(svmresult)
    mydata <- data.frame(group=fspacenames[redlist[l2]], Real=wclasses[l2]+rnorm(length(l2))*0.1,Pred=pred2)
    gpic <-  ggplot(mydata,aes(Real,Pred,color=group,fill=group))+geom_point()+
        guides(colour = guide_legend(override.aes = list(size = pltsz)))+
                         theme(text = element_text(size=pltsz*2)) +
                     scale_size(range=c(pltsz/2, pltsz))
    ggsave("myqplot_svm.pdf",height=8,width=12)
    }
###########################################################
  sentspace2<-cbind(  log( sentspace - min(sentspace) + 1 ) ) #  sentspace2<-sentspace 
  # multivariate correlation between global bold features and eigensentences
  nccavecs<-nv
  perword<-1 # length(locwordlist)
  sccanBdictionary<-matrix( rep(0,ncol(featspace)*perword*nccavecs),nrow=ncol(featspace))
  sccanWdictionary<-matrix( rep(0,ncol(sentspace2)*perword*nccavecs),nrow=ncol(sentspace2))
  wct<-1
#  for ( myw in locwordlist ) {
#    blulist<-sort(unique(c(grep(myw, fspacenames[l1] ))))
    blulist<-redlist[l1]
    ccamats1<-list( whiten( ccafeatspace[ blulist, ] ) , whiten( sentspace2[ blulist, ] ) )
#    ccamats1<-list(       ( ccafeatspace[ blulist, ] ) ,       ( sentspace2[ blulist, ] ) )
    antsSetSpacing(mask4d, c(rep(0.5,3),0.5) )
    fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, inmask = c(mask4d, NA), robust=0, smooth=0., cthresh = c(0, 0) ) #, nboot=50 )  # subaal
   for ( j in 1:nccavecs ) {
      pmat<-timeseries2matrix( fcca1$eig1[[j]], subaal )
      pmat<-timeserieswindow2matrix( data.matrix( pmat ), subaal, 4, responselength, 0 )$eventmatrix
      sccanBdictionary[,wct+(j-1)]<-pmat[1,]
      sccanWdictionary[,wct+(j-1)]<-fcca1$eig2[,j]
    }
    wct<-wct+nccavecs
#  }
  fcca1$eig1<-sccanBdictionary
  fcca1$eig2<-sccanWdictionary
  pj1<-  ccafeatspace[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig1)
  pj2<- sentspace2[ redlist[l2]  , ]  %*%  as.matrix(fcca1$eig2)
  par(mfrow=c(2,1))
  brainregionlist<-list()
  for ( k in 1:1 ) {
    plot( pj2[,k], pj1[,k] )
    mm<-matrix(fcca1$eig1[,k],nrow=responselength)
    myestimatedhrf<-rowMeans(mm)
    whtimes<-which( abs(myestimatedhrf) > 0 )
    plot(myestimatedhrf,type='l')
    mmmag<-sqrt( mm^2 )
    myestimatedbrainregionsval<-apply(mmmag[whtimes,],FUN=sum,MARGIN=2)
    myestimatedbrainregionsval<-myestimatedbrainregionsval/sum(myestimatedbrainregionsval)
    myccavec<-antsImageClone( subaal )
    myccavec[ subaal > 0 ]<-myestimatedbrainregionsval
    antsImageWrite(myccavec,'temp.nii.gz')
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
  mydfc <-data.frame( dx=sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l1], ]   %*% decodemat  )
  myudfc<-data.frame( dx=sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= ccafeatspace[ redlist[l2], ]   %*% decodemat )
  myrf<-RRF( dx ~ . , data=mydfc ) # , ntry=2000 , mtry=5 )
  predc<-predict( myrf, newdata=myudfc )
  mydata <- data.frame(group=fspacenames[redlist[l2]], Real=myudfc$dx,Pred=predc)
  eigSz<-apply(sentspace2[ redlist[l2]  , ],FUN=max,MARGIN=1)*1.5
  chart_title<-locwordlist
  pltsz<-8
  gpic <-  ggplot(mydata,aes(Real,Pred,size=eigSz,color=group,fill=group))+geom_point()+
      guides(colour = guide_legend(override.aes = list(size = pltsz)))+
                         theme(text = element_text(size=pltsz*2)) +
                     scale_size(range=c(pltsz/2, pltsz))
  ggsave("myqplot.pdf",height=8,width=12)

  mydf <-data.frame( dx=wclassesf[l1],  # sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= (ccafeatspace[ redlist[l1], ]   %*% decodemat  ) )
  myudf<-data.frame( dx=wclassesf[l2], # sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= ( ccafeatspace[ redlist[l2], ]   %*% decodemat ) )
  myrf<-svm(dx ~ . ,mydf , kernel='linear',type='C-classification',probability=TRUE)
#  myrf<-randomForest( dx ~ . ,  mydf  ,  ntree=5000 )
  pred<-predict( myrf, newdata=myudf )
  ccaerr<-sum(wclassesf[l2]==pred)/length(pred)
  randerr<-1.0/length(wclasslevs)
  ccaresult<-paste("CCA-PredErr:",ccaerr*100,"%, vs random",randerr*100,"%")
  print(ccaresult)
  source(paste(srcdir,"setup.R",sep=''))
  sentencesubset<- sentencedf$sentences %in% unique(fspacenames[redlist[l2]])
  nodedf<-data.frame( nodename=sentencedf[sentencesubset,1], nodeid=sentencedf[sentencesubset,2] )
  ww<-  misclassnetwork( nodesIn=nodedf, wclassesf[l2], pred ,outfile='temp2.html', mycharge=-2066,zoom=T) 
  ww<-  misclassnetwork( nodesIn=nodedf, wclassesf[l2], pred ,outfile='temp.html',zoom=F, mycharge=-2066, whichviz="Simple")
  
  pltsz<-10
  mydata <- data.frame(group=fspacenames[redlist[l2]], Real=wclasses[l2]+rnorm(length(l2))*0.1,Pred=pred)
  gpic <-  ggplot(mydata,aes(Real,Pred,color=group,fill=group))+geom_point()+
      guides(colour = guide_legend(override.aes = list(size = pltsz)))+
          theme(text = element_text(size=pltsz)) +
              scale_size(range=c(pltsz/2, pltsz))
  ggsave("myqplot_cca.pdf",height=8,width=12)

  
  mydata <- data.frame(group=fspacenames[redlist[l2]], Real=myudfc$dx,Predc=predc, Pred=pred)
  mydata$group <- as.character(mydata$group)
  m1 <- rPlot(Predc ~ Real  | group  ,type = 'point', color='group', data = mydata)
  m2 <- rPlot(Predc ~ Real           ,type = 'point', color='group', data = mydata)
  m2$addParams(width = 1200, height = 600, dom = 'chart1',
    title = "SCCAN-Decode")
#  m2$publish('r1', id='17af51c7e4b8c809b096' ) # host='gist' )
  

  
  # myrf<-bgp(sX=mydf[,2:ncol(mydf)],Z=mydf[,1],XX=myudf[,2:ncol(myudf)])
  # cbind(myudf$dx,pred) )
# print( sum( myudf$dx == pred )/length(pred) )
if ( typeof(pred) != "integer" ) print( cor.test(myudf$dx,pred) )
}
# print(brainregionlist)
#print("Final result of this script: we can predict the low-dimensional representation of the eigenwords as seen by bold.")
# print(unique(fspacenames[redlist]))


