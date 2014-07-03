#########################################
# parameters
#########################################
dosvd<-F
docca<-T
nv<-30; its<-41 # cca params
nvsvm<-8      # svd params
mysparse<-c(  -1/(nv-1),  -1/(nv-1) )
###########constant params below########
longc<-0
nperm<-0
myrob<-0
#########################################
#########################################
#########################################
print("########basic validation########")
#########################################
if ( ! exists("eventdata") ) {
print("########basic speech parts########")
sentenceids<-1:length(unique(sentences))
sentencedf<-data.frame(  sentences=sentences, sentenceids=sentenceids )
wordids<-1:length(unique(words))
worddf<-data.frame( words=words, wordids=wordids )
eventdata<-data.frame(eventtimes=eventtimes,sentences=fspacenames)
enouns1<-rep(NA,nrow(eventdata))
enouns1lab<-rep(NA,nrow(eventdata))
everbs1<-rep(NA,nrow(eventdata))
everbs1lab<-rep(NA,nrow(eventdata))
enouns2<-rep(NA,nrow(eventdata))
enouns2lab<-rep(NA,nrow(eventdata))
sentlab<-rep(NA,nrow(eventdata))
sent_token_annotator <- Maxent_Sent_Token_Annotator()
word_token_annotator <- Maxent_Word_Token_Annotator()
pos_tag_annotator <- Maxent_POS_Tag_Annotator()
wordcols<-colnames(dmatw)
for ( i in 1:nrow(eventdata) ) {
    s<-eventdata$sentences[i]
    s<-gsub('[.]',' ',s)
    realwords<-unlist( strsplit(s,' ') )
    rpl<-as.character(sentencedf$sentences) == eventdata$sentences[i]
    if ( sum(rpl) > 0 ) sentlab[i]<-sentencedf$sentenceids[ rpl ]
    locwords<-wordcols[ dmatw[ eventdata$eventtimes[i], ] == 1 ]
    # find position of locwords in realwords
    locwordpos<-rep(0,length(locwords))
    ct<-1
    for ( k in locwords ) {
      pos<-1
      for ( j in realwords ) {
        if( length( grep(k,j) ) > 0 ) locwordpos[ct]<-pos;
        pos<-pos+1
        }
      ct<-ct+1
    }
    a2 <- annotate(s, list(sent_token_annotator, word_token_annotator))
    a3 <- annotate(s, pos_tag_annotator, a2)
    a3w <- subset(a3, type == "word")
    tags <- sapply(a3w$features, `[[`, "POS")
    grepverb<-grep("V",tags[locwordpos])
    if (length(grepverb)>0){
    everbs1[i]<-locwords[grepverb]
    rpl<-as.character(worddf$words) == everbs1[i]
    if ( sum(rpl) > 0 ) everbs1lab[i]<-worddf$wordids[ rpl ]
    }
    grepnouns<-grep("NN",tags[locwordpos])
    enouns1[i]<-locwords[grepnouns[1]]
    rpl<-as.character(worddf$words) == enouns1[i]
    if ( sum(rpl) > 0 ) enouns1lab[i]<-worddf$wordids[ rpl ]
    nextnoun<-grepnouns[1]
    if ( length(grepnouns) > 1 ) nextnoun<-grepnouns[2]
    enouns2[i]<-locwords[nextnoun]
    rpl<-as.character(worddf$words) == enouns2[i]
    if ( sum(rpl) > 0 ) enouns2lab[i]<-worddf$wordids[ rpl ]
}
eventdata<-cbind( eventdata, enouns1=enouns1, enouns1lab=enouns1lab, everbs1=everbs1, everbs1lab=everbs1lab, enouns2=enouns2, enouns2lab=enouns2lab, sentlab=sentlab )
} # check if data exists
################################################################################################################################
################################################################################################################################
if ( !exists("ccafeatspace") ) ccafeatspace<-residuals(lm(featspace~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ]  ))
nl<-nrow(  ccafeatspace )
inds1<-seq(1,(nl-1),by=2)
inds2<-inds1+1
inds1<-c(inds1,inds2[1:(length(inds2)*5/6)])
inds2<-inds2[(length(inds2)*5/6):length(inds2)]
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
locwordlist<-c('politician','scientist','judge','doctor','artist')
locwordlist<-c('yellow','white','blue','black','green','red') #
locwordlist<-c('lake','mountain','stone','beach','river','tree')
locwordlist<-c('judge','criminal')
locwordlist<-'coffee'
####################################  
wordcounts<-rep(0,length(words))
wct<-1 ; l1<-length(fspacenames)/2
for ( w in words ) {
    myct<-length( unique( grep(w,unique(fspacenames[1:l1])) ) )
    wordcounts[wct]<-myct
    wct<-wct+1
}
#  locwordlist<-words[ wordcounts > 5 ]
designmat<-dmats # dmats/w
redlist<-which( !is.na(eventdata$sentences) )
wclassesf<-as.factor( eventdata$sentlab[redlist] )
whichcols<- colnames(designmat) %in%  eventdata$sentences[redlist]
wclasslevs<-( levels( wclassesf ) )
l1<-seq(1,length(redlist)-1,2)
l2<-l1+1
if ( dosvd )
    {
    print(paste("SVD",length(wclasslevs)))
    svdred<-svd( ccafeatspace[ redlist[l1], ] , nv=nvsvm, nu=0 )
    mydf <-data.frame( dx=wclassesf[l1],  fsp=ccafeatspace[ redlist[l1], ]  %*% svdred$v )
    myudf <-data.frame( dx=wclassesf[l2], fsp=ccafeatspace[ redlist[l2], ]  %*% svdred$v )
    myrf<-svm(dx ~ . ,mydf, kernel='linear',type='C-classification',probability=TRUE)
    pred2<-predict( myrf , newdata=myudf )
    svmerr<-sum(wclassesf[l2]==pred2)/length(pred2)
    randerr<-1.0/length(wclasslevs)
    svmresult<-paste("SVM-PredErr:",svmerr*100,"%, vs random",randerr*100,"%")
    print(svmresult)
    }
###########################################################
#  Great!  Now do some cca based dimensionality reduction #
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
    classmatrix<-data.matrix( designmat[ eventtimes , whichcols ] )
    classmatrix<-classmatrix[ blulist , ]
    if ( T ) {
    classmatrix<-interleaveMatrixWithItself( classmatrix, eigsentbasislength )
    for ( i in 1:nrow(classmatrix) )
      {
      esent<-sentspace[blulist[i],]
      classmatrix[i, (classmatrix[i,] > 0 ) ]<-esent
      }
    }
    ccamats1<-list( ( ccafeatspace[ blulist, ] ) , (classmatrix) )
#    mysparse<-c(mysparse[1], -1.0/length(wclasslevs))
#    ccamats1<-list(       ( ccafeatspace[ blulist, ] ) ,       ( sentspace2[ blulist, ] ) )
    antsSetSpacing(mask4d, c(rep(0.5,3),0.5) )
if ( docca == T ) {
    print(paste("CCA",length(wclasslevs)))
    if ( ! exists("fcca1") )
        fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, ell1=10 , perms=nperm, robust=0, smooth=0., cthresh = c(1250, 0) ,  inmask = c(mask4d, NA)) #, nboot=50 )  # subaal
    if ( typeof(fcca1$eig1[[1]]) != "double" )  {
      vislist<-list()
      for ( j in 1:nccavecs ) {
        viewimg<-projectImageAlongAxis(  fcca1$eig1[[j]], subaal )
        viewimg[ subaal > 0 ]<-viewimg[ subaal > 0 ]/max( viewimg[ subaal > 0 ] )
        vislist<-lappend( vislist, viewimg )
        pmat<-timeseries2matrix( fcca1$eig1[[j]], subaal )
        pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=subaal, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
        sccanBdictionary[,j]<-pmat[1,]
#        sccanWdictionary[,wct+(j-1)]<-fcca1$eig2[,j]
      }
    plotANTsImage( ref, vislist , slices='12x56x2' , thresh="0.25x1", color=rainbow(length(vislist)), outname="eanatviz.jpg" )
    fcca1$eig1<-sccanBdictionary
    wct<-wct+nccavecs
  }
#  }
  decodemat<-as.matrix( fcca1$eig1 )
  if ( TRUE ) {
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
}


if ( TRUE  ) {
  if ( docca == F ) {
    eanat1<-sparseDecom(  ( ccafeatspace[ blulist, ] )   , sparseness=mysparse[1], nvecs=nv, its=1, mycoption=1, cthresh = 250 , inmask=mask4d )
    if ( typeof(eanat1$eig[[1]]) == "double" ) decodemat<-as.matrix( eanat1$eig ) else {
        for ( j in 1:nv ) {
          pmat<-timeseries2matrix( eanat1$eig[[j]], subaal )
          pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=subaal, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
          sccanBdictionary[,j]<-pmat[1,]
          decodemat<-as.matrix( sccanBdictionary )
      }}
  }
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
  print(svmresult)

  sentencesubset<- sentencedf$sentences %in% unique(fspacenames[redlist[l2]])
  nodedf<-data.frame( nodename=sentencedf[sentencesubset,1], nodeid=sentencedf[sentencesubset,2] )
  ww<-  misclassnetwork( nodesIn=nodedf, wclassesf[l2], pred ,outfile='temp2.html', mycharge=-2066,zoom=T)

  mydata <- data.frame(group=fspacenames[redlist[l2]], Real=myudf$dx,Pred=pred)
  eigSz<-apply(sentspace2[ redlist[l2]  , ],FUN=max,MARGIN=1)*1.5
  chart_title<-locwordlist
  pltsz<-8
  gpic <-  ggplot(mydata,aes(Real,Pred,size=eigSz,color=group,fill=group))+geom_point()+
      guides(colour = guide_legend(override.aes = list(size = pltsz)))+
                         theme(text = element_text(size=pltsz*2)) +
                     scale_size(range=c(pltsz/2, pltsz))
  ggsave("myqplot.pdf",height=8,width=12)

}
#
# fglobsig<-apply(featspace,FUN=mean,MARGIN=1)
# ccafeatspace<-feat2$amplitudeTransform
# ccafeatspace<-residuals(lm(feat2$frequencyTransform~ 1+as.numeric( nchar[ eventsw > 0 ] ) + eventss[ eventsw > 0 ]  ))
# ccafeatspace<-featspace
#
