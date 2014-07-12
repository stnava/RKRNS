sccanBasedDecoder <- function( eventdata, designmat, boldFeatureMatrix, sentenceSpace, mysparse=c(-0.1,-0.1), nvecs=5, its=1, smooth=0, cthresh=0, mask=NA, strategy=NA, doEanat=F, joinEanat=F, outputfileprefix='sccanBasedDecoder', interleave=FALSE , sentenceTransformation="none", trainset=NA , locwordlist="" )
{
#########################################
# parameters for dimensionality reduct.
#########################################
nv<-nvecs
###########constant params below########
longc<-0
nperm<-0
myrob<-0
if ( ! is.na( mask ) )
  {
  if ( mask@dimension == 4 )
    {
    myarr<-as.array( mask )
    arr3d<-myarr[,,,1]
    mask3d<-as.antsImage( arr3d )
    }
  }
#########################################
################################################################################################################################
# train / test split
################################################################################################################################
redlist<-c()
for ( w in locwordlist ) redlist<-sort(unique(c(redlist,grep(w, fspacenames ))))
wclassesf<-as.factor( eventdata$sentlab[redlist] )
whichcols<- colnames(designmat) %in%  eventdata$sentences[redlist]
wclasslevs<-( levels( wclassesf ) )
if ( identical(trainset,rep(NA,length(trainset))) ) l1<-seq(1,length(redlist)-1,2) 
fcca1<-0
###########################################################
#  Great!  Now do some cca based dimensionality reduction #
###########################################################
sentspace2<-sentenceSpace 
if ( sentenceTransformation == "log" ) 
  sentspace2<-cbind(  log( sentenceSpace - min(sentenceSpace) + 1 ) )
if ( sentenceTransformation == "sim" )
    {
    }
  nccavecs<-nv
  perword<-1 
  sccanBdictionary<-matrix( rep(0,ncol(featspace)*perword*nccavecs),nrow=ncol(featspace))
  sccanWdictionary<-matrix( rep(0,ncol(sentspace2)*perword*nccavecs),nrow=ncol(sentspace2))
  blulist<-redlist[l1]
  classmatrix<-data.matrix( designmat[ eventdata$eventtimes , whichcols ] )
  classmatrixTrain<-classmatrix[  blulist , ]
  if ( interleave ) {
    classmatrixTrain<-interleaveMatrixWithItself( classmatrixTrain, ncol(sentenceSpace) )
    for ( i in 1:nrow(classmatrixTrain) )
      {
      esent<-sentenceSpace[blulist[i],]
      tmp<-(classmatrixTrain[i,] > 0 )
      if (length(tmp)==length(esent))classmatrixTrain[i, tmp ]<-esent
      }
    ccamatsTrain<-list( ( ccafeatspace[ blulist, ] ) , classmatrixTrain )
    } else ccamatsTrain<-list( ( ccafeatspace[ blulist, ] ) , sentenceSpace[blulist,] )
  print(paste("CCA",length(wclasslevs),its))
  fcca1<-sparseDecom2( inmatrix=(ccamatsTrain), nvecs=nv, sparseness=mysparse, its=its, mycoption=1, perms=nperm, robust=0, smooth=smooth, cthresh = c(cthresh, 0) ,  inmask = c(mask, NA), ell1=10 ) #, nboot=50 )
  if ( typeof(fcca1$eig1[[1]]) != "double" )
    {
    for ( j in 1:nccavecs )
      {
      img<-fcca1$eig1[[j]]
      kk<-spatioTemporalProjectionImage( img, sum, mask3d )
      myestimatedhrf<-kk$timefunction
      plot(myestimatedhrf,type='l')
      Sys.sleep( 0.25 )
      antsImageWrite( kk$spaceimage , paste(outputfileprefix,j,'cca.nii.gz',sep=''))
      pmat<-timeseries2matrix( img , mask3d )
      pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=mask3d, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
      sccanBdictionary[,j]<-pmat[1,]
      }
    } else sccanBdictionary <- fcca1$eig1
  decodemat<-as.matrix(sccanBdictionary)
  fcca1$eig2<-sccanWdictionary
  if ( doEanat  )
    {
    if ( is.na( mask ) ) 
      eanat1<-sparseDecom(  ccamatsTrain[[1]] , sparseness=mysparse[1], nvecs=nv, its=1, mycoption=0, cthresh = cthresh , smooth=smooth  )
    if (!is.na( mask ) ) 
      eanat1<-sparseDecom(  ccamatsTrain[[1]] , sparseness=mysparse[1], nvecs=nv, its=1, mycoption=0, cthresh = cthresh , smooth=smooth , inmask=mask )
    sccanBdictionary2<-sccanBdictionary*0
    if ( typeof(eanat1$eig[[1]]) == "double" ) sccanBdictionary2<-as.matrix( eanat1$eig )
    if ( typeof(eanat1$eig[[1]]) != "double" ) 
      {
        for ( j in 1:nv )
          {
          img<-eanat1$eig[[j]]
          kk<-spatioTemporalProjectionImage( img, sum, mask3d )
          myestimatedhrf<-kk$timefunction
          plot(myestimatedhrf,type='l')
          Sys.sleep( 0.25 )
          antsImageWrite( kk$spaceimage , paste(outputfileprefix,j,'cca.nii.gz',sep=''))
          pmat<-timeseries2matrix( img , mask3d )
          pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=mask3d, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
          sccanBdictionary2[,j]<-pmat[1,]
          antsImageWrite( kk$spaceimage , paste(outputfileprefix,j,'eanat.nii.gz',sep=''))
          }
    }
  decodemat<-cbind( decodemat, as.matrix(sccanBdictionary2) )
  }
  if ( joinEanat & nvecs > 3 )
    {
    decodemat2<-decodemat
    kk<-joinEigenanatomy( ccamatsTrain[[1]], mask=NA, decodemat2 , c(1:10)/100 )
    decodemat<-t( kk$fusedlist )
    }
  mydf <-data.frame( dx=wclassesf[l1],  # sentspace2[ redlist[l1]  , ] %*% decodemat2[,1],
                    fsp= scale(ccafeatspace[ redlist[l1], ]   %*% decodemat  ) )
  myudf<-data.frame( dx=wclassesf[-l1], # sentspace2[ redlist[l2]  , ] %*% decodemat2[,1],
                    fsp= scale( ccafeatspace[ redlist[-l1], ]   %*% decodemat ) )
  myrf<-svm(dx ~ . ,mydf , kernel='linear',type='C-classification',probability=TRUE)
#  myrf<-RRF( dx ~ . ,  mydf  ,  ntree=5000 )
  pred<-predict( myrf, newdata=myudf )
  ccaerr<-sum(wclassesf[-l1]==pred)/length(pred)
  randerr<-1.0/length(wclasslevs)
  ccaresult<-paste("CCA-PredErr:",ccaerr*100,"%, vs random",randerr*100,"%")
  print(ccaresult)
  mydata <- data.frame(group=eventdata$sentences[redlist[-l1]], Real=myudf$dx,Pred=pred)
  locsents<-which( duplicated( eventdata$sentences ) == FALSE )
  sentencedf<-data.frame( sentences=eventdata$sentences[locsents],
                           ids=eventdata$sentlab[locsents]  )
  sentencesubset<-sentencedf$sentences %in% unique(eventdata$sentences[redlist[-l1]])
  nodedf<-data.frame( nodename=sentencedf$sentences[sentencesubset], nodeid=sentencedf$ids[sentencesubset] )
  ww<-  classificationNetwork( nodesIn=nodedf, mydata$Real, mydata$Pred ,outfile=paste(outputfileprefix,".html",sep=''), mycharge=-2066,zoom=T)
#  mydata<-cbind( mydata, nodedf )
  return( list(ccaresult=ccaresult,ccapredictions=mydata, nodedf=nodedf, ccaDictionary=decodemat,  ccaobject=fcca1 ) )

  if ( FALSE ) {
#    eigSz<-apply(sentenceSpace[ redlist[-l1]  , ],FUN=max,MARGIN=1)*1.5
    chart_title<-"SCCAN-Decode"
    pltsz<-8
    gpic <-  ggplot(mydata,aes(Real,Pred,color=group,fill=group))+geom_point()+
      guides(colour = guide_legend(override.aes = list(size = pltsz)))+
                         theme(text = element_text(size=pltsz*2)) +
                     scale_size(range=c(pltsz/2, pltsz))
    ggsave(paste(outputfileprefix,".pdf",sep=''),height=8,width=12)
    }
}
