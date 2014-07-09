sccanBasedDecoder <- function( eventdata, designmat, boldFeatureMatrix, sentenceSpace, mysparse=c(-0.1,-0.1), nvecs=5, its=1, smooth=0, cthresh=0, mask4d=NA, strategy=NA, doEanat=F, joinEanat=F, outgraphfile='temp2.html'  )
{
#########################################
# parameters for dimensionality reduct.
#########################################
nv<-nvecs
###########constant params below########
longc<-0
nperm<-0
myrob<-0
#########################################
################################################################################################################################
# train / test split
################################################################################################################################
nl<-nrow(  ccafeatspace )
inds1<-seq(1,(nl-1),by=2)
inds2<-inds1+1
inds1<-c(inds1,inds2[1:(length(inds2)*5/6)])
inds2<-inds2[(length(inds2)*5/6):length(inds2)]
redlist<-c()
redlist<-which( !is.na(eventdata$sentences) )
wclassesf<-as.factor( eventdata$sentlab[redlist] )
whichcols<- colnames(designmat) %in%  eventdata$sentences[redlist]
wclasslevs<-( levels( wclassesf ) )
l1<-seq(1,length(redlist)-1,2)
l2<-l1+1
###########################################################
#  Great!  Now do some cca based dimensionality reduction #
###########################################################
#  sentspace2<-cbind(  log( sentspace - min(sentspace) + 1 ) )
  sentspace2<-sentenceSpace 
  # multivariate correlation between global bold features and eigensentences
  nccavecs<-nv
  perword<-1 # length(locwordlist)
  sccanBdictionary<-matrix( rep(0,ncol(featspace)*perword*nccavecs),nrow=ncol(featspace))
  sccanWdictionary<-matrix( rep(0,ncol(sentspace2)*perword*nccavecs),nrow=ncol(sentspace2))
  wct<-1
  blulist<-redlist[l1]
  classmatrix<-data.matrix( designmat[ eventtimes , whichcols ] )
  classmatrix<-classmatrix[ blulist , ]
  if ( TRUE ) {
    classmatrix<-interleaveMatrixWithItself( classmatrix, ncol(sentenceSpace) )
    for ( i in 1:nrow(classmatrix) )
      {
      esent<-sentenceSpace[blulist[i],]
      classmatrix[i, (classmatrix[i,] > 0 ) ]<-esent
      }
    }
  ccamats1<-list( ( ccafeatspace[ blulist, ] ) , (classmatrix) )
  print(paste("CCA",length(wclasslevs),its))
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=2, perms=nperm, robust=0, smooth=smooth, cthresh = c(cthresh, 0) ,  inmask = c(mask4d, NA), ell1=0.1, z=-1 ) #, nboot=50 )  # subaal mask4d
  if ( typeof(fcca1$eig1[[1]]) != "double" )
    {
    for ( j in 1:nccavecs )
      {
      pmat<-timeseries2matrix( fcca1$eig1[[j]], subaal )
      pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=subaal, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
      sccanBdictionary[,j]<-pmat[1,]
      }
    wct<-wct+nccavecs
    } else sccanBdictionary <- fcca1$eig1
  decodemat<-as.matrix(sccanBdictionary)
  fcca1$eig2<-sccanWdictionary
  for ( k in 1:ncol(decodemat)) 
    {
    kk<-spatioTemporalProjectionImage( decodemat[,k], responselength, sum, subaal )
    myestimatedhrf<-kk$timefunction
    plot(myestimatedhrf,type='l')
    Sys.sleep(0.1)
    antsImageWrite( kk$spaceimage , paste('temp',k,'.nii.gz',sep=''))
    }
  if ( doEanat  )
    {
    if ( is.na( mask4d ) ) 
      eanat1<-sparseDecom(  ccamats1[[1]] , sparseness=mysparse[1], nvecs=nv, its=1, mycoption=0, cthresh = cthresh , smooth=smooth  )
    if (!is.na( mask4d ) ) 
      eanat1<-sparseDecom(  ccamats1[[1]] , sparseness=mysparse[1], nvecs=nv, its=1, mycoption=0, cthresh = cthresh , smooth=smooth , inmask=mask4d )
    sccanBdictionary2<-sccanBdictionary*0
    if ( typeof(eanat1$eig[[1]]) == "double" ) sccanBdictionary2<-as.matrix( eanat1$eig )
    if ( typeof(eanat1$eig[[1]]) != "double" ) 
      {
        for ( j in 1:nv )
          {
          pmat<-timeseries2matrix( eanat1$eig[[j]], subaal )
          pmat<-timeserieswindow2matrix( data.matrix( pmat ), mask=subaal, eventlist=1, timewindow=responselength, zeropadvalue=0 )$eventmatrix
          sccanBdictionary2[,j]<-pmat[1,]
          kk<-spatioTemporalProjectionImage( sccanBdictionary2[,j], responselength, sum, subaal )
          myestimatedhrf<-kk$timefunction
          plot(myestimatedhrf,type='l'); Sys.sleep(1)
          antsImageWrite( kk$spaceimage , paste('temp',j,'.nii.gz',sep=''))
          }
    }
  decodemat<-cbind( decodemat, as.matrix(sccanBdictionary2) )
  }
  if ( joinEanat & nvecs > 3 )
    {
    decodemat2<-decodemat
    kk<-joinEigenanatomy( ccamats1[[1]], mask=NA, decodemat2 , c(1:10)/100 )
    decodemat<-t( kk$fusedlist )
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
  locsents<-which( duplicated( eventdata$sentences ) == FALSE )
  sentencedf<-data.frame( sentences=eventdata$sentences[locsents],
                           ids=eventdata$sentlab[locsents]  )
  sentencesubset<-sentencedf$sentences %in% unique(fspacenames[redlist[l2]])
  nodedf<-data.frame( nodename=sentencedf$sentences[sentencesubset], nodeid=sentencedf$ids[sentencesubset] )
  ww<-  classificationNetwork( nodesIn=nodedf, wclassesf[l2], pred ,outfile=outgraphfile, mycharge=-2066,zoom=T)

  mydata <- data.frame(group=fspacenames[redlist[l2]], Real=myudf$dx,Pred=pred)
  eigSz<-apply(sentspace2[ redlist[l2]  , ],FUN=max,MARGIN=1)*1.5
  chart_title<-locwordlist
  pltsz<-8
  gpic <-  ggplot(mydata,aes(Real,Pred,size=eigSz,color=group,fill=group))+geom_point()+
      guides(colour = guide_legend(override.aes = list(size = pltsz)))+
                         theme(text = element_text(size=pltsz*2)) +
                     scale_size(range=c(pltsz/2, pltsz))
  ggsave("myqplot.pdf",height=8,width=12)
  return( mydata )
}
