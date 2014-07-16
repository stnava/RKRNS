templateBasedClassification <- function( exemplarmat, labels, newmat,
                                        method="corr", mask=NA  )
{
# create the dictionary
classlabels<-sort(unique(labels))
nclasses<-length(classlabels)
featuretemplate<-matrix( rep(0,ncol(exemplarmat)*nclasses), ncol=ncol(exemplarmat) )
for ( i in 1:nclasses )
  {
  whichclass<-classlabels[i]
  featuretemplate[i,]<-apply(  exemplarmat[labels==whichclass,] , FUN=mean, MARGIN=2)
  }
rownames(featuretemplate)<-classlabels
#svdfeaturetemplate<-t(svd( exemplarmat )$v[,1:10])
#rownames(svdfeaturetemplate)<-paste("SVD",1:nrow(svdfeaturetemplate))
#print(cor(t(svdfeaturetemplate),t(featuretemplate)))
if ( method == "dict" )
  {
  testmat<-featuretemplate
  votes<-cor(t(testmat),(newmat[1,]))
  return(list(votes=votes, featuretemplate=featuretemplate ) )
  }
if ( method == "eanat" )
  {
  initlist<-list()
  if ( is.na(mask) ) {
    maskmat<-newmat*0
    maskmat[1,]<-1
    mask<-as.antsImage( maskmat )
  }
  eanatnames<-rep(as.character("A"),nclasses*2)
  ct<-1
  for ( i in 1:nclasses ) {
    vecimg<-antsImageClone( mask )
    vecimg[ mask == 1 ]<-featuretemplate[i,]
    initlist<-lappend(initlist,vecimg)
    initlist<-lappend(initlist,vecimg)
    eanatnames[ct]<-toString(classlabels[i])
    eanatnames[ct+1]<-eanatnames[ct]
    ct<-ct+2
  }
  eanat<-sparseDecom(exemplarmat,inmask=mask, nvecs=length(initlist),
                     sparseness=0.01,  mycoption=0, z=0.5, smooth=0.1,
                     initializationList=initlist, cthresh=250, its=1 )
  eanatmat<-imageListToMatrix( eanat$eigenanatomyimages, mask )
  rownames(eanatmat)<-eanatnames
  mycor<-cor( t(newmat) , t(eanatmat))
  print(mycor[1,])
  pheatmap( rbind(mycor[1,],mycor[1,]))
  fclass<-which.max(abs(mycor[1,]))
  mycor[1,fclass]<-0
  sclass<-which.max(abs(mycor[1,]))
  return(paste(eanatnames[c(fclass,sclass)]))
  }
testmat<-exemplarmat
votes<-rep(0,nclasses)
for ( j in 1:nrow(testmat) )
  {
  locvotes<-rep(0,nclasses)
  for ( i in 1:nclasses )
    {
#    locvotes[i]<-mean(abs(testmat[i,]-newmat[1,]))
    locvotes[i]<-(cor(testmat[i,],newmat[1,]))
#    locvotes[i]<-cosineDist(testmat[i,],newmat[1,])
    }
  votes<-votes+locvotes/sum(abs(locvotes))
  }
votes<-votes/sum(abs(votes))
names(votes)<-classlabels
print( summary( lm( newmat[1,] ~ t(featuretemplate) ) ) )
return(list(votes=votes, featuretemplate=featuretemplate ) )
# rownames(featuretemplate)<-sentsoi
# pheatmap(cor( t(ccafeatspace[grep(wordoi,eventdata$sentences),]), t(featuretemplate )))
# j<-44; i<-grep(wordoi,eventdata$sentences)[j]; print( summary( lm( ccafeatspace[i,] ~ t(featuretemplate) ) ) ); eventdata$sentences[i]
}
