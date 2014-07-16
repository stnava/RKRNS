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
  nreps<-2
  eanatnames<-rep(as.character("A"),nclasses*nreps)
  ct<-1
  for ( i in 1:nclasses ) {
    vecimg<-antsImageClone( mask )
    vecimg[ mask == 1 ]<-featuretemplate[i,]
    for (  nr in 1:nreps )
      {
      initlist<-lappend(initlist,vecimg)
      eanatnames[ct+nr-1]<-toString(classlabels[i])
      }
    ct<-ct+nreps
  }
  eanat<-sparseDecom(exemplarmat,inmask=mask, nvecs=length(initlist),
                     sparseness=0.01,  mycoption=1,  smooth=0.0, z=1/nclasses,
                     initializationList=initlist, cthresh=0, its=5 )
  print(eanatnames)
  eanatmat<-imageListToMatrix( eanat$eigenanatomyimages, mask )
  print(dim(eanatmat))
  rownames(eanatmat)<-eanatnames
  mycor<-rep(0,length(initlist))
  for ( i in 1:length(initlist) )
    {
    x<-newmat[1,]
    y<-eanatmat[i,]
#    ww<-which( x != 0 & y!= 0 )
#    mycor[i]<-cor.test( x, y , method="spearman" )$est
#    mycor[i]<-sqrt( sum( ( x[ww] - y[ww] )^2 ) )*(-1)
        mycor[i]<- tempclassrobcor( x, y )
#    mycor[i]<-robcosineSim( x, y )
    }
  names(mycor)<-eanatnames
  print(mycor)
#  pheatmap( rbind(mycor,mycor))
  fclass<-which.max(abs(mycor))
  mycor[fclass]<-0
  sclass<-which.max(abs(mycor))
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


robcosineSim<-function (xin, yin) 
{
    ww<-which( xin != 0 & yin!= 0 )
    x <- t(as.matrix(xin[ww]))
    y <- t(as.matrix(yin[ww]))
    return(as.numeric(1 - x %*% t(y)/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))))
}


tempclassrobcor <-function(  x, y )
{
ww<-which( x != 0 & y!= 0 )
# plot( x[ww] , y[ww] )
return( cor( x[ww] , y[ww] ) )
}
