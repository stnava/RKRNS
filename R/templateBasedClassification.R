templateBasedClassification <- function( exemplarmat, labels, newmat,
                                        method="corr", mask=NA, eigsents=NA , sparval=c(-0.2,-0.9) )
{
mylocaldistfun <- robcosineSim
mylocaldistfun <- tempclassrobcor
mylocaldistfun <- overlapper

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
                     sparseness=sparval[1],  mycoption=1,  smooth=0.0, #z=-1/nclasses,
                     cthresh=250, its=3, #, nsamp=1 )
                     initializationList=initlist ) #, nsamp=1 )
  eanatmat<-imageListToMatrix( eanat$eigenanatomyimages, mask )
  rownames(eanatmat)<-eanatnames
  mycor<-rep(0,length(initlist))
  for ( i in 1:length(initlist) )
    {
    x<-newmat[1,]
    y<-eanatmat[i,]
#    ww<-which( x != 0 & y!= 0 )
#    mycor[i]<-cor.test( x, y , method="spearman" )$est
#    mycor[i]<-sqrt( sum( ( x[ww] - y[ww] )^2 ) )*(-1)
        mycor[i]<- mylocaldistfun( x, y )
#    mycor[i]<-robcosineSim( x, y )
    }
  names(mycor)<-eanatnames
  print(mycor)
#  pheatmap( rbind(mycor,mycor))
  fclass<-which.max(abs(mycor))
  mycor[fclass]<-0
  sclass<-which.max(abs(mycor))
  return( list(class=paste(eanatnames[c(fclass,sclass)]),
               patternimages=list(eanat$eig[[fclass]], eanat$eig[[sclass]]),
               featuretemplate=featuretemplate,
               eanatmat=eanatmat ) )
  }
if ( method == "sccan" & !is.na(eigsents) )
  {
  initlist<-list()
  if ( is.na(mask) ) {
    maskmat<-newmat*0
    maskmat[1,]<-1
    mask<-as.antsImage( maskmat )
  }
  nreps<-1
  eanatnames<-rep(as.character("A"),nclasses*nreps)
  ct<-1
  for ( i in 1:nclasses ) {
    vecimg<-antsImageClone( mask )
    vecimg[ mask == 1 ]<-featuretemplate[i,]
    vecimg[ mask == 1 ]<-eanatsparsify( featuretemplate[i,] , -0.25 )
    for (  nr in 1:nreps )
      {
      initlist<-lappend(initlist,vecimg)
      eanatnames[ct+nr-1]<-toString(classlabels[i])
      }
    ct<-ct+nreps
  }
  # build eigsent maps
  
  eanat<-sparseDecom2(list(exemplarmat,eigsents),
                      inmask=c(mask,NA), # z=-1/nclasses, 
                      nvecs=length(initlist),
                      sparseness=sparval,  mycoption=1,
                      smooth=0.0, cthresh=c(250,0), its=10, ell1=1.0,
                      initializationList=initlist )
  eanatmat<-imageListToMatrix( eanat$eig1, mask )
  rownames(eanatmat)<-eanatnames
  mycor<-rep(0,length(initlist))
  for ( i in 1:length(initlist) )
    {
    x<-newmat[1,]
    y<-eanatmat[i,]
#    ww<-which( x != 0 & y!= 0 )
#    mycor[i]<-cor.test( x, y , method="spearman" )$est
    ww<-which( abs(x)/max(abs(x)) > 1.e-6 & abs(y)/max(abs(y)) > 1.e-6 )
    mycor[i]<-sqrt( sum( ( x[ww] - y[ww] )^2 ) )*(-1)
#    mycor[i] <- mylocaldistfun( x, y )
    }
  names(mycor)<-eanatnames
  print(mycor)
#  pheatmap( rbind(mycor,mycor))
  fclass<-which.max((mycor))
  mycor[fclass]<-0
  sclass<-which.max((mycor))
  return( list(class=paste(eanatnames[c(fclass,sclass)]),
               patternimages=list(eanat$eig1[[fclass]], eanat$eig1[[sclass]]),
               featuretemplate=featuretemplate,
               eanatmat=eanatmat ) )
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

robcosineSim<-function (xin, yin, eps=0) 
{
    ww<-which( xin != 0 & yin!= 0 )
    ww<-which( abs(xin) > eps & abs(yin) > eps )
    x <- t(as.matrix(xin[ww]))
    y <- t(as.matrix(yin[ww]))
    return(as.numeric(1 - x %*% t(y)/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))))
}

tempclassrobcor <-function(  x, y )
  {
  ww<-which( x != 0 & y!= 0 )
  return( cor( x[ww] , y[ww] ) )
  }

overlapper <-function(  x, y , eps=0.1 )
{
ww<-which( abs(x)/max(abs(x)) > eps & abs(y)/max(abs(y)) > eps )
ww<-(cor( rank(x[ww]), rank(y[ww]) )) # ( abs(x) > eps & abs(y) > eps )
return( ww )
}
