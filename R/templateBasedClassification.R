templateBasedClassification <- function( exemplarmat, labels, newmat,
                                        method="corr", mask=NA, eigsents=NA , sparval=c(-0.2,-0.9) )
{
mylocaldistfun <- tempclassrobcor
mylocaldistfun <- robcosineSim
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
rownames(featuretemplate)<-classlabels[1:nclasses]
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
  nreps<-1
  eanatnames<-rep(as.character("A"),nclasses*nreps)
  ct<-1
  for ( i in 1:nclasses ) {
    vecimg<-antsImageClone( mask )
    initf<-featuretemplate[i,] + rnorm( length(featuretemplate[i,]), mean=0, sd=sd(featuretemplate[i,]) )*0.1
    vecimg[ mask == 1 ]<-eanatsparsify( initf , sparval[1] )
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
    ww<-which( abs(x)/max(abs(x)) > 1.e-6 & abs(y)/max(abs(y)) > 1.e-6 )
    mycor[i]<-sqrt( sum( ( x[ww] - y[ww] )^2 ) )*(-1)
#    mycor[i]<- mylocaldistfun( x, y )
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
  initlist2<-list()
  if ( is.na(mask) ) {
    maskmat<-newmat*0
    maskmat[1,]<-1
    mask<-as.antsImage( maskmat )
  }
  sentmaskmat<-eigsents*0
  sentmaskmat[1,]<-1
  sentmask<-as.antsImage( sentmaskmat )
  nreps<-1
  eanatnames<-rep(as.character("A"),nclasses*nreps)

#  mylm<-bigLMStats( lm( exemplarmat ~ es[,1] ) )
  ct<-1
  for ( i in 1:nclasses ) {
    vecimg<-antsImageClone( mask )
    initf<-featuretemplate[i,] + rnorm( length(featuretemplate[i,]), mean=0, sd=sd(featuretemplate[i,]) )*0.0
    vecimg[ mask == 1 ]<-initf # eanatsparsify( initf , abs(sparval[1]) )
    for (  nr in 1:nreps )
      {
      initlist<-lappend(initlist,vecimg)
      eanatnames[ct+nr-1]<-toString(classlabels[i])
      }
    sentimg<-antsImageClone( sentmask )
    if ( (i %% nclasses) == 1 ) sentimg[ sentmask == 1 ]<-c( rep(1,(ncol(eigsents)/2)), rep(0,(ncol(eigsents)/2)) )
    if ( (i %% nclasses) == 0 ) sentimg[ sentmask == 1 ]<-c( rep(0,(ncol(eigsents)/2)), rep(1,(ncol(eigsents)/2)) )
    if ( (i %% nclasses) == 1 ) sentimg[ sentmask == 1 ]<-c( rnorm((ncol(eigsents)/2), 1,0.1), rnorm((ncol(eigsents)/2),  0,0.1) )
    if ( (i %% nclasses) == 0 ) sentimg[ sentmask == 1 ]<-c( rnorm((ncol(eigsents)/2), 0,0.1), rnorm((ncol(eigsents)/2),  1,0.1 ) )
    for (  nr in 1:nreps )
      {
      initlist2<-lappend(initlist2,sentimg)
      }
    ct<-ct+nreps
  }
  # build eigsent maps
  antsImageWrite( as.antsImage( exemplarmat ), 'exemplarmat.mha' )
  antsImageWrite( as.antsImage( eigsents )   , 'eigsents.mha' )
  antsImageWrite( mask        , 'exmask.mha' )
  antsImageWrite( sentmask    , 'semask.mha' )
  eanat<-sparseDecom2( list((exemplarmat),eigsents),
                       inmask=c(mask,NA), #  z=-1/nclasses, 
                       nvecs=length(initlist), perms=0,
                       sparseness=sparval,  mycoption=1,
                       smooth=0.0, cthresh=c(250,0), its=25, ell1=10 ) # ) # ,
#                       initializationList=initlist,
#                       initializationList2=initlist2 )
  eanatmat<-imageListToMatrix( eanat$eig1, mask )

  # measure similarity of eanatmat to featuretemplate
  sentpred<-0
  # sentpred<-imageListToMatrix(eanat$eig2, sentmask )
  rownames(eanatmat)<-eanatnames
  mycor<-rep(0,length(initlist))
  e1<-0
  e2<-0
  e1rows<-which( labels ==  rownames(eanatmat)[1]  ) 
  e2rows<-which( labels ==  rownames(eanatmat)[2]  ) 
  y1<-eanatmat[1,]    
  y2<-eanatmat[2,]
  y1b<-abs(y1)/max(abs(y1))
  y2b<-abs(y2)/max(abs(y2))
  nz1<-which( y1b > 1.e-3 )
  nz2<-which( y2b > 1.e-3 )
  nzall<-unique( c(nz1,nz2) )
#  print(paste("11", cor.test(eanatmat[1,nz1], featuretemplate[1,nz1] )$est ))
#  print(paste("22", cor.test(eanatmat[2,nz2], featuretemplate[2,nz2] )$est ))
#  print(paste("21", cor.test(eanatmat[2,nzall], featuretemplate[2,nzall] )$est ))
#  print(paste("12", cor.test(eanatmat[1,nzall], featuretemplate[1,nzall] )$est ))
#  print(paste("t1", cor.test(eanatmat[1,nz1], newmat[1,nz1] )$est ))
#  print(paste("t2", cor.test(eanatmat[2,nz2], newmat[1,nz2] )$est ))
#  print(paste("t1b", robcosineSim(eanatmat[1,nz1], newmat[1,nz1] ) ))
#  print(paste("t2b", robcosineSim(eanatmat[2,nz2], newmat[1,nz2] ) ))
  t1<-0#(cor.test(eanatmat[1,nz1], newmat[1,nz1] )$est)*(-1)
  t2<-0# cor.test(eanatmat[2,nz2], newmat[1,nz2] )$est
  pp1<-as.matrix( exemplarmat )  %*%   t(eanatmat)
  pp2<-newmat  %*%   t(eanatmat)
  mydf1<-data.frame( labels=rownames(exemplarmat), bold=pp1 )
  mydf2<-data.frame( bold=pp2 )
  mdl<-svm( labels ~ . , mydf1 )
  pp<-( predict( mdl, newdata=mydf2 ) )
  if ( abs(t1) > abs(t2) )  robcospred<-rownames(featuretemplate)[1] else   robcospred<-rownames(featuretemplate)[2]
#  print(paste("robcospred",robcospred))

  print(paste("L(w)",length(nz1),length(nz2),length(nzall)))
  mydf<-data.frame( feats=newmat[,nzall] )
  for ( i in e1rows )
    {
    x<-newmat[1,]
    z<-exemplarmat[i,]
#    ww<-which( abs(x)/max(abs(x)) > 1.e-2 & abs(y)/max(abs(y)) > 1.e-2 )
    e1<-e1+sqrt( sum( ( x[nz1]/mean(abs(x[nz1])) - z[nz1]/mean(abs(z[nz1])) )^2 ) )*(-1)/length(nz1)
    }
  for ( i in e2rows )
    {
    x<-newmat[1,]
    z<-exemplarmat[i,]
    y<-eanatmat[2,]
#    ww<-which( abs(x)/max(abs(x)) > 1.e-2 & abs(y)/max(abs(y)) > 1.e-2 )
    e2<-e2+sqrt( sum( ( x[nz2]/mean(abs(x[nz2])) - z[nz2]/mean(abs(z[nz2])) )^2 ) )*(-1)/length(nz2)
    }
  mycor[1]<-e1/length(e1rows)
  mycor[2]<-e2/length(e2rows)
  names(mycor)<-eanatnames
#  pheatmap( rbind(mycor,mycor))
  fclass<-which.max((mycor))
  mycor[fclass]<-0
  sclass<-which.max((mycor))
  # paste(eanatnames[c(fclass,sclass)])
  return( list(class=pp, svmclass=robcospred,  pp1=pp1, pp2=pp2, 
               patternimages=list(eanat$eig1[[fclass]], eanat$eig1[[sclass]]),
               featuretemplate=featuretemplate,
               eanatmat=eanatmat, sentpred=sentpred, eanat=eanat ) )
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
ww<-abs(cor( rank(x[ww]), rank(y[ww]) )) # ( abs(x) > eps & abs(y) > eps )
return( ww )
}
