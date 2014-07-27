glmDenoiseR <- function( boldmatrix, designmatrixIn , hrfbasislength=50, kfolds=4, whichbase = NA, selectionthresh=0.25, maxnoisepreds=12, collapsedesign=TRUE , reestimatenoisepool=FALSE, debug=FALSE, polydegree=6 )
{
nvox<-ncol(boldmatrix)
designmatrix<-as.matrix( designmatrixIn[,colMeans(designmatrixIn)>0 ] )
groups<-c()
grouplength<-round(nrow(boldmatrix)/kfolds)-1
for ( k in 1:kfolds ) groups<-c(groups,rep(k,grouplength))

getnoisepool<-function( x, frac = selectionthresh ) {
  xord<-sort(x)
  l<-round(length(x)*frac)
  val<-xord[l]
  return( x < val )
}

crossvalidatedR2<-function( residmat, designmathrf, groups , noiseu=NA, p=NA ) {
  nvox<-ncol(residmat)
  kfo<-max(groups)
  R2<-matrix(rep(0, nvox * kfo ),nrow=kfo)
  for ( k in 1:kfo )
    {
    selector <- groups!=k
    mydf<-data.frame( designmathrf[selector,] )
    if ( ! all( is.na(noiseu) ) )
      mydf<-data.frame( mydf, noiseu[selector,1:i] )  
    if ( ! all( is.na(p) ) )
      mydf<-data.frame( mydf, p[selector,] )
    mylm1<-lm( residmat[selector,]   ~  . , data=mydf )
    selector <- groups==k
    mydf<-data.frame( designmathrf[selector,] )
    if ( ! all( is.na(noiseu) ) )
      mydf<-data.frame( mydf, noiseu[selector,1:i] )  
    if ( ! all( is.na(p) ) )
      mydf<-data.frame( mydf, p[selector,] )
    predmat<-predict(mylm1,newdata=mydf)
    realmat<-residmat[selector,]
    for ( v in 1:nvox ) R2[k,v]<-100*( 1 -  sum( ( predmat[,v] - realmat[,v] )^2 ) / sum(  (mean(realmat[,v]) - realmat[,v] )^2 )  )
#    plot( ts( rowMeans(realmat) ) )
    }
  return(R2)
}

mygamma <- function(x, a1 = 6.,   a2 = 12., b1 = 0.9, b2 = 0.99, cc = 0.35) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- (x/d1)^a1
    c2 <- cc * (x/d2)^a2
    res <- c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)
    res
  }

# start with a reasonable default set of bases 
b1<-mygamma(c(1:hrfbasislength),  1,  5 ,0.9,0.9,0.05) 
b2<-mygamma(c(1:hrfbasislength),  5, 10 ,0.9,0.9,0.05) 
b3<-mygamma(c(1:hrfbasislength), 10, 15 ,0.9,0.9,0.05) 
b4<-mygamma(c(1:hrfbasislength), 15, 20 ,0.9,0.9,0.05) 
b5<-mygamma(c(1:hrfbasislength), 20, 25 ,0.9,0.9,0.05) 
b6<-mygamma(c(1:hrfbasislength), 25, 30 ,0.9,0.9,0.05)
basismat<-cbind( b2, b3, b4, b5, b6 )
if ( ! all(is.na(whichbase)) ) basismat<-basismat[,whichbase[ whichbase <= ncol(basismat)]]
basismat<-as.matrix( basismat )
#################################################
####### convolutions on basis image matrix ######
nbase<-ncol(basismat) # number of basis functions
if ( ncol(designmatrix) == 1 ) {
  designmatrix<-rep( designmatrix, nbase )
  designmatrix<-t(matrix( designmatrix, nrow=nbase ))
}
if ( ncol(designmatrix) > 1 ) designmatrixext<-interleaveMatrixWithItself( designmatrix, nbase )
k<-1
if (debug) print('init conv')
for ( i in 1:ncol(designmatrixext) )
  {
  if ( k > ncol(basismat) ) k<-1
  designmatrixext[,i]<-conv( designmatrixext[,i]  , basismat[,k] )[1:nrow(designmatrixext)]
  if (debug) print(paste(k,i))
  k<-k+1
  }
if (debug) print('init conv done')
# overall description of the method
# 0. estimate hrf
# 1. regressors include: design + trends + noise-pool
# 2. find noise-pool by initial cross-validation without noise regressors
# 3. cross-validate predictions using different numbers of noise regressors
# 4. select best n for predictors from noise pool
# 5. return the noise mask and the value for n
# make regressors
p<-stats::poly( 1:nrow(designmatrixext) ,degree=polydegree )
rawboldmat<-data.matrix(boldmatrix)
svdboldmat<-residuals(lm(rawboldmat~p))
if (debug) print('lm')
mylm<-lm( svdboldmat  ~  designmatrixext  )
mylm<-bigLMStats( mylm, 0.01 )
betas<-mylm$beta
if (debug) print('meanmax')
meanmax<-function( x ) {  return( mean(sort(x,decreasing=T)[1:50]) ) }
betamax<-apply(abs(betas),FUN=meanmax,MARGIN=1)
betamax<-betamax/sum(betamax)
if ( debug ) print(betamax)
hrf<-basismat[,1]*0
k<-1
for ( i in 1:length(betamax) )
  {
  if ( k > ncol(basismat) ) k<-1
  hrf<-hrf+basismat[, k ]*betamax[i]
  k<-k+1
  }
hrf<-hrf/max(hrf)
if ( debug ) plot( ts( hrf ) )
################### now redo some work w/new hrf
if ( collapsedesign ) designmatrix<-as.matrix( rowSums( designmatrix ) )
designmatrixext<-designmatrix
if (debug) print('hrf conv')
for ( i in 1:ncol(designmatrixext) )
  {
  designmatrixext[,i]<-conv( designmatrix[,i]  , hrf )[1:nrow(designmatrix)]
  }
R2base<-crossvalidatedR2(  svdboldmat, designmatrixext, groups , p=NA )
R2base<-apply(R2base,FUN=min,MARGIN=2)
noisepool<-getnoisepool( R2base )
if ( all( noisepool==TRUE ) )
  {
  print("all voxels meet your pvalthresh - try increasing the value")
  return(NA)
  } 
if ( all( noisepool==FALSE ) )
  {
  print("zero voxels meet your pvalthresh - try decreasing the value")
  return(NA)
  } else print(paste("Noise pool has nvoxels=",sum(noisepool)))
svdboldmat<-scale(svdboldmat)
noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=maxnoisepreds )$u
R2summary<-rep(0,maxnoisepreds)
for ( i in 1:maxnoisepreds )
  {
  if ( reestimatenoisepool )
    {
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=maxnoisepreds )$u
    R2<-crossvalidatedR2(  svdboldmat, designmatrixext, groups , noiseu, p=NA  )
    R2<-apply(R2,FUN=median,MARGIN=2)
    noisepool<-getnoisepool( R2 )
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=maxnoisepreds )$u
    }
  R2<-crossvalidatedR2(  svdboldmat, designmatrixext, groups , noiseu, p=NA  )
  R2<-apply(R2,FUN=min,MARGIN=2)
  if ( reestimatenoisepool ) noisepool<-getnoisepool( R2 )
  if ( i == 1 ) R2perNoiseLevel<-R2 else R2perNoiseLevel<-cbind(R2perNoiseLevel,R2)
  print(paste("Noise pool has nvoxels=",sum(noisepool)))
  R2summary[i]<-mean(R2)
  print(paste("NoiseU:",i,"MeanRSqrd",  R2summary[i] ))
  }
scl<-0.95
if (max(R2summary)<0) scl<-1.05
bestn<-which( R2summary > scl*max(R2summary) )[1]
hrfdesignmat<-designmatrixIn
for ( i in 1:ncol(hrfdesignmat) )
  {
  hrfdesignmat[,i]<-conv( hrfdesignmat[,i]  , hrf )[1:nrow(hrfdesignmat)]
  }
# glmdenoisedataframe<-data.frame(  hrfdesignmat=hrfdesignmat, noiseu=noiseu[,1:bestn], polys=p )
# return( glmdenoisedataframe )
return(list( n=bestn, hrf=hrf, noisepool=noisepool, R2base=R2base, R2final=R2perNoiseLevel, hrfdesignmat=hrfdesignmat, noiseu=noiseu[,1:bestn], polys=p ))
}


#
# mylm<-lm( data.matrix(boldmatrix)  ~  designmatrixext + p )
# biglm1<-bigLMStats( mylm, 0.01 )
# collapsep<-biglm1$beta.pval[1,]
# if ( ncol(designmatrixext) > 1 )
#  collapsep<-apply( biglm1$beta.pval[1:ncol(designmatrixext),], FUN=mean , MARGIN=2)
# noisepool<-collapsep > pvalthresh
#
# mylm<-lm( data.matrix(boldmatrix)  ~  designmatrixext + p + noiseu[,1:which.min(bestnoiselev)]  )
# biglm2<-bigLMStats( mylm, 0.01 )
# return( list( biglm1=biglm1, biglm2=biglm2 ) )
#
