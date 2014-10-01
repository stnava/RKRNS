glmDenoiseR <- function( boldmatrix, designmatrixIn , hrfBasis=NA, hrfShifts=4, selectionthresh=0.25, maxnoisepreds=1:12, collapsedesign=TRUE , reestimatenoisepool=FALSE, debug=FALSE, polydegree=6 , crossvalidationgroups=4, timevals=NA, runfactor=NA,  tr=1, baseshift=0, auxiliarynuisancevars=NA )
{
nvox<-ncol(boldmatrix)
designmatrix<-as.matrix( designmatrixIn[,colMeans(designmatrixIn)>0 ] )
groups<-crossvalidationgroups
if ( length(groups) == 1 ) {
  kfolds<-groups
  groups<-c()
  grouplength<-round(nrow(boldmatrix)/kfolds)-1
  for ( k in 1:kfolds ) groups<-c(groups,rep(k,grouplength))
}

getnoisepool<-function( x, frac = selectionthresh ) {
  xord<-sort(x)
  l<-round(length(x)*frac)
  val<-xord[l]
  return( x < val )
}

crossvalidatedR2<-function( residmat, designmathrf, groups , noiseu=NA, p=NA, howmuchnoise ) {
  nvox<-ncol(residmat)
  kfo<-max(groups)
  R2<-matrix(rep(0, nvox * kfo ),nrow=kfo)
  for ( k in 1:kfo )
    {
    selector <- groups!=k
    mydf<-data.frame( designmathrf[selector,] )
    if ( ! all( is.na(noiseu) ) )
      mydf<-data.frame( mydf, noiseu[selector,1:howmuchnoise] )  
    if ( ! all( is.na(p) ) )
      mydf<-data.frame( mydf, p[selector,] )
    mylm1<-lm( residmat[selector,]   ~  . , data=mydf )
    selector <- groups==k
    mydf<-data.frame( designmathrf[selector,] )
    if ( ! all( is.na(noiseu) ) )
      mydf<-data.frame( mydf, noiseu[selector,1:howmuchnoise] )  
    if ( ! all( is.na(p) ) )
      mydf<-data.frame( mydf, p[selector,] )
    predmat<-predict(mylm1,newdata=mydf)
    realmat<-residmat[selector,]
    for ( v in 1:nvox ) R2[k,v]<-100*( 1 -  sum( ( predmat[,v] - realmat[,v] )^2 ) / sum(  (mean(realmat[,v]) - realmat[,v] )^2 )  )
    }
  return(R2)
}

#################################################
# overall description of the method
# 0. estimate hrf
# 1. regressors include: design + trends + noise-pool
# 2. find noise-pool by initial cross-validation without noise regressors
# 3. cross-validate predictions using different numbers of noise regressors
# 4. select best n for predictors from noise pool
# 5. return the noise mask and the value for n
# make regressors
if ( all(is.na(timevals)) ) timevals<-1:nrow(designmatrix)
p<-stats::poly( timevals ,degree=polydegree )
if ( all( !is.na(auxiliarynuisancevars) ) ) p<-cbind(  p, data.matrix(auxiliarynuisancevars) ) 
if ( all( !is.na(runfactor) ) ) p<-cbind(p,runfactor)
rawboldmat<-data.matrix(boldmatrix)
rawboldmatsd<-apply( rawboldmat , FUN=sd, MARGIN=2 )
rawboldmat[ , rawboldmatsd==0 ]<-rowMeans( rawboldmat[ , rawboldmatsd>0 ] )
svdboldmat<-residuals(lm(rawboldmat~0+p))
if (debug) print('lm')
if ( !all(is.na(hrfBasis)) ) { # use shifted basis functions
  if ( hrfShifts > 1 ) {
    fir<-finiteImpuleResponseDesignMatrix( designmatrix,
            n=hrfShifts, baseshift=baseshift )
  } else fir<-designmatrix
  for ( i in 1:ncol(fir) )
      fir[,i]<-conv( fir[,i]  , hrfBasis )[1:nrow(designmatrix)]
  mylm<-lm( svdboldmat  ~  fir )
  mylm<-bigLMStats( mylm, 0.01 )
  betas<-mylm$beta.t[1:hrfShifts,]
  if (debug) print('meanmax')
  meanmax<-function( x ) {  return( mean(sort((x),decreasing=T)[1:50]) ) }
  if ( hrfShifts <= 1 ) {
      betamax<-meanmax( betas )
  } else {
      betamax<-apply( (betas),FUN=meanmax,MARGIN=1)
  }
  betamax<-betamax/sum(abs(betamax))
  if ( debug ) print(betamax)
  hrf<-hrfBasis*0
  k<-1
  for ( i in 1:length(betamax) )
    {
    hrf<-hrf+shift(hrfBasis,baseshift+i-1)*betamax[i]
    k<-k+1
    }
} else { # new way below
  ldes<-matrix(rowMeans(designmatrix),ncol=1)
  fir<-finiteImpuleResponseDesignMatrix( ldes,
            n=hrfShifts, baseshift=baseshift )
  mylm<-lm( rawboldmat  ~  fir + p )
  mylm<-bigLMStats( mylm, 0.01 )
  betablock<-mylm$beta.t[1:ncol(fir),]
  sumbetablock<-betablock[1:hrfShifts,]*0
  j<-1
  for ( i in 1:ncol(ldes) ) {
    sumbetablock<-betablock[j:(j+hrfShifts-1),]
    j<-j+hrfShifts
  }
  betablock<-sumbetablock
  temp<-apply( (betablock) , FUN=sum, MARGIN=2)
  tempord<-sort(temp,decreasing=TRUE)
  bestvoxnum<-50
  bestvoxels<-which( temp > tempord[bestvoxnum]  )
  hrf<-rowSums( (betablock[,bestvoxels] ) )
  meanhrfval<-mean(hrf)
  mxdf<-abs(max(hrf)-meanhrfval)
  mndf<-abs(min(hrf)-meanhrfval)
  if ( mndf > mxdf  ) hrf<-hrf*(-1)
  if ( abs(min(hrf)) > max(hrf) ) hrf<-hrf*(-1)
#  hrf<-data.frame(stl(ts(hrf, frequency = 4),"per")$time.series)$trend
}
hrf<-hrf/max(hrf)
if ( debug ) plot( ts( hrf ) )
################### now redo some work w/new hrf
# reset designmatrix
designmatrix<-as.matrix( designmatrixIn[,colMeans(designmatrixIn)>0 ] )
if ( collapsedesign ) designmatrix<-as.matrix( as.numeric( rowSums( designmatrix ) > 0 ) )
if (debug) print('hrf conv')
hrfdesignmat<-designmatrix
for ( i in 1:ncol(hrfdesignmat) )
  {
  hrfdesignmat[,i]<-conv( hrfdesignmat[,i]  , hrf )[1:nrow(hrfdesignmat)]
  }

R2base<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , p=NA )
R2base<-apply(R2base,FUN=min,MARGIN=2)
noisepool<-getnoisepool( R2base )
if ( max(maxnoisepreds) == 0 )
  {
  return(list( n=0, R2atBestN=NA, hrf=hrf, noisepool=noisepool, R2base=R2base, R2final=NA, hrfdesignmat=hrfdesignmat, noiseu=rep(1,nrow(hrfdesignmat)), polys=p ))
  }
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
noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u
R2summary<-rep(0,length(maxnoisepreds))
ct<-1
for ( i in maxnoisepreds )
  {
  svdboldmat<-residuals(lm(rawboldmat~0+p+noiseu[,1:i]))
  if ( reestimatenoisepool )
    {
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u
    R2<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , noiseu, howmuchnoise=i, p=NA  )
    R2<-apply(R2,FUN=min,MARGIN=2)
    noisepool<-getnoisepool( R2 )
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u
    }
  R2<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , noiseu=NA, howmuchnoise=i, p=NA  )
  if ( reestimatenoisepool ) {
      R2min<-apply(R2,FUN=min,MARGIN=2)
      noisepool<-getnoisepool( R2min )
  }
  R2max<-apply(R2,FUN=max,MARGIN=2)
  if ( ct == 1 ) R2perNoiseLevel<-R2max else R2perNoiseLevel<-cbind(R2perNoiseLevel,R2max)
  # print(paste("Noise pool has nvoxels=",sum(noisepool)))
  R2pos<-R2max[ R2max > 0 ]
  R2summary[ct]<-median(R2pos)
  print(paste("NoiseU:",i,"MeanRSqrd",  R2summary[ct] ))
  ct<-ct+1
  }
scl<-0.95
if (max(R2summary)<0) scl<-1.05
bestn<-maxnoisepreds[which( R2summary > scl*max(R2summary) )[1]]
hrf<-hrf/max(hrf)
return(list( n=bestn, R2atBestN=R2summary[bestn], hrf=hrf, noisepool=noisepool, R2base=R2base, R2final=R2perNoiseLevel, hrfdesignmat=hrfdesignmat, noiseu=noiseu[,1:bestn], polys=p ))
}
