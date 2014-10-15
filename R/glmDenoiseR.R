glmDenoiseR <- function( boldmatrix, designmatrixIn , hrfBasis=NA, hrfShifts=4, 
    selectionthresh=0.1, maxnoisepreds=1:12, collapsedesign=TRUE, 
    reestimatenoisepool=FALSE, debug=FALSE, polydegree=4 , crossvalidationgroups=4, 
    timevals=NA, runfactor=NA,  tr=1, baseshift=0, auxiliarynuisancevars=NA, 
    svdonallruns=FALSE, noisepoolfun=max, myintercept=0 )
{
nvox<-ncol(boldmatrix)
designmatrix<-as.matrix( designmatrixIn[,colMeans(abs(designmatrixIn))>0 ] )
groups<-crossvalidationgroups
if ( length(groups) == 1 ) {
  kfolds<-groups
  groups<-c()
  grouplength<-round(nrow(boldmatrix)/kfolds)-1
  for ( k in 1:kfolds ) groups<-c(groups,rep(k,grouplength))
  groups<-c( rep(1,nrow(boldmatrix)-length(groups)) , groups)
}

getnoisepool<-function( x, frac = selectionthresh ) {
  xord<-sort(x)
  l<-round(length(x)*frac)
  val<-xord[l]
  return( x < val & x < 0 )
}

crossvalidatedR2<-function( residmat, designmathrf, groups , noiseu=NA, p=NA, howmuchnoise ) {
  nvox<-ncol(residmat)
  kfo<-unique( groups )
  R2<-matrix(rep(0, nvox * length(kfo) ), nrow=length(kfo) )
  for ( k in kfo )
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
    for ( v in 1:nvox ) R2[k,v]<-100*( 1 -  sum( ( predmat[,v] - realmat[,v] )^2 ) / 
                                      sum(  (mean(realmat[,v]) - realmat[,v] )^2 )  )
    }
  # TODO write some kendrick like graphics that show the R2 
  # plotted over the bold image - need a mask as input
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
# make polynomial regressors per run / cv group
if ( all(is.na(timevals)) ) {
  timevals<-rep(0,nrow(designmatrix))
  for ( run in unique(groups)  ) {
    timeinds<-which( groups == run )
    timevals[ timeinds ]<-1:length(timeinds)
  }
}
p<-stats::poly( timevals ,degree=polydegree )
if ( all( !is.na(auxiliarynuisancevars) ) ) p<-cbind(  p, data.matrix(auxiliarynuisancevars) )
if ( all( !is.na(runfactor) ) ) p<-cbind(p,runfactor)
rawboldmat<-data.matrix(boldmatrix)
rawboldmatsd<-apply( rawboldmat , FUN=sd, MARGIN=2 )
rawboldmat[ , rawboldmatsd==0 ]<-rowMeans( rawboldmat[ , rawboldmatsd>0 ] )
svdboldmat<-rawboldmat
for ( run in unique(groups)  )
  {
  # FIXME - should clarify that this is done in the documentation
  # TODO - add option of single model
  timeinds<-( groups == run )
  svdboldmat[timeinds,]<-residuals( lm( rawboldmat[timeinds,] ~ myintercept + p[ timeinds, ]  ) )
  }
# FIXME - consider residualizing nuisance against design matrix
if (debug) print('lm')
# FIXME - factor out both HRF estimation approaches as functions
# FIXME - implement HRF library and just loop over library
if ( !all(is.na(hrfBasis)) ) { # use shifted basis functions
  if ( hrfShifts > 1 ) {
    fir<-finiteImpulseResponseDesignMatrix( designmatrix,
            n=hrfShifts, baseshift=baseshift )
  } else fir<-designmatrix
  for ( i in 1:ncol(fir) )
      fir[,i]<-conv( fir[,i]  , hrfBasis )[1:nrow(designmatrix)]
  mylm<-lm( svdboldmat  ~  fir )
  mylm<-bigLMStats( mylm, 0.01 )
  # here call crossvalidatedR2 to allow us to select best voxels by R2
  betas<-mylm$beta.t[1:hrfShifts,]
  if (debug) print('meanmax')
  meanmax<-function( x ) {  return( mean(sort((x),decreasing=T)[1:50]) ) }
  if ( hrfShifts <= 1 ) {
    # Old-school VCR format.
    betamax<-meanmax( betas )
  } else {
    betamax<-apply( (betas),FUN=meanmax,MARGIN=1)
  }
  betamax<-betamax/sum(abs(betamax))
  if ( debug ) print(betamax)
  hrf<-hrfBasis*0
  for ( i in 1:length(betamax) )
    {
    hrf<-hrf+shift(hrfBasis,baseshift+i-1)*betamax[i]
    }
} else { # use FIR / deconvolution
  # Q: What's the important difference with this new way?
  # After looking through it, looks like this is deconvolution rather than
  # convolving event onset with an assumed HRF. Is that right?
  # A: Yes
  ldes<-matrix(rowMeans(abs(designmatrix)),ncol=1)
  ldes<-ldes/ldes; ldes[is.nan(ldes)]<-0
  fir<-finiteImpulseResponseDesignMatrix( ldes,
            n=hrfShifts, baseshift=baseshift )
  mylm<-lm( rawboldmat  ~  fir + p )
  mylm<-bigLMStats( mylm, 0.01 )
  betablock<-mylm$beta.t[1:ncol(fir),]
  sumbetablock<-betablock[1:hrfShifts,]*0
  j<-1
  for ( i in 1:ncol(ldes) ) {
    sumbetablock<-sumbetablock+betablock[j:(j+hrfShifts-1),]
    j<-j+hrfShifts
  }
  # Q: Not sure what's being summed here. Is this summing the betas for all the HRF shifts tested?
  # If so, why? Or are you summing betas in an FIR model to get area-under-the-curve?
  # A: estimate the HRF from the "best fit" set of predictors. "best fit" defined by high beta values.
  betablock<-sumbetablock
  temp<-apply( (betablock) , FUN=sum, MARGIN=2)
  tempord<-sort(temp,decreasing=TRUE)
  bestvoxnum<-50
  # Q: Finding the voxels with the 50 highest betas?
  # A: Yes
  bestvoxels<-which( temp >= tempord[bestvoxnum]  )
  # Q: Deriving an HRF model from the voxels with the highest summed betas in an FIR model?
  # A: Yes
  # FIXME : define bestvoxels based on crossvalidatedR2 based on residualized data
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
designmatrix<-as.matrix( designmatrixIn[,colMeans(abs(designmatrixIn))>0 ] )
if ( collapsedesign ) designmatrix<-as.matrix( as.numeric( rowSums( designmatrix ) > 0 ) )
if (debug) print('hrf conv')
hrfdesignmat<-designmatrix
for ( i in 1:ncol(hrfdesignmat) )
  {
  hrfdesignmat[,i]<-conv( hrfdesignmat[,i]  , hrf )[1:nrow(hrfdesignmat)]
  }
R2base<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , p=NA )
R2base<-apply(R2base,FUN=noisepoolfun,MARGIN=2)
noisepool<-getnoisepool( R2base )
if ( max(maxnoisepreds) == 0 )  # TODO: denoise from polynomials by run
  {
  return(list( n=0, R2atBestN=NA, hrf=hrf, noisepool=noisepool, R2base=R2base, R2final=NA, hrfdesignmat=hrfdesignmat, noiseu=rep(1,nrow(hrfdesignmat)), polys=p ))
  }
print(paste("Noise pool has nvoxels=",sum(noisepool)))
# Step 5. [Calculate noise regressors using PCA on time-series of voxels
# in the noise pool] For each run, we extract the time-series of the
# voxels in the noise pool, project out the polynomial regressors from
# each time-series, normalize each time-series to unit length, and
# perform principal components analysis (PCA) (Behzadi et al., 2007;
# Bianciardi et al., 2009b). The resulting principal components
# constitute candidate noise regressors.
svdboldmat<-scale(svdboldmat) # z-score 
if ( svdonallruns ) {
  noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u
} else {
  for ( run in unique(groups)  ) {
    locmat<-rawboldmat[  groups == run ,noisepool]
    locmat<-scale( residuals( lm( locmat ~ myintercept + p[ groups == run, ] ) ) )
    locsvd<-svd( locmat, nv=0, nu=max(maxnoisepreds) ) # TODO: should this be pca?
    if ( run == unique(groups)[1]  ) noiseu<-locsvd$u else noiseu<-rbind( noiseu, locsvd$u )
  }
}
R2summary<-rep(0,length(maxnoisepreds))
ct<-1
for ( i in maxnoisepreds )
  {
  svdboldmat<-residuals(lm(rawboldmat~myintercept+p+noiseu[,1:i]))
  # The noise pool could change as one finds a better model so allow user to optionally re-estimate it
  if ( reestimatenoisepool )
    {
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u
    R2<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , noiseu, howmuchnoise=i, p=NA  )
    R2<-apply(R2,FUN=noisepoolfun,MARGIN=2)
    noisepool<-getnoisepool( R2 )
    noiseu<-svd( svdboldmat[,noisepool], nv=0, nu=max(maxnoisepreds) )$u # FIXME: Re-estimate by run
    }
  R2<-crossvalidatedR2(  svdboldmat, hrfdesignmat, groups , noiseu=NA, howmuchnoise=i, p=NA  )
  if ( reestimatenoisepool ) {
      R2min<-apply(R2,FUN=noisepoolfun,MARGIN=2)
      noisepool<-getnoisepool( R2min )
  }
  R2max<-apply(R2,FUN=max,MARGIN=2)
  if ( ct == 1 ) R2perNoiseLevel<-R2max else R2perNoiseLevel<-cbind(R2perNoiseLevel,R2max)
  R2pos<-R2max[ R2max > 0 ]
  R2summary[ct]<-median(R2pos)
  print(paste("NoiseU:",i,"MeanRSqrd",  R2summary[ct] ))
  ct<-ct+1
  }
scl<-0.95
if (max(R2summary)<0) scl<-1.05
bestn<-maxnoisepreds[which( R2summary > scl*max(R2summary) )[1]]
hrf<-hrf/max(hrf)
for ( run in unique(groups)  ) {
  locmat<-rawboldmat[  groups == run ,noisepool]
  locmat<-scale( residuals( lm( locmat ~ myintercept + noiseu[ groups == run, 1:bestn ]
                              + p[ groups == run, ] ) ) )
  if ( run == unique(groups)[1]  ) denoisedBold<-locmat else denoisedBold<-rbind( denoisedBold, locmat )
}
return(list( denoisedBold=denoisedBold, n=bestn, R2atBestN=R2summary[bestn],
            hrf=hrf, noisepool=noisepool, R2base=R2base, R2final=R2perNoiseLevel,
            hrfdesignmat=hrfdesignmat, noiseu=noiseu[,1:bestn], polys=p ))
}
