assembleSessionBlocks <- function( bmask, aalimg, labs , datadir, imagepostfix, assembledDesignOutPrefix,
    assembledImageOutPrefix, dmat, usedesignrow, hrf, ncompcor=6, zscore=TRUE, spatialsmooth=0 )
{
# print("#########assemble image blocks, potentially event-specifically#########")
maskdim<-dim( bmask )
aalmask<-antsImageClone( aalimg )
aalmask[ aalmask > 1 ]<-1
subaal<-antsImageClone( aalimg )
subaal[ aalimg > 0 ]<-0
for ( lab in labs ) subaal[ aalimg == lab ]<-1
print(paste("assemble blocks for",afn))
if ( all( is.na(as.numeric(imat)) ) ) { # assembly begin ......
mysessions<-sort( unique( dmat$session) )
for ( session in mysessions ) {
  imat<-matrix()
  betafn<-paste(datadir,"/sessions/",session,"session.mha",sep='')
  cat(paste(betafn,"*"))
  if ( ! file.exists( betafn ) ) {
  localblocks<-sort( unique( dmat$blockNumb[ dmat$session == session ] ) )
  for ( block in localblocks ) {
    # reconstruct imagefn from dmat
    blocknum<-sprintf("%03d", block )
    blockindices<-which( dmat$session == session & dmat$blockNumb == block )
    partsofblocktouse<-rep( TRUE, length(blockindices) )
    imagefn<-paste(imagedir,subject,"_",session,"_",blocknum,imagepostfix,sep='')
    if ( file.exists(imagefn) ) img<-antsImageRead( imagefn ,4)
    if ( isTRUE( all.equal( maskdim, dim(img)[1:3] ) ) ) # or identical(...)
    {
    if ( spatialsmooth > 0 ) {
        simg<-antsImageClone( img )
        SmoothImage(4,img,spatialsmooth,simg)
        img<-simg
    }
    locmat<-timeseries2matrix( img , subaal )
    locmatmean<-apply( locmat[(throwaway+1):nrow(locmat),], FUN=mean, MARGIN=2 )
    for( j in 1:throwaway ) locmat[j,]<-locmatmean
    if ( ncompcor > 0 ) {
        locmatfull<-timeseries2matrix( img , bmask )
        locmatfullmean<-apply( locmatfull[(throwaway+1):nrow(locmat),], FUN=mean, MARGIN=2 )
        for( j in 1:throwaway ) locmatfull[j,]<-locmatfullmean
        if ( ! exists("mycompcorv") ) {
          mycompcorv<-compcor( locmatfull, ncompcor=ncompcor, variance_extreme = compcorvarval, returnv=TRUE )
          highvarmatinds<-compcor( locmatfull, ncompcor=ncompcor,
                                  variance_extreme = compcorvarval, returnhighvarmatinds=TRUE )
        }
#        mycompcor<-scale( locmatfull[,highvarmatinds]  %*% mycompcorv )
        mycompcor<-compcor( locmatfull, ncompcor, variance_extreme = compcorvarval )
        locmat<-residuals(lm(locmat~0+mycompcor))
    }
    if ( length( blockindices ) != nrow(locmat) ) stop(" length( blockindices ) != nrow(locmat) ")
    locmat<-locmat[partsofblocktouse,]
    usedesignrow[blockindices]<-partsofblocktouse
    if ( dim(imat)[1] == 1 ) imat<-locmat
    if ( dim(imat)[1] > 1  ) imat<-rbind(imat,locmat)
    } # identical check
    cat("*")
    } # block in localblocks
    #### now do glmdenoise on imat
    # cat("glmdenoise")
    # dd2<-glmDenoiseR( imat, bb$desmat[,1:4], hrfBasis=hrf, hrfShifts = 4 ,
    #    crossvalidationgroups=runs,  maxnoisepreds=4 , selectionthresh=0.1,
    #    collapsedesign=T, polydegree=4, verbose=T )
    } # beta exists
  cat("write")
  if ( zscore ) {
    imatmean<-colMeans( imat )
    imatsd<-apply( imat, FUN=sd, MARGIN=2 )
    imatsd[ imatsd == 0 ]<-1
    imat<-( imat - imatmean ) / imatsd
  }
  antsImageWrite( as.antsImage(imat), betafn  )
  cat(paste(session,"done",dim(imat)[1],"by",dim(imat)[2],"assembled so far\n"))
} # session
dmat<-dmat[usedesignrow,]
} # existence
return( list( imat=imat, dmat=dmat, usedesignrow=usedesignrow, subaal=subaal ) )
}
