#' compute betas from image blocks of bold
#' 
#' Computes event-wise betas for each image block and writes to an image
#' 
#' 
#' @param bmask brain mask - compcor is done within this mask
#' @param aalimg a label image - does not have to be aal
#' @param labs which labels of the above label image to use e.g. c(1,2,3)
#' @param datadir the input data directory
#' @param imagepostfix the image blocks postfix e.g. motioncorr.nii.gz
#' @param assembledBlocksOutPrefix output file name for assembled design file -
#' will be read if it already exists (rather than reassembled)
#' @param assembledImageOutPrefix output file name for assembled bold file -
#' will be read if it already exists (rather than reassembled)
#' @param usedesignrow organizational boolean tracking which parts of the
#' design file have concurrent BOLD images, useful for sanity checks
#' @param imat the current estimate of the assembled image matrix
#' @param zscore should we zscore each voxel?
#' @param spatialsmoothing scalar value greater than or equal to zero that
#' controls voxel-wise smoothing in 4D
#' @return list is output containing organizational variables for design matrix
#' @author Avants BB, Phillips JS
#' @examples
#' 
#' \dontrun{
#' hrf<-ts( hemodynamicRF( 30, onsets=2, durations=1, rt=0.5,cc=0.1,a1=6,a2=12,b1=0.6, b2=0.6 ) )
#' assembly2<-assembleSessionBlocks( bmask, aalimg, labs, datadir, imagepostfix,  dfn, afn, dmat, usedesignrow, imat )
#' }
#' 
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
