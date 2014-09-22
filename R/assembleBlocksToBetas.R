assembleBlocksToBetas <- function( bmask, aalimg, labs , datadir, imagepostfix, assembledDesignOutPrefix, 
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
imat<-matrix()
for ( session in mysessions ) {
  localblocks<-sort( unique( dmat$blockNumb[ dmat$session == session ] ) )
  for ( block in localblocks ) {
    # reconstruct imagefn from dmat
    blocknum<-sprintf("%03d", block )
    blockpad<-paste(block,sep='')
    if ( as.numeric(block) < 1000 ) blockpad<-paste('0',block,sep='') 
    if ( as.numeric(block) < 100 )  blockpad<-paste('00',block,sep='') 
    if ( as.numeric(block) < 10 )   blockpad<-paste('000',block,sep='') 
    betafn<-paste(datadir,"/betas/",session,"_",blockpad,"_betas.mha",sep='')
    cat(paste(betafn,"*"))
    blockindices<-which( dmat$session == session & dmat$blockNumb == block )
    partsofblocktouse<-rep( TRUE, length(blockindices) )
    imagefn<-paste(imagedir,subject,"_",session,"_",blocknum,imagepostfix,sep='')
    if ( ! file.exists( betafn ) ) {
    if ( file.exists(imagefn) ) img<-antsImageRead( imagefn ,4)
    if ( isTRUE( all.equal( maskdim, dim(img)[1:3] ) ) ) # or identical(...)
    {
    if ( spatialsmooth > 0 ) {
        simg<-antsImageClone( img )
        SmoothImage(4,img,0.5,simg)
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
    if ( zscore ) {
        locmatmean<-colMeans( locmat )
        locmatsd<-apply( locmat, FUN=sd, MARGIN=2 )
        locmatsd[ locmatsd == 0 ]<-1
        locmat<-( locmat - locmatmean ) / locmatsd
    }
    if ( length( blockindices ) != nrow(locmat) ) stop(" length( blockindices ) != nrow(locmat) ")
    locmat<-locmat[partsofblocktouse,]
#### now go bold2betas
      btsc<-bold2betas( boldmatrix=locmat, 
        designmatrix=dmat[blockindices, 281:ncol(dmat) ], baseshift=0, verbose=F,
        blockNumb=rep(1,nrow(locmat)), maxnoisepreds=4, hrfBasis=hrf,
        hrfShifts=6, polydegree=4, selectionthresh=0.2 )
      antsImageWrite( as.antsImage( data.matrix( btsc$eventbetas ) ), betafn )
    } # beta exists
    if ( file.exists(imagefn) )
    {
    usedesignrow[blockindices]<-partsofblocktouse
    }
    } # identical check
    } # localblocks
  cat(paste(session,"done",dim(imat)[1],"by",dim(imat)[2],"assembled so far\n"))
} # session 
dmat<-dmat[usedesignrow,]
} # existence
return( list( dmat=dmat, usedesignrow=usedesignrow, subaal=subaal ) )
}
