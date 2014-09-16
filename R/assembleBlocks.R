assembleBlocks <- function( bmask, aalimg, labs , datadir, imagepostfix, assembledDesignOutPrefix,
    assembledImageOutPrefix, dmat, usedesignrow, imat=NA, ncompcor=6, zscore=TRUE, spatialsmooth=0 )
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
    imagefn<-paste(imagedir,subject,"_",session,"_",blocknum,imagepostfix,sep='')
    if ( file.exists(imagefn) ) img<-antsImageRead( imagefn ,4)
    if ( isTRUE( all.equal( maskdim, dim(img)[1:3] ) ) ) # or identical(...)
    {
    if ( spatialsmooth > 0 ) SmoothImage(4,img,0.5,img)
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
    ###########################################################################################
    # partsofblocktouse IS WHERE YOU MIGHT USE SOMETHING DIFFERENT TO FILTER BOTH DMAT & IMAT #
    # ... e.g. selected events                                                                #
    ###########################################################################################
    blockindices<-which( dmat$session == session & dmat$blockNumb == block )
    if ( length( blockindices ) != nrow(locmat) ) stop(" length( blockindices ) != nrow(locmat) ")
    partsofblocktouse<-rep( TRUE, length(blockindices) )
    locmat<-locmat[partsofblocktouse,]
    usedesignrow[blockindices]<-partsofblocktouse
    if ( nrow(imat)==1 & ncol(imat)==1 ) imat<-locmat else imat<-rbind( imat, locmat )
    }
    cat(paste(session,block,"*"))
    } # localblocks
  cat(paste(session,"done",dim(imat)[1],"by",dim(imat)[2],"assembled so far\n"))
} # session 
imat<-data.frame(imat)
antsImageWrite( as.antsImage( data.matrix( imat ) ) , afn )
if ( sum( usedesignrow ) != nrow(imat) ) stop("  sum( usedesignrow ) != nrow(imat) ")
dmat<-dmat[usedesignrow,]
if ( nrow(dmat) != nrow(imat) ) print("CHECK DIMENSIONS MATCH!")
write.csv(dmat,dfn,row.names=F)
} # existence
imat<-data.frame(imat)
return( list( dmat=dmat, imat=imat, usedesignrow=usedesignrow, subaal=subaal ) )
}
