print("#########assemble image blocks, potentially event-specifically#########")
maskdim<-dim( bmask )
aalmask<-antsImageClone( aalimg )
aalmask[ aalmask > 1 ]<-1
subaal<-antsImageClone( aalimg )
subaal[ aalimg > 0 ]<-0
for ( lab in labs ) subaal[ aalimg == lab ]<-lab
print(paste("assemble blocks for",afn))
if ( ! exists("imat") ) { # assembly begin ......
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
    locmat<-timeseries2matrix( img , subaal )
    for( j in 1:throwaway ) locmat[j,]<-apply( locmat[(throwaway+1):nrow(locmat),], FUN=mean, MARGIN=2 )
    if ( ncompcor > 0 ) {
        locmatfull<-timeseries2matrix( img , bmask )
        if ( ! exists("mycompcorv") ) mycompcorv<-compcor( locmatfull, ncompcor=ncompcor, variance_extreme = 0.97, returnv=TRUE )
        highvarmat<-compcor( locmatfull, ncompcor, variance_extreme = 0.97, returnhighvarmat=TRUE )
        mycompcor<- scale(highvarmat %*% t(mycompcorv) )
#        mycompcor<-compcor( locmatfull, ncompcor, variance_extreme = 0.97 )
        locmat<-residuals(lm(locmat~0+mycompcor))
    }
    locmat<-locmat/mean(locmat)
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
colnames(imat)<-aal$label_name[labs]
write.csv(imat,afn,row.names=F)
if ( sum( usedesignrow ) != nrow(imat) ) stop("  sum( usedesignrow ) != nrow(imat) ")
dmat<-dmat[usedesignrow,]
if ( nrow(dmat) != nrow(imat) ) print("CHECK DIMENSIONS MATCH!")
write.csv(dmat,dfn,row.names=F)
} # existence
imat<-data.frame(imat)
colnames(imat)<-aal$label_name[labs]
print("assembly done")
