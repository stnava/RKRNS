randomClassificationRate <- function( myclasses, nperm=100 ) {
    nlab<-length( myclasses )
    randcorr<-0
    for ( i in 1:nperm )
        randcorr<-( 100.0 * sum(as.numeric( myclasses == sample(myclasses) ) )/nlab )
    return( randcorr )
}


