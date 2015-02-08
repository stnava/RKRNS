#' samples an array of labels nperm times to estimate the random classification
#' rate
#' 
#' simple empirical estimate of random classification rate
#' 
#' 
#' @examples
#' 
#' mylabs<-sample( rep( c("A","B","C","D"), 100 ) )
#' randomClassificationRate( mylabs , 100 )
#' 
randomClassificationRate <- function( myclasses, nperm=100 ) {
    nlab<-length( myclasses )
    randcorr<-0
    for ( i in 1:nperm )
        randcorr<-( 100.0 * sum(as.numeric( myclasses == sample(myclasses) ) )/nlab )
    return( randcorr )
}


