spatioTemporalProjectionImage <- function( spacetimematrixIn, summaryFUN=mean, nDminus1ref=NA ) {
  if ( class(spacetimematrixIn)[[1]] == "antsImage" )
    {
    mm<-timeseries2matrix( spacetimematrixIn, nDminus1ref )
    } else mm<-spacetimematrixIn
  myestimatedhrf<-apply(mm,FUN=summaryFUN,MARGIN=1)
  mmmag<-sqrt( mm^2 )
  myestimatedbrainregionsval<-apply(mmmag,FUN=summaryFUN,MARGIN=2)
  spaceimage<-NA
  if ( !is.na(nDminus1ref) )
    {
    spaceimage<-antsImageClone( nDminus1ref )
    indexvec<-spaceimage > 0
    if ( sum(indexvec) == length( myestimatedbrainregionsval) )
      {
      spaceimage[ indexvec ]<-0
      spaceimage[ nDminus1ref > 0 ]<-myestimatedbrainregionsval
      }
  }
  return( list(spacefunction=myestimatedbrainregionsval, timefunction=myestimatedhrf, spaceimage=spaceimage ) )
} 
