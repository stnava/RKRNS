#' Projects a 4D image to a 3D and 1D image (ND to 1 and D minus 1)
#' 
#' Applies a summary function to a spatiotemporal matrix to get lower
#' dimensional representations.
#' 
#' 
#' @param spacetimevec input space time matrix or antsImage of length equal to
#' nvoxels in a mask times the number of time points
#' @param summaryFUN how to summarize the spatial part across time
#' @param nDminus1ref reference spatial image
#' @return list(spacefunction=myestimatedbrainregionsval,
#' timefunction=myestimatedhrf, spaceimage=spaceimage )
#' @author Avants BB
#' @examples
#' 
#' \dontrun{
#' img<-antsImageRead("TEST0View1vec000.nii.gz",4)
#' ff<-spatioTemporalProjectionImage( img, sum, mask )
#' antsImageWrite(ff$spaceimage,'temp.nii.gz')
#' plot(ff$timefunction,type='l')
#' }
#' 
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
    spaceimage[ spaceimage < 0 ]<-0
    indexvec<-spaceimage > 0
    if ( sum(indexvec) == length( myestimatedbrainregionsval) )
      {
      spaceimage[ indexvec ]<-0
      spaceimage[ nDminus1ref > 0 ]<-myestimatedbrainregionsval
      }
  }
  return( list(spacefunction=myestimatedbrainregionsval, timefunction=myestimatedhrf, spaceimage=spaceimage ) )
} 
