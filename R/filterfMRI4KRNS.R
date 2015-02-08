#' Filtering options specific to KRNS
#' 
#' Filters the input matrix using butterworth filtering and time-series trends.
#' 
#' 
#' @param mat input bold time series matrix
#' @param filterlowfrequency get rid of variation that takes many seconds
#' @param filterhighfrequency get rid of variation that takes few seconds
#' @param trendfrequency time series trend option
#' @param winsorval remove spikes in signal given this winsorization value
#' @param removeSentLengthEffects nuisance vector of size nrow of matrix
#' @param removeEventOverlap nuisance vector of size nrow of matrix
#' @return matrix is output
#' @author Avants BB
#' @examples
#' 
#' fmat<-filterfMRI4KRNS( mat ) 
#' 
filterfMRI4KRNS <- function( imat,
    tr=as.numeric(0.5),
    filterlowfrequency=0.1, # 0.05 if you multiply by TR
    filterhighfrequency=1.0, # 0.4 # because of expected bold response < 25secs, > 5 seconds
    trendfrequency=3,
    trendfrequency2=NA,
    winsorval=0.01,
    removeSentLengthEffects=NA,
    removeEventOverlap=NA 
    )
{
################## filter the fmri ##################
  imatf<-imat
  if ( !is.na(removeSentLengthEffects) &  is.na(removeEventOverlap) ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+removeSentLengthEffects   ) ) )
  if (  is.na(removeSentLengthEffects) & !is.na(removeEventOverlap) ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+removeEventOverlap   ) ) )
  if ( !is.na(removeSentLengthEffects) & !is.na(removeEventOverlap) ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+removeSentLengthEffects +removeEventOverlap  ) ) )
  colnames(imatf)<-colnames(imat)
  fh<-filterhighfrequency
  fl<-filterlowfrequency
#  imatb<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" )) # , opt="trig" ))
  colnames(imatf)<-colnames(imat)
  for ( i in 1:ncol(imatf) ) {
      if (!is.na(trendfrequency)) imatf[,i]<-data.frame(stl( ts( winsor(imatf[,i],winsorval),frequency=trendfrequency)  , "per" )$time.series)$trend
      if (!is.na(trendfrequency2)) {
          temp<-data.frame(stl( ts(        imatf[,i]          ,frequency=trendfrequency2) , "per" )$time.series)$trend
          imatf[,i]<-imatf[,i]-temp
      }
  }
  return( imatf ) 
} 
