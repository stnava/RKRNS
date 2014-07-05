filterfMRI4KRNS <- function( imat,
    tr=as.numeric(0.5),
    responselength=16, # e.g. 15 seconds div by 0.5 tr => 30 volumes
    filterlowfrequency=0.1, # 0.05 if you multiply by TR
    filterhighfrequency=1.0, # 0.4 # because of expected bold response < 25secs, > 5 seconds
    trendfrequency=3,
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
  imatb<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" )) # , opt="trig" ))
  colnames(imatf)<-colnames(imat)
  for ( i in 1:ncol(imatf) ) {
      imatf[,i]<-data.frame(stl( ts( winsor(imatb[,i],winsorval),frequency=trendfrequency) , "per" )$time.series)$trend
  }
  return( imatf ) 
} 
