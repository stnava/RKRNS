timeserieswindow2freqamp <- function( timeseriesmatrix , eventlist, timewindow, f, wl )
{
  require( seewave )
  if ( length(dim(timeseriesmatrix)) != 2 )
    {
    print("Mask should be of dimensionality 3")
    return(NA)
    }
  nevents<-length(eventlist)
  if ( min(eventlist) < 0 ) {
    print("ERROR: min(eventlist) < 0 ")
    return(NA)
  }
  if ( max(eventlist) > nrow(timeseriesmatrix) ) {
    print("ERROR: max(eventlist) > nrow(timeseriesmatrix) ")
    return(NA)
  }
  nvox<-ncol( timeseriesmatrix )
  outmatf<-matrix( rep(0, nvox*nevents) , nrow = nevents )
  outmata<-outmatf
  outmath<-outmatf
  locmatTransform <- function( xin , lf=f, lwl=wl )
    {
    x<-oscillo(xin,f, plot=F)
    xentropy<-H(x,f=lf,wl=lwl)
#    return( c(0,0,xentropy) )
    spec <- meanspec( x, f=lf, wl=lwl, plot=F)
    fanda<-fpeaks(spec,plot=F)
    if ( anyNA(fanda) ) return( c(0,0,xentropy) )
    return(  c( fanda[1,], xentropy ) )
    }
  locSpec <- function( xin )
    {
    s1<-spectrum( xin , plot=FALSE )
    return( c(s1$freq[ which.max(s1$spec) ], max(s1$spec) ) )
    }
  ct<-1
  progress <- txtProgressBar(min = 0, max = length(eventlist), style = 3)
  for ( i in eventlist ) {
    maxrow<-i+(timewindow-1)
    if ( maxrow > nrow(timeseriesmatrix) ) maxrow<-nrow(timeseriesmatrix)
    locmat<-( timeseriesmatrix[ i:maxrow, ] )
#    fff<-apply( locmat, FUN=locmatTransform, MARGIN=2 )
    fff<-apply( locmat, FUN=locSpec, MARGIN=2 )
    outmatf[i,]<-fff[1,]
    outmata[i,]<-fff[2,]
#    outmath[i,]<-fff[3,]
    ct<-ct+1
    setTxtProgressBar( progress, ct )
  }
  close(progress)
  return( list(frequencyTransform=outmatf, amplitudeTransform=outmata, entropyTransform=outmath ) )
} 
