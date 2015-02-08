#' Simple finiteImpulseResponseDesignMatrix function.
#' 
#' Convert a design matrix (or any matrix) to a time-shifted version of itself
#' as needed by finite impulse response models used when deconvolving
#' hemodynamic response functions.
#' 
#' 
#' @param mat input matrix
#' @param n shift amount
#' @param baseshift baseline shift amount so total shift is n+baseshift
#' @return matrix is output
#' @author Avants BB
#' @examples
#' 
#' mat<-diag(10) 
#' wmat<-finiteImpulseResponseDesignMatrix( mat , 5 ) 
#' 
finiteImpulseResponseDesignMatrix <- function(x,n=1,baseshift=0) {
  if (nargs() == 0) {
    print( args( finiteImpulseResponseDesignMatrix ) )
    return(1)
  }
  nc<-ncol(x)
  nr<-nrow(x)
  xi<-matrix( rep(0,nr*nc*n), nrow=nr )
  colnames(xi)<-paste("V",1:ncol(xi))
  colnamesxi<-colnames(xi)
  colnamesx<-colnames(x)
  j<-1
  for ( i in 1:ncol(x) )
    {
    basecol<-shift(x[,i],baseshift)
    shiftmat<-basecol
    for ( k in 2:n ) shiftmat<-cbind(shiftmat,shift(basecol,k-1))
    xi[,j:(j+n-1)]<-shiftmat
    colnamesxi[j:(j+n-1)]<-paste(colnamesx[i],1:n,sep='')
    j<-j+n
    }
  colnames(xi)<-colnamesxi
  return(xi)
} 
