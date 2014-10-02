\name{finiteImpulseResponseDesignMatrix}
\alias{finiteImpulseResponseDesignMatrix}

\title{
Simple finiteImpulseResponseDesignMatrix function.
}

\description{
Convert a design matrix (or any matrix) to a time-shifted version of itself 
as needed by finite impulse response models used when deconvolving hemodynamic
response functions.
}

\usage{
imat<-finiteImpulseResponseDesignMatrix( mat , 5 )
}

\arguments{
\item{mat}{
  input matrix
}
\item{n}{
  shift amount
}
\item{baseshift}{
  baseline shift amount so total shift is n+baseshift
}
}

\value{
matrix is output
}

\author{
Avants BB
}

\examples{
mat<-diag(10) 
wmat<-finiteImpulseResponseDesignMatrix( mat , 5 ) 
}
