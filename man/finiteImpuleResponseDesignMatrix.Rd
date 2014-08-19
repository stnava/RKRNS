\name{finiteImpuleResponseDesignMatrix}
\alias{finiteImpuleResponseDesignMatrix}

\title{
Simple finiteImpuleResponseDesignMatrix function.
}

\description{
Convert a design matrix (or any matrix) to a time-shifted version of itself 
as needed by finite impulse response models used when deconvolving hemodynamic
response functions.
}

\usage{
imat<-finiteImpuleResponseDesignMatrix( mat , 5 )
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
wmat<-finiteImpuleResponseDesignMatrix( mat , 5 ) 
}
