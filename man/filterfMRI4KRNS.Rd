\name{filterfMRI4KRNS}
\alias{filterfMRI4KRNS}

\title{
Filtering options specific to KRNS
}

\description{
Filters the input matrix using butterworth filtering and time-series trends.
}

\usage{
fmat<-filterfMRI4KRNS( inputMatrix,
    tr=as.numeric(0.5),
    filterlowfrequency=0.1, 
    filterhighfrequency=1.0, 
    trendfrequency=3,
    winsorval=0.01,
    removeSentLengthEffects=NA,
    removeEventOverlap=NA  )
}

\arguments{
\item{mat}{
  input bold time series matrix
}
\item{filterlowfrequency}{
  get rid of variation that takes many seconds
}
\item{filterhighfrequency}{
  get rid of variation that takes few seconds
}
\item{trendfrequency}{
  time series trend option 
}
\item{winsorval}{
  remove spikes in signal given this winsorization value
}
\item{removeSentLengthEffects}{
  nuisance vector of size nrow of matrix 
}
\item{removeEventOverlap}{
  nuisance vector of size nrow of matrix  
}
}

\value{
matrix is output
}

\author{
Avants BB
}

\examples{
fmat<-filterfMRI4KRNS( mat ) 
}
