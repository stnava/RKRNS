\name{spatioTemporalProjectionImage}
\alias{spatioTemporalProjectionImage}

\title{
Projects a 4D image to a 3D and 1D image (ND to 1 and D minus 1)
}

\description{
Applies a summary function to a spatiotemporal matrix to get lower dimensional representations.
}

\usage{
wmat<-spatioTemporalProjectionImage( spacetimevec, timelength, summaryFUN=mean, nDminus1ref=NA ) 
}

\arguments{
\item{spacetimevec}{
  input space time vector of length equal to nvoxels in a mask times the number of time points
}
\item{timelength}{
  the number of time points
}
\item{summaryFUN}{
  how to summarize the spatial part across time
}
\item{nDminus1ref}{
  reference spatial image
}
}

\value{
list(spacefunction=myestimatedbrainregionsval, timefunction=myestimatedhrf, spaceimage=spaceimage )
}

\author{
Avants BB
}

\examples{
\dontrun{
out<-spatioTemporalProjectionImage( vec, ) 
}
}
