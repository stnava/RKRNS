\name{bold2betasCCA}
\alias{bold2betasCCA}
\title{Convert a bold time series and design matrix to event-wise cca-based betas}
\description{
  Uses cca and nuisance variables to estimate multivariate betas per event.
}
\usage{
  eventbetas<-bold2betasCCA(  boldmatrix, designmatrix, blockNumb, bl=12, baseshift=5, mask=NA, sparseness=c(0,0), multievents=FALSE, polydegree=10, bestvoxnum=50, uselm=FALSE , nvecs=5  )
}

\arguments{
\item{boldmatrix}{
  input raw bold data in time by space matrix 
}
\item{designmatrix}{
  input design matrix - binary/impulse entries for event related design, blocks otherwise
}
\item{blockNumb}{
  numbers for the rows that should be treated together as runs 
}
\item{nvecs}{
  number of cca predictors to explore e.g. 5
}
\item{polydegree}{
  number of polynomial predictors
}
\item{bl}{
  basis length for hrf estimation
}
\item{baseshift}{
  number of time points to ignore post-event onset
}
}

\value{
returns a list with relevant output 
}

\author{
Avants BB
}

\examples{
 subaal<-NA
 mysp<-c(-0.01,0.9) 
 ee1<-bold2betasCCA( i1 , dmats[samp1,], dmat$blockNumb[samp1], bl=12, baseshift=5, sparseness=mysp, bestvoxnum=50, mask=subaal ) 
}
