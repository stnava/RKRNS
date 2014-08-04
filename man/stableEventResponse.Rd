\name{stableEventResponse}
\alias{stableEventResponse}
\title{Denoise and estimate response variability across runs}
\description{
Use extended glmdenoise and cross-validation to summarize 
response variability across runs for an event.  Returns a
matrix of betas where there is a beta vector for each left 
out run.  The betas are t-tstats for the selected event.
}
\usage{
eventbetas<-stableEventResponse( boldmatrix, designmatrix, 
    blockNumb, eventselect='coffee', polydegree=10, bl=50,
    crossvalidationgroups=4, selectionthresh=0.1, maxnoisepreds=4 )
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
\item{eventselect}{
  selection criteria determining which events are estimated
}
\item{maxnoisepreds}{
  number of noise predictors to explore e.g. a list of values 1 to 10 or a few values to try
}
\item{polydegree}{
  number of polynomial predictors
}
\item{bl}{
  basis length for hrf estimation
}
\item{crossvalidationgroups}{
  number of data splits in glm cross-validation
}
}

\value{
returns a list of cross-validated betas 
}

\author{
Avants BB
}

\examples{
ff<-stableEventResponse(imat,dmats, dmat$blockNumb, 
  eventselect='coffee.was.hot',maxnoisepreds=14)
ff<-stableEventResponse(imat,dmats, dmat$blockNumb, 
  eventselect=grep('coffee.was.hot',colnames(dmats)),
  maxnoisepreds=14)
}
