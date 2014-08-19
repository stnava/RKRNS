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
# get example image
fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
eximg<-antsImageRead(fn,3)
fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="") 
mask<-antsImageRead(fn,3)
bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
mat<-timeseries2matrix( bb$simbold, bb$mask )
runs<-bb$desmat$Run; 
### cca betas
mysp<-c( -0.01, -0.9 )
btsc<-bold2betasCCA( data.matrix(mat)  , bb$desmat[,1:4], blockNumb=runs,
     bl=25, baseshift=0, sparseness=mysp, bestvoxnum=10, mask=mask, 
     polydegree=1, mycoption=1, its=12, onlyhrf=FALSE )
plot(ts(t(btsc$runhrfs)))
# from here see the help for bold2betas
}
