\name{glmDenoiseR}
\alias{glmDenoiseR}
\title{Optimize a regression model for BOLD based on cross-validated model fitting}
\description{
Inspired by discussion with Kendrick Kay regarding his glm denoise tool 
http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00247/abstract
 0. estimate hrf using assumed function or finite impulse response (FIR)
 1. regressors include: design + trends + noise-pool
 2. find noise-pool by initial cross-validation without noise regressors
 3. cross-validate predictions using different numbers of noise regressors
 4. select best n for predictors from noise pool
 5. return the noise mask and the value for n
}
\usage{
  regressionmodeldataframe<-glmDenoiseR( boldmatrix, designmatrix , hrfBasis=NA, crossvalidationgroups=4, hrfShifts = 4, selectionthresh=0.1, maxnoisepreds=6 )
}

\arguments{
\item{boldmatrix}{
  input raw bold data in time by space matrix 
}
\item{designmatrix}{
  input design matrix - binary/impulse entries for event related design, blocks otherwise
}
\item{hrfBasis}{
  basis function for assumed HRF otherwise use FIR
}
\item{hrfShifts}{
  n-shifts of assumed hrf - shifts by 1 or, for FIR, length of estimated HRF
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
boldImage<-bb$simbold
mat<-timeseries2matrix( bb$simbold, bb$mask )
runs<-bb$desmat$Run; 
# finite impulse response
hrfbasislength<-20
dd<-glmDenoiseR( mat, bb$desmat[,1:4], hrfBasis=NA, hrfShifts = hrfbasislength, 
  crossvalidationgroups=runs, maxnoisepreds=c(0,1,4,6,10,14) , selectionthresh=0.1 , 
  collapsedesign=F, polydegree=4, reestimatenoisepool = F )
# average of assumed HRFs
tr<-1
a1<-4
a2<-10
hrf<-hemodynamicRF( hrfbasislength, onsets=2, 
  durations=tr, rt=tr,cc=0.1,a1=a1,a2=a2,b1=0.9, b2=0.9 )
plot(ts(hrf))
dd2<-glmDenoiseR( mat, bb$desmat[,1:4], hrfBasis=hrf, hrfShifts = 4 , crossvalidationgroups=runs, 
  maxnoisepreds=4 , selectionthresh=0.1 , collapsedesign=T, polydegree=4 )
# or refine FIR 
dd3<-glmDenoiseR( mat, bb$desmat[,1:4], hrfBasis=shift(dd$hrf,-2), hrfShifts = 4 , crossvalidationgroups=runs, 
  maxnoisepreds=0:2 , selectionthresh=0.1 , collapsedesign=T, polydegree=4 )
}
