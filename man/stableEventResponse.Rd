\name{stableEventResponse}
\alias{stableEventResponse}
\title{Estimate effect sizes per-voxel, per stimulus, cross-validated}
\description{
Use leave-a-run-out to summarize response variability across runs for
a stimulus class.  Returns a matrix of beta effect sizes per stimulus 
class where the beta effect size are defined as mean of the beta response
across left out runs divided by standard deviation of the same.
}
\usage{
eventbetas<-stableEventResponse(  boldmatrix, designmatrixIn, runIDs,
                                polydegree=4, baseshift=0, hrf=NA, timevals=NA )
}

\arguments{
\item{boldmatrix}{
  input raw bold data in time by space matrix 
}
\item{designmatrix}{
  input design matrix - binary/impulse entries for event related design, blocks otherwise
}
\item{runIDs}{
  numbers for the rows that should be treated together as runs 
}
\item{polydegree}{
  number of polynomial predictors
}
\item{hrf}{
  input hrf to use in regression
}
\item{baseshift}{
  basis shift for design matrix
}
}

\value{
returns a list of cross-validated effect sizes for different stimuli classes
}

\author{
Avants BB
}

\examples{
fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
eximg<-antsImageRead(fn,3)
fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="") 
mask<-antsImageRead(fn,3)
bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
boldImage<-bb$simbold
mat<-timeseries2matrix( bb$simbold, bb$mask )
runs<-bb$desmat$Run; 
hrf<-hemodynamicRF( 20, onsets=1, durations=1, rt=1,cc=0.1 )
stb<-stableEventResponse(mat,  bb$desmat[,1:4], runs,  hrf=hrf )
var1i<-antsImageClone( mask ) 
var1i[mask==1]<-stb[1,] # then write this out ...
# or run eigseg 
stb2<-stb
stb2[ stb < 6 ]<-0
ee<-eigSeg(mask, matrixToImages( stb2, mask) )
ImageMath(3,ee,'ClusterThresholdVariate',ee,mask,5)
# antsImageWrite ...
}
