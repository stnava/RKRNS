\name{simulateBOLD}
\alias{simulateBOLD}
\title{simulated bold study with controllable design and noise parameters}
\usage{
  bb<-simulateBOLD(option="basic")
  bb<-simulateBOLD(option="henson",eximg=ff)
}
\description{
  the henson example produces a mask, an image and a design matrix.
  the basic example produces a noisy data matrix and design matrix. 
  ... see args for some of the options and reasonable defaults.
}
\example{
# get example image
fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
eximg<-antsImageRead(fn,3)
bb<-simulateBOLD(option="henson",eximg=eximg)
# now test glmDenoise
mat<-timeseries2matrix( bb$simbold, bb$mask )
dd<-glmDenoiseR( mat, bb$desmat, crossvalidationgroups=8 , maxnoisepreds=1:8,
  hrfbasislength=50, selectionthresh=0.05 ,collapsedesign=T, polydegree=4 )
}
