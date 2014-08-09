\name{bold2betas}
\alias{bold2betas}
\title{Convert a bold time series and design matrix to event-wise betas}
\description{
Inspired by discussion with Kendrick Kay regarding his glm denoise tool 
http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00247/abstract.
Here, we apply glm denoise run by run to get event-wise betas.
}
\usage{
  eventbetas<-bold2betas( boldmatrix, designmatrix, blockNumb, maxnoisepreds, polydegree=10, bl=40, crossvalidationgroups=4  )
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
runs<-bb$desmat$Run; runs[]<-1
btsc<-bold2betas( boldmatrix=data.matrix(mat) , designmatrix=bb$desmat[,1:4], 
     blockNumb=runs, maxnoisepreds=0,
     bl=12, polydegree=2, selectionthresh=0.2 )
zz<-apply(btsc$eventbetas,FUN=mean,MARGIN=2)
zze<-zz/apply(btsc$eventbetas,FUN=sd,MARGIN=2)
inds<-sample(1:nrow(btsc$eventbetas),size=round(nrow(btsc$eventbetas)*1./2.))
ff<-which( abs(zze) > 0.5 )
mylabs<-rep("",nrow(btsc$eventbetas))
for ( i in 1:nrow(btsc$eventbetas) ) mylabs[i]<-substr( rownames(btsc$eventbetas)[i],1,2)
mylabs<-as.factor(mylabs)
mydf<-data.frame( lab=mylabs,  vox=data.matrix(btsc$eventbetas)[,ff] )
mdl<-svm( lab ~., data=mydf[inds,])
err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
print(paste("NPredVox",length(ff),"Correct",err*100))
}
