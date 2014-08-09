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
runs<-bb$desmat$Run; runs[]<-1
### cca betas
btsc<-bold2betasCCA( data.matrix(mat)  , bb$desmat[,1:4], blockNumb=runs,
     bl=15, baseshift=0, sparseness=mysp, bestvoxnum=50, mask=mask, 
     polydegree=2, mycoption=0, its=12, uselm=1 )
zz<-apply(btsc$eventbetas,FUN=mean,MARGIN=2)
zze<-zz/apply(btsc$eventbetas,FUN=sd,MARGIN=2)
inds<-sample(1:nrow(btsc$eventbetas),size=round(nrow(btsc$eventbetas)*1./2.))
ff<-which( abs(zze) > 0.35 )
mylabs<-rep("",nrow(btsc$eventbetas))
for ( i in 1:nrow(btsc$eventbetas) ) mylabs[i]<-substr( rownames(btsc$eventbetas)[i],1,2)
mylabs<-as.factor(mylabs)
mydf<-data.frame( lab=mylabs,  vox=data.matrix(btsc$eventbetas)[,ff] )
mdl<-svm( lab ~., data=mydf[inds,])
err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
print(paste("NPredVox",length(ff), "Correct",err*100))
}
