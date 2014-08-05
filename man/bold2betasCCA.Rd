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
 subaal<-antsImageRead("subaal.nii.gz",3)
 mysp<-c(-0.01,0.9) 
 samp1<-1:1000
 ee1<-bold2betasCCA( data.matrix(imat[samp1,]) , dmats[samp1,], dmat$blockNumb[samp1], bl=12, baseshift=5, sparseness=mysp, bestvoxnum=50, mask=subaal ) 

mysp<-c(-0.02,-0.9) 
whichcols<-grep("coffe",colnames(dmats))
samp1<-1000:11000
ee1<-bold2betasCCA( data.matrix(imat[samp1,]) , dmats[samp1,], dmat$blockNumb[samp1], bl=12, baseshift=4, sparseness=mysp, bestvoxnum=10, mask=subaal, polydegree=4, whichcols=whichcols ) 

# example of what you might do for decoding
  ee2<-bold2betasCCA( data.matrix(imat[samp1,]) , dmats[samp1,], dmat$blockNumb[samp1], bl=12, baseshift=4, sparseness=mysp, bestvoxnum=10, mask=subaal, polydegree=8, whichcols=whichcols ) 
  eanat<-sparseDecom( data.matrix(ee2$eventbetas), nvecs=10, sparseness=-0.01, mycoption=1 )
  rownames(eanat$projections)<-rownames(ee2$eventbetas)
  pheatmap(cor(t(eanat$projections)))
  classes<-rep("",nrow(eanat$projections))
  for ( i in 1:nrow(eanat$projections)) classes[i]<-strsplit(rownames(eanat$projections)[i],"[.][.]")[[1]][1]
  classes<-as.factor(classes)
  inds<-sample(1:nrow(ee2$eventbetas))[1:8]
  mydf<-data.frame( lab=classes, p=(data.matrix(ee2$eventbetas)) %*% data.matrix(eanat$eigenanatomyimages) )
  mdl<-svm(lab~., mydf[inds,] )
  print( predict( mdl, newdata=mydf[-inds,]) )

}
