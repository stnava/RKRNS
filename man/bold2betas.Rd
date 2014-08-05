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
  betas<-bold2betas( imat, dmats, dmat$blockNumb, maxnoisepreds=c(12,14,18) )
# example of what you might do for decoding
  ee2<-bold2betas( data.matrix(imat[samp1,]) , dmats[samp1,], dmat$blockNumb[samp1], bl=30, polydegree=8, whichcols=whichcols ) 
  eanat<-sparseDecom( data.matrix(ee2$eventbetas), nvecs=10, sparseness=0.01, mycoption=1, inmask=subaal, cthresh=50 )
  rownames(eanat$projections)<-rownames(ee2$eventbetas)
  pheatmap(cor(t(eanat$projections)))
  classes<-rep("",nrow(eanat$projections))
  for ( i in 1:nrow(eanat$projections)) classes[i]<-strsplit(rownames(eanat$projections)[i],"[.][.]")[[1]][1]
  classes<-as.factor(classes)
  inds<-sample(1:nrow(ee2$eventbetas))[1:8]
  mydf<-data.frame( lab=classes, p=(data.matrix(ee2$eventbetas)) %*% t(data.matrix(imageListToMatrix(eanat$eigenanatomyimages,subaal)) ))
  mdl<-randomForest(lab~., mydf[inds,] )
  print( predict( mdl, newdata=mydf[-inds,]) )
}
