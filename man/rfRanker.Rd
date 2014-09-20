\name{rfRanker}
\alias{rfRanker}
\title{use random forest to predict and rank sentences from features}
\usage{
  rfRanker( mydf, trainingRows, usesparse=F )
}
\description{
  takes a dataframe that contains a column named lab such that mydf$lab gives the 
  true labeling of the data. trainingRows is a vector identifying which rows
  should be training data.  success is defined by a prediction probability being 
  within the top targetSuccessRate*100 percent.  Returns 
}
\examples{
mylabs<-rep(c("a","b","c"),5)
voxeldata<-replicate(100, rnorm(length(mylabs)))
mydf<-data.frame( lab=mylabs, vox=voxeldata)
rfr<-rfRanker( mydf, 1:round(nrow(mydf)/2) )
print( rfr$successPercent )
}
