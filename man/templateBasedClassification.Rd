\name{templateBasedClassification}
\alias{templateBasedClassification}

\title{
Use an expected function shape to classify new instances
}

\description{
Applies a function to a matrix representation of a sentence to get an eigensentence map.
}

\usage{
newlabels<-templateBasedClassification( exemplarmat, labels, newmat, method="corr"  )
}

\arguments{
\item{exemplarmat}{
  each row is an example of a specific class - no specific ordering assumed
}
\item{labels}{
  labels for each row of the exemplar
}
\item{newmat}{
  newmatrix instances to classify
}
\item{method}{
  one of a selection of options to drive the classification
}
}

\value{
class label is output
}

\author{
Avants BB
}

\examples{
  mat<-replicate(100, rnorm(21)) 
  labels<-rep(c(1:3),7)
  newmat<-replicate(100, rnorm(21))[1:2,]
  myclass<-templateBasedClassification( mat, labels, newmat, method="corr" )
}
