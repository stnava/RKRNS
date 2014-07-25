\name{glmDenoiseR}
\alias{glmDenoiseR}
\title{Optimize a regression model for BOLD based on cross-validated model fitting}
\usage{
  glmDenoiseR( boldmatrix, designmatrix )
}
\description{
Inspired by discussion with Kendrick Kay regarding his glm denoise tool.
}

\usage{
regressionmodeldataframe<-glmDenoiseR( boldmatrix , designmatrix )
}

\arguments{
\item{boldmatrix}{
  input raw bold data in time by space matrix 
}
\item{boldmatrix}{
  input design matrix - binary/impulse entries for event related design, blocks otherwise
}
}

\value{
returns a data frame
}

\author{
Avants BB
}

\examples{
mat<-replicate(25, rnorm(100)) 
dmat<-mat[,1:2]*0
dmat[,1]<-as.numeric( mat[ ,1] > 0.99)
dmat[,2]<-as.numeric( mat[ ,2] > 0.99)
regressionmodeldataframe<-glmDenoiseR( mat , dmat )
}
