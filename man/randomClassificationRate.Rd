\name{randomClassificationRate}
\alias{randomClassificationRate}
\title{samples an array of labels nperm times to estimate the random classification rate}
\usage{
  randomClassificationRate( labs, 100 )
}
\description{
  simple empirical estimate of random classification rate
}
\examples{
mylabs<-sample( rep( c("A","B","C","D"), 100 ) )
randomClassificationRate( mylabs , 100 )
}
