\name{eigenSentences}
\alias{eigenSentences}

\title{
Simple eigenSentence estimation function.
}

\description{
Applies a function to a matrix representation of a sentence to get an eigensentence map.
}

\usage{
esent<-eigenSentence(  wordembed, functiontoapply = mean , normalize=T, eigsentbasislength = 5,  sentencesIn=NA  )
}

\arguments{
\item{wordembed}{
  the words and their vector embedding - a data frame with n basis length columns 
}
\item{functiontoapply}{
  e.g. mean, max, min, a custom function --- needs to work with apply
}
\item{normalize}{
  normalize the magnitude of the eigsentences to be 1
}
\item{eigsentbasislength}{
  less than or equal to the full basis length (100)
}
}

\value{
eigensentence matrix is output
}

\author{
Avants BB
}

\examples{
  data(reuters_words,package="RKRNS")
  esent<-eigenSentences( reuters_words )
}
