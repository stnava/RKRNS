\name{sccanBasedDecoder}
\alias{sccanBasedDecoder}

\title{
Simple eigenSentence estimation function.
}

\description{
Applies a function to a matrix representation of a sentence to get an eigensentence map.
}

\usage{
sccanBasedDecoder( eventdata, designmat, boldFeatureMatrix, sentenceSpace, mysparse=c(-0.1,-0.1), nvecs=5, its=1, smooth=0, cthresh=0, mask4d=NA, strategy=NA, doEanat=F, joinEanat=F, outgraphfile='temp2.html'  )
}

\arguments{
\item{eventdata}{
  output from annotateEvents
}
\item{designmat}{
  sentences part of the assembled design matrix
}
\item{boldFeatureMatrix}{
  a feature matrix for bold data, currently spatiotemporal
}
\item{sentenceSpace}{
  the sentence feature space
}
\item{doEanat}{
  also run a sparse PCA on the data
}
\item{other parameters}{
  see sparseDecom2 documentation
}
}

\value{
a data frame with annotation and ground truth vs prediction annotated ...
}

\author{
Avants BB
}

\examples{
\dontrun{
 ccaresults<-sccanBasedDecoder( eventdata, dmats, ccafeatspace , sentspace )
}
}
