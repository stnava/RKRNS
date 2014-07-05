\name{classificationNetwork}
\alias{classificationNetwork}

\title{
Visualize a (mis)classification network.
}

\description{
Input a prediction and an outcome vector that are in correspondence, get out a graph.
}

\usage{
wmat<-classificationNetwork( mat  )
}

\arguments{
\item{mat}{
  input matrix
}
\item{k}{
  rank to use 
}
\item{reducex}{
  reduce the input matrix to k-size subspace
}
}

\value{
matrix is output
}

\author{
Avants BB
}

\examples{
\dontrun{
  sentencesubset <- sentencedf$sentences %in% unique(fspacenames[redlist[l2]])
  nodedf <- data.frame( nodename=sentencedf[sentencesubset,1], nodeid=sentencedf[sentencesubset,2] )
  ww <- misclassnetwork( nodesIn=nodedf, wclassesf[l2], pred, outfile='temp2.html', mycharge=-2066, zoom=T)
}
}
