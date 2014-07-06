\name{assembleDesign}
\alias{assembleDesign}

\title{
Simple assembleDesigning function.
}

\description{
Assembles the study design from raw block data
}

\usage{
assembleDesign( blocksCSVlist , datadir, assembledDesignOutPrefix,  assembledImageOutPrefix )
}

\arguments{
\item{blocksCSVlist}{
  a list of block-wise csv files detailing the design
}
\item{datadir}{
  root data directory containing subject design, bold, etc
}
\item{assembledDesignOutPrefix}{
  output file name for assembled design file - will be read if it already exists (rather than reassembled)
}
\item{assembledImageOutPrefix}{
  output file name for assembled bold file - will be read if it already exists (rather than reassembled)
}
}

\value{
list is output containing organizational variables for design matrix
}

\author{
Avants BB, Phillips JS
}

\examples{
\dontrun{
assembly<-assembleDesign( blocksCSVlist, datadir, dfn, afn )
}
}
