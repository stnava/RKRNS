\name{assembleBlocks}
\alias{assembleBlocks}

\title{
assemble image blocks of bold
}

\description{
Assembles the bold study from raw block data
}

\usage{
assembly2<-assembleBlocks( bmask, aalimg, labs, datadir, imagepostfix,  assembledBlocksOutPrefix, afn, assembledImageOutPrefix, usedesignrow, imat, zscore=TRUE )
}

\arguments{
\item{bmask}{
  brain mask - compcor is done within this mask
}
\item{aalimg}{
  a label image - does not have to be aal
}
\item{labs}{
  which labels of the above label image to use e.g. c(1,2,3)
}
\item{datadir}{
  the input data directory
}
\item{imagepostfix}{
  the image blocks postfix e.g. motioncorr.nii.gz
}
\item{assembledBlocksOutPrefix}{
  output file name for assembled design file - will be read if it already exists (rather than reassembled)
}
\item{assembledImageOutPrefix}{
  output file name for assembled bold file - will be read if it already exists (rather than reassembled)
}
\item{usedesignrow}{
  organizational boolean tracking which parts of the design file have concurrent BOLD images, useful for sanity checks
}
\item{imat}{
  the current estimate of the assembled image matrix 
}
\item{zscore}{
  should we zscore each voxel?
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
assembly2<-assembleBlocks( bmask, aalimg, labs, datadir, imagepostfix,  dfn, afn, dmat, usedesignrow, imat )
}
}
