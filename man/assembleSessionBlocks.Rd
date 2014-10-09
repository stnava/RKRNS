\name{assembleSessionBlocks}
\alias{assembleSessionBlocks}

\title{
compute betas from image blocks of bold
}

\description{
Computes event-wise betas for each image block and writes to an image
}

\usage{
assembly2<-assembleSessionBlocks( bmask, aalimg, labs, datadir, imagepostfix,
   assembledBlocksOutPrefix, afn, assembledImageOutPrefix,
   usedesignrow, hrf, zscore=TRUE, spatialsmoothing=0 )
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
\item{spatialsmoothing}{
  scalar value greater than or equal to zero that controls voxel-wise smoothing in 4D
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
hrf<-ts( hemodynamicRF( 30, onsets=2, durations=1, rt=0.5,cc=0.1,a1=6,a2=12,b1=0.6, b2=0.6 ) )
assembly2<-assembleSessionBlocks( bmask, aalimg, labs, datadir, imagepostfix,  dfn, afn, dmat, usedesignrow, imat )
}
}
