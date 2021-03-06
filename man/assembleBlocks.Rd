% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembleBlocks.R
\name{assembleBlocks}
\alias{assembleBlocks}
\title{assemble image blocks of bold}
\usage{
assembleBlocks(bmask, aalimg, labs, datadir, imagepostfix,
  assembledDesignOutPrefix, assembledImageOutPrefix, dmat, usedesignrow,
  imat = NA, ncompcor = 6, zscore = TRUE, spatialsmooth = 0)
}
\arguments{
\item{bmask}{brain mask - compcor is done within this mask}

\item{aalimg}{a label image - does not have to be aal}

\item{labs}{which labels of the above label image to use e.g. c(1,2,3)}

\item{datadir}{the input data directory}

\item{imagepostfix}{the image blocks postfix e.g. motioncorr.nii.gz}

\item{assembledImageOutPrefix}{output file name for assembled bold file -
will be read if it already exists (rather than reassembled)}

\item{usedesignrow}{organizational boolean tracking which parts of the
design file have concurrent BOLD images, useful for sanity checks}

\item{imat}{the current estimate of the assembled image matrix}

\item{zscore}{should we zscore each voxel?}

\item{assembledBlocksOutPrefix}{output file name for assembled design file -
will be read if it already exists (rather than reassembled)}

\item{spatialsmoothing}{scalar value greater than or equal to zero that
controls voxel-wise smoothing in 4D}
}
\value{
list is output containing organizational variables for design matrix
}
\description{
Assembles the bold study from raw block data
}
\examples{

\dontrun{
assembly2<-assembleBlocks( bmask, aalimg, labs, datadir, imagepostfix,  dfn, afn, dmat, usedesignrow, imat )
}

}
\author{
Avants BB, Phillips JS
}
