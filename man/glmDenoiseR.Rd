\name{glmDenoiseR}
\alias{glmDenoiseR}
\title{Optimize a regression model for BOLD based on cross-validated model fitting}
\usage{
  glmDenoiseR( boldmatrix, designmatrix , hrfbasislength=50, kfolds=10, whichbase = NA, pvalthresh=0.25, maxnoisepreds=6 )
}
\description{
Inspired by discussion with Kendrick Kay regarding his glm denoise tool 
http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00247/abstract
 0. estimate hrf 
 1. regressors include: design + trends + noise-pool
 2. find noise-pool by initial cross-validation without noise regressors
 3. cross-validate predictions using different numbers of noise regressors
 4. select best n for predictors from noise pool
 5. return the noise mask and the value for n
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
returns a list with relevant output 
}

\author{
Avants BB
}

\examples{
mat<-replicate(250, rnorm(1000)) 
dmat<-mat[,1:2]*0
for ( i in 1:ncol(dmat) ) dmat[,i]<-as.numeric( mat[ ,i] > 0.9*max(mat[,i]))
regressionmodeldataframe<-glmDenoiseR( mat , dmat )
dd<-glmDenoiseR( mat, dmat, whichbase=c(1:2) , kfolds=4 , maxnoisepreds=10 , selectionthresh=0.1 ,collapsedesign=T, reestimatenoisepool=F )
# note - this may not be stable b/c the designmatrix is ill-conditioned
# dd<-glmDenoiseR( mat, dmat, whichbase=c(1:2) , kfolds=4 , maxnoisepreds=10 , selectionthresh=0.1 ,collapsedesign=F, reestimatenoisepool=F )
# now use its results in regression
glmdf<-data.frame( dd$hrfdesignmat, noiseu=dd$noiseu, polys=dd$polys )
mylm<-lm(  data.matrix(mat) ~ . , data=glmdf )
mylm<-bigLMStats( mylm , 0.001 )
rownames(mylm$beta.pval)
min(p.adjust(mylm$beta.pval[1,],method='n'))
min(p.adjust(mylm$beta.pval[2,],method='n'))
min(p.adjust(mylm$beta.pval[1,],method='BH'))
min(p.adjust(mylm$beta.pval[2,],method='BH'))
}
