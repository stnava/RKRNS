\name{bold2betas}
\alias{bold2betas}
\title{Convert a bold time series and design matrix to event-wise betas}
\description{
Inspired by discussion with Kendrick Kay regarding his glm denoise tool 
http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00247/abstract.
Here, we apply glm denoise run by run to get event-wise betas.
}
\usage{
  eventbetas<-bold2betas( boldmatrix, designmatrix, blockNumb, maxnoisepreds, polydegree=10, bl=40, crossvalidationgroups=4  )
}

\arguments{
\item{boldmatrix}{
  input raw bold data in time by space matrix 
}
\item{designmatrix}{
  input design matrix - binary/impulse entries for event related design, blocks otherwise
}
\item{blockNumb}{
  numbers for the rows that should be treated together as runs 
}
\item{maxnoisepreds}{
  number of noise predictors to explore e.g. a list of values 1 to 10 or a few values to try
}
\item{polydegree}{
  number of polynomial predictors
}
\item{bl}{
  basis length for hrf estimation
}
\item{crossvalidationgroups}{
  number of data splits in glm cross-validation
}
}

\value{
returns a list with relevant output 
}

\author{
Avants BB
}

\examples{
# get example image
fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
eximg<-antsImageRead(fn,3)
fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="") 
mask<-antsImageRead(fn,3)
bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
boldImage<-bb$simbold
mat<-timeseries2matrix( bb$simbold, bb$mask )
runs<-bb$desmat$Run; 
# use glm denoise to estimate a FIR HRF
glmda<-glmDenoiseR( data.matrix(mat), bb$desmat[,1:4], hrfBasis=NA, 
  hrfShifts = 12, maxnoisepreds=0 , selectionthresh=0.1 , 
  collapsedesign=T, polydegree=12, baseshift=0 )
# use glm denoise to refine an assumed HRF
b1<-hemodynamicRF( 12, onsets=1, durations=1, rt=1,
    cc=0.1,a1=3,a2=8,b1=0.5, b2=0.75 )
glmda<-glmDenoiseR( data.matrix(mat), bb$desmat[,1:4], hrfBasis=b1, 
  hrfShifts = 12, maxnoisepreds=0 , selectionthresh=0.1, 
  collapsedesign=T, polydegree=4 )
# looks ok so do the same w/in bold2betas
btsc<-bold2betas( boldmatrix=data.matrix(mat),  verbose=T,
      designmatrix=bb$desmat[,1:4], baseshift=0,
      blockNumb=runs, maxnoisepreds=2, hrfBasis=glmda$hrf,
      hrfShifts=0, polydegree=4, selectionthresh=0.1 )
mylabs<-rep("",nrow(btsc$eventbetas))
for ( i in 1:nrow(btsc$eventbetas) ) mylabs[i]<-substr( rownames(btsc$eventbetas)[i],1,2)
mylabs<-as.factor(mylabs)
# sample for training
inds<-sample(1:nrow(btsc$eventbetas),size=round(nrow(btsc$eventbetas)*3./4.))
# basic voxel selection & classification
zz<-apply(btsc$eventbetas,FUN=mean,MARGIN=2)
zze<-zz/apply(btsc$eventbetas,FUN=sd,MARGIN=2)
th<-0.5
ff<-which( abs(zze) > th )
mydf<-data.frame( lab=mylabs,  vox=data.matrix(btsc$eventbetas)[,ff]) #, runs=runs[btsc$eventrows] )
mdl<-svm( lab ~., data=mydf[inds,])
err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
print(paste("NPredVox",length(ff),"Correct",err*100))
# cca voxel selection & classification
ccamats<-list( data.matrix(btsc$eventbetas)[inds,] , 
               data.matrix(bb$desmat[btsc$eventrows,1:4])[inds,] )
##
## this particular cca configuration selects orthogonal information
## with initialization similar to that above.  however, the voxel
## information is regularized and clustered first. each component 
## is also associated with particular "contrasts" (loosely speaking).  
## the result is 8 predictors associated with connected anatomical
## regions where each predictor is an averaged beta map.
## .... 2 init ideas below.
zze1<-zze2<-zze3<-zze4<-zze
zze1[ zze < th ]<-0
zze2[ zze < th*0.5 ]<-0
zze3[ zze > (th*(-1)) ]<-0
zze4[ zze > (th*(-0.5)) ]<-0
initcca<-t(cbind(zze1,zze2,zze3,zze4))
initcca<-t( svd( btsc$eventbetas, nu=0, nv=5 )$v )
initcca<-initializeEigenanatomy( initcca, mask=mask, nreps=2 )$initlist
nv<-length(initcca)
mycca<-sparseDecom2( inmatrix=ccamats, initializationList=initcca,
  sparseness=c( -0.001, -0.95 ), nvecs=nv, its=10, cthresh=c(250,0),
  uselong=0, smooth=0.0, mycoption=1, inmask=c(mask,NA) )
ccaout<-(data.matrix(imageListToMatrix( mycca$eig1, mask )))
ff<-which(colSums(abs(ccaout))>1.e-4)
mydf<-data.frame( lab=mylabs,  
  vox=data.matrix(btsc$eventbetas) %*% t(ccaout) )
mdl<-svm( lab ~., data=mydf[inds,])
err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
print(paste("CCA:NPredVox",length(ff),"Correct",err*100))
vox<-antsImageClone(mask)
vox[mask==1]<-colSums(abs(ccaout))
antsImageWrite(vox,'temp.nii.gz')
}
