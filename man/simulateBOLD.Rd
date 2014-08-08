\name{simulateBOLD}
\alias{simulateBOLD}
\title{simulated bold study with controllable design and noise parameters}
\usage{
  bb<-simulateBOLD(option="basic")
  bb<-simulateBOLD(option="henson",eximg=ff)
}
\description{
  the henson example produces a mask, an image and a design matrix.
  the basic example produces a noisy data matrix and design matrix. 
  ... see args for some of the options and reasonable defaults.
}
\example{
# get example image
fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
eximg<-antsImageRead(fn,3)
fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="") 
mask<-antsImageRead(fn,3)
bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
# now test glmDenoise
mat<-timeseries2matrix( bb$simbold, bb$mask )
dd<-glmDenoiseR( mat, bb$desmat, crossvalidationgroups=2 , maxnoisepreds=0,
   hrfbasislength=35, selectionthresh=0.1 ,collapsedesign=T, polydegree=6 )
plot(ts(dd$hrf))
glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
glmdf<-data.frame( dd$hrfdesignmat, glmdfnuis )
mylm<-lm(  data.matrix(mat) ~ . , data=glmdf )
mylm<-bigLMStats( mylm , 1.e-4 )
print(rownames(mylm$beta.t))
for ( i in 1:nrow(mylm$beta.t)  ) 
  print(paste(rownames(mylm$beta.t)[i],":",max(abs(mylm$beta.t[i,]))))
vizimg<-antsImageClone(mask)
vizimg[mask==1]<-mylm$beta.t[3,]
plotANTsImage( eximg, functional=list(vizimg),
   slices='1x60x1', threshold="3x11", color='red')
}
