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
#  betas<-bold2betas( imat, dmats, dmat$blockNumb, maxnoisepreds=c(12,14,18) )
# example of what you might do for decoding
# for real data
#  ee2<-bold2betas( data.matrix(imat[samp1,]) , dmats[samp1,], dmat$blockNumb[samp1], bl=30, polydegree=8, whichcols=whichcols ) 
# but here make some simulated data
  bsim<-simulateBOLD()
  simbold<-bsim$simbold
  desmat<-bsim$desmat
  ee2<-bold2betas( simbold , desmat, rep(1,nrow(simbold)), bl=20, polydegree=1, maxnoisepreds=1:10 ) 
  ee2$eventbe[nrow(ee2$eventbe),]<-colMeans(  ee2$eventbe[1:(nrow(ee2$eventbe)-1),] )	   	   
  plot(ts(colMeans(ee2$eventbe[,1:2000])))
  classes<-rep("",nrow(ee2$eventbetas))
  for ( i in 1:nrow(ee2$eventbetas)) classes[i]<-strsplit(rownames(ee2$eventbetas)[i],"[.]")[[1]][1]
  classes<-as.factor(classes)

  subvox<-( abs(colMeans(ee2$eventbetas)) > mean(abs(colMeans(ee2$eventbetas)))*2 )
  submat<-data.matrix(ee2$eventbetas[,])
  mysvd<-svd(submat,nu=4)$u 
  eanat<-sparseDecom( submat, nvecs=30, its=1, sparseness=-0.01, mycoption=0, cthresh=0) # 0, z=-0.5 )
  projmat<-data.matrix( eanat$eigenanatomyimages )
  rownames(eanat$projections)<-rownames(ee2$eventbetas)
  inds<-sample(1:nrow(ee2$eventbetas))[1:round(length(classes)*3/4)]
#  mydf<-data.frame( lab=classes, p=submat[,abs(projmat[,1])>0|abs(projmat[,2])>0|abs(projmat[,3])>0] )
#  mydf<-data.frame( lab=classes, p=submat %*% projmat )
  mydf<-data.frame( lab=classes, p=mysvd )
  errrr<-0; nt<-20
  for ( test in 1:nt ) {
    inds<-sample(1:nrow(ee2$eventbetas))[1:round(length(classes)*3/4)]
    mdl<-randomForest(lab~., mydf[inds,] )
    predres<-predict( mdl, newdata=mydf[-inds,]) 
    names(predres)<-paste(mydf$lab[-inds],'true')
    err<-sum(predres==mydf$lab[-inds])/length(predres)
    errrr<-errrr+err/nt
  }
  print(errrr)
# now cca 

  inds<-sample(1:nrow(betmat))[1:round(length(classes)*7/10)]

  betmat<-t(scale( t(data.matrix(ee2$eventbetas) ) ))
  submat<-t(scale(t(submat)))
  classmat<-data.matrix(cbind(rep(0,nrow(submat)),rep(0,nrow(submat)),rep(0,nrow(submat))))
  for (k in 1:length(levels(classes))) classmat[,k]<-classes[inds]==levels(classes[inds])[k]
  eanat<-sparseDecom2( list(submat,classmat), nvecs=5, its=20, sparseness=c( -0.01, -0.1 ), mycoption=0 ) 
  projmat<-(data.matrix( eanat$eig1 ))
  whvox <- abs(projmat[,1])>0 | abs(projmat[,2])>0  | abs(projmat[,3])>0| abs(projmat[,4])>0
  whvox <- abs(projmat[,1])>0  # | abs(projmat[,2])>0 # | abs(projmat[,3])>0
  mydfcca<-data.frame( lab=classes, p=betmat[,whvox] )
  mydfcca<-data.frame( lab=classes, p=betmat %*% projmat )
  mdlcca<-svm(lab~., mydfcca[inds,] )
  predres<-predict( mdlcca, newdata=mydfcca[-inds,]) 
  names(predres)<-paste(mydfcca$lab[-inds],'true')
  err1<-sum(predres==mydfcca$lab[-inds])/length(predres)
  print(err1)

# compare to a sort of ground truth 
  mydf<-data.frame( lab=classes, p=svd(ee2$eventbetas[,c(100:600)],nu=10 )$u  )
  mdl<-svm(lab~., mydf[inds,] )
  predres<-predict( mdl, newdata=mydf[-inds,]) 
  names(predres)<-paste(mydf$lab[-inds],'true')
  err2<-sum(predres==mydf$lab[-inds])/length(predres)
  print(paste(err1,err2))

}
