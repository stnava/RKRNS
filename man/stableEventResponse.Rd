\name{stableEventResponse}
\alias{stableEventResponse}
\title{Denoise and estimate response variability across runs}
\description{
Use extended glmdenoise and leave-a-run-out to summarize 
response variability across runs for an event.  Returns a
matrix of betas where there is a beta vector for each left 
out run.  The betas are t-tstats for the selected event.
}
\usage{
eventbetas<-stableEventResponse( boldmatrix, designmatrix, 
    blockNumb, eventselect='coffee', polydegree=10, bl=50,
    crossvalidationgroups=4, selectionthresh=0.1, maxnoisepreds=4 )
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
\item{eventselect}{
  selection criteria determining which events are estimated
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
returns a list of cross-validated betas 
}

\author{
Avants BB
}

\examples{
  bsim<-simulateBOLD(signalscale=5)
  imat<-bsim$simbold
  dmats<-bsim$desmat
  runl<-nrow(imat)/4
  fakeruns<-c( rep(1,runl) , rep(2,runl), rep(3,runl) )
  fakeruns<-c( fakeruns, rep(4,nrow(imat)-length(fakeruns)))
  st1<-stableEventResponse(imat,dmats, fakeruns, 
    eventselect="a", maxnoisepreds=0, polydegree=1, bl=30 )
  st2<-stableEventResponse(imat,dmats, fakeruns, 
    eventselect="b", maxnoisepreds=0, polydegree=1, bl=30 )
  plot(ts(colMeans(st1$cveventbetas[,1:1000])))
  plot(ts(colMeans(st2$cveventbetas[,1:1000])))

  nrows<-nrow(imat)
  samp1<-1:round(nrows*0.6)
  samp2<-(max(samp1)+1):nrows
  soldier<-grep('soldier.kicked',colnames(dmats))
  imatd<-data.matrix(imat)
  rm(imat)
  imatd<-winsor(imatd,0.01)
  imatd<-imatd/colMeans(imatd)*100.0-100.0
  i1<-imatd[samp1,] # now % signal change
  i2<-imatd[samp2,] # now % signal change
  rm(imatd)
  ee1<-stableEventResponse(i1,dmats[samp1,], dmat$blockNumb[samp1], 
    eventselect=soldier, maxnoisepreds=0, polydegree=4 )
  dmf<-data.matrix(ee1)
  eanat<-sparseDecom( inmatrix=dmf, inmask=subaal, 
    sparseness=0.01, nv, cthresh=20, its=5 )
  eanatmat<-imageListToMatrix( eanat$eig, subaal )
  for ( k in 1:nv) 
    {
    whichvox<-( abs(eanatmat[k,])>1.e-4 )
    locor<-cor(t(dmf[,whichvox]))
    mm<-apply(dmf[,whichvox],FUN=mean,MARGIN=2)
    ms<-apply(dmf[,whichvox],FUN=sd,MARGIN=2)
    mn<-apply(abs(dmf[,whichvox]),FUN=min,MARGIN=2)
    effsz<-mean( mm/ms )
    print(paste(k,":",sum(whichvox),mean(locor),'e:',effsz,'t:',mean(mm),'min:',mean(mn)))
#   pheatmap(cor(t(dmf[,whichvox])))
    }

ee<-stableEventResponse(imat,dmats, dmat$blockNumb, 
  eventselect=soldier,
  maxnoisepreds=14, polydegree=4 )
ff<-stableEventResponse(imat,dmats, dmat$blockNumb, 
  eventselect='coffee.was.hot',
  maxnoisepreds=14, polydegree=4 )
nv<-20
eanatf<-sparseDecom( inmatrix=data.matrix(ff), inmask=subaal, 
  sparseness=0.01, nv, cthresh=20, its=5 )
eanatg<-sparseDecom( inmatrix=data.matrix(ee), inmask=subaal, 
  sparseness=0.01, nv, cthresh=20, its=5 )
eanatmatf<-imageListToMatrix( eanatf$eig, subaal )
eanatmatg<-imageListToMatrix( eanatg$eig, subaal )
for ( k in 1:nv) 
  {
  whichvox<-( abs(eanatmat[k,])>1.e-4 )
  locor<-cor(t(dmf[,whichvox]))
  mm<-apply(dmf[,whichvox],FUN=mean,MARGIN=2)
  ms<-apply(dmf[,whichvox],FUN=sd,MARGIN=2)
  mn<-apply(abs(dmf[,whichvox]),FUN=min,MARGIN=2)
  effsz<-mean( mm/ms )
  print(paste(k,":",sum(whichvox),mean(locor),'e:',effsz,'t:',mean(mm),'min:',mean(mn)))
  }
}
