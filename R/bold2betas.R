bold2betas <- function( boldmatrix, designmatrix, blockNumb, maxnoisepreds, polydegree=10, bl=40, crossvalidationgroups=4 )
{
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
allruns<-unique( blockNumb ) 
for ( runs in allruns[1] ) 
  {
  print(paste("ct",rct,"run",runs,":",rct/length(allruns)*100,"%"))
  kkt<-which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  denoisedes<-as.matrix( denoisedes[,colMeans(denoisedes)>0 ] )
  submat<-boldmatrix[kkt,]
  dd<-glmDenoiseR( submat, denoisedes, whichbase=1:6, 
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, selectionthresh=0.2, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree ) 
  glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
  glmdf<-data.frame( dd$hrfdesignmat, glmdfnuis )
  eventbetas<-data.frame(matrix( rep(0,sum(denoisedes)*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
  eventclass<-rep(0,sum(denoisedes))
  ct<-1
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,] ) > 0 )
          {
          col<-which( denoisedes[row,] > 0 )
          denoisematmod1<-denoisedes*0
          denoisematmod2<-denoisedes
          denoisematmod1[row,col]<-1
          denoisematmod2[row,col]<-0 
          twoeventmat<-cbind( denoisematmod1[,col], rowSums(denoisematmod2))
          twoeventmat[,1]<-conv(twoeventmat[,1],dd$hrf)[1:nrow(twoeventmat)]
          twoeventmat[,2]<-conv(twoeventmat[,2],dd$hrf)[1:nrow(twoeventmat)]
          glmdf<-data.frame( twoeventmat, glmdfnuis )
          mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
          mylm<-bigLMStats( mylm , 0.001 )
          eventbetas[ct,]<-mylm$beta.t[1,]
          eventclass[ct]<-col
          rownames(eventbetas)[ct]<-paste(colnames(denoisedes)[col],".",rownames(boldmatrix)[row],sep='')
          print(paste(rownames(eventbetas)[ct],"Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          ct<-ct+1
      }  
    }
  if ( is.na( fulleventbetas[1,1] ) ) fulleventbetas<-eventbetas else fulleventbetas<-rbind(fulleventbetas,eventbetas)
  rct<-rct+1
  }# runs
return(fulleventbetas)
}
