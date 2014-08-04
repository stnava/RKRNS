bold2betas <- function( boldmatrix, designmatrix, blockNumb, maxnoisepreds, polydegree=10, bl=40, crossvalidationgroups=4, selectionthresh=0.1 )
{
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
allruns<-unique( blockNumb ) 
for ( runs in allruns[1:4] ) 
  {
  print(paste("ct",rct,"run",runs,":",rct/length(allruns)*100,"%"))
  rct2<-rct+1
  rct3<-rct+2
  if ( rct2 == (length(allruns)+1) ) { rct2<-1; rct3<-2 }
  if ( rct3 == (length(allruns)+1) ) { rct3<-1 }
  kkt<-which( blockNumb == runs )
  denoisedes1<-designmatrix[kkt,]
#  denoisedes1<-as.matrix( denoisedes1[,colMeans(denoisedes1)>0 ] )
  submat<-boldmatrix[kkt,]
  dd<-glmDenoiseR( submat, denoisedes1, whichbase=1:6, 
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds,
                  selectionthresh=selectionthresh, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree ) 
  glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
  glmdf<-data.frame( dd$hrfdesignmat, glmdfnuis )

  kkt<-which( blockNumb == allruns[rct2]  )
  denoisedes<-designmatrix[kkt,]
#  denoisedes<-as.matrix( denoisedes[,colMeans(denoisedes)>0 ] )
  submat<-boldmatrix[kkt,]
  ee<-glmDenoiseR( submat, denoisedes, whichbase=1:6, 
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, 
                  selectionthresh=selectionthresh, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree ) 
  glmdfnuis2<-data.frame( noiseu=ee$noiseu, polys=ee$polys )
  glmdf2<-data.frame( ee$hrfdesignmat, glmdfnuis2 )

  kkt<-which( blockNumb == allruns[rct3]  )
  denoisedes<-designmatrix[kkt,]
#  denoisedes<-as.matrix( denoisedes[,colMeans(denoisedes)>0 ] )
  submat<-boldmatrix[kkt,]
  ee<-glmDenoiseR( submat, denoisedes, whichbase=1:6, 
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, 
                  selectionthresh=selectionthresh, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree ) 
  glmdfnuis3<-data.frame( noiseu=ee$noiseu, polys=ee$polys )
  glmdf3<-data.frame( ee$hrfdesignmat, glmdfnuis3 )

  glmdf<-rbind( glmdf, glmdf2, glmdf3 )
  glmdfnuis<-rbind( glmdfnuis, glmdfnuis2, glmdfnuis3 )
  
  kkt<-which( blockNumb == runs | blockNumb == allruns[rct2] | blockNumb == allruns[rct3] )
  denoisedes<-designmatrix[kkt,]
#  denoisedes<-as.matrix( denoisedes[,colMeans(denoisedes)>0 ] )
  eventbetas<-data.frame(matrix( rep(0,sum(denoisedes1)*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
  submat<-boldmatrix[kkt,]
  ct<-1
  for ( row in 1:nrow(denoisedes1) )
      {
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
          print(paste(ct,mean(mylm$beta.t[1,])))
          rownames(eventbetas)[ct]<-paste(colnames(denoisedes1)[col],".",rownames(boldmatrix)[row],sep='')
          print(paste(rownames(eventbetas)[ct],"Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          ct<-ct+1
          }  
      }
  if ( is.na( fulleventbetas[1,1] ) ) fulleventbetas<-eventbetas else fulleventbetas<-rbind(fulleventbetas,eventbetas)
  rct<-rct+1
  }# runs
return(fulleventbetas)
}
