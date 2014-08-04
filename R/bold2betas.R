bold2betas <- function( boldmatrix, designmatrix, blockNumb, maxnoisepreds, polydegree=10, bl=50, crossvalidationgroups=4, selectionthresh=0.1, multievents=FALSE )
{
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
allruns<-unique( blockNumb )
for ( runs in allruns ) 
  {
  print(paste("ct",rct,"run",runs,":",rct/length(allruns)*100,"%"))
  kkt<-which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  denoisedes<-as.matrix( denoisedes[,colMeans(denoisedes)>0 ] )
  submat<-boldmatrix[kkt,]
  dd<-glmDenoiseR( submat, denoisedes, whichbase=3:3, selectionthresh=selectionthresh,
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree ) 
  glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
  glmdf<-data.frame( dd$hrfdesignmat, glmdfnuis )
  nevents<-sum(denoisedes)
  eventbetas<-data.frame(matrix( rep(0,nevents*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
  ct<-1
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,] ) > 0 )
          {
          col<-which( denoisedes[row,] > 0 )
          denoisematmod1<-denoisedes*0
          denoisematmod2<-denoisedes
          denoisematmod1[row,col]<-1
          denoisematmod2[row,col]<-0
          if ( multievents ) {
            rows<-which(denoisedes[,col]==1)
            rowsx<-which( (rows-row) > 0 )
            rows<-c( row, rows[rowsx[1]] )
            denoisematmod1[rows,col]<-1
            denoisematmod2[-rows,col]<-0
          }
          twoeventmat<-cbind( denoisematmod1[,col], rowSums(denoisematmod2))
          twoeventmat[,1]<-conv(twoeventmat[,1],dd$hrf)[1:nrow(twoeventmat)]
          twoeventmat[,2]<-conv(twoeventmat[,2],dd$hrf)[1:nrow(twoeventmat)]
          glmdf<-data.frame( twoeventmat, glmdfnuis )
          mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
          mylm<-bigLMStats( mylm , 0.001 )
          eventbetas[ct,]<-mylm$beta.t[1,]
#          rownames(eventbetas)[ct]<-paste(colnames(denoisedes)[col],".",rownames(boldmatrix)[row],sep='')
#          print(paste(rownames(eventbetas)[ct],"Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          print(paste(ct/nevents*100,"Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          ct<-ct+1
      }  
    }
  if ( is.na( fulleventbetas[1,1] ) ) fulleventbetas<-eventbetas else fulleventbetas<-rbind(fulleventbetas,eventbetas)
  rct<-rct+1
  }# runs
return(fulleventbetas)
}
