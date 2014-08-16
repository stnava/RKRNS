bold2betas <- function( boldmatrix, designmatrixIn, blockNumb, maxnoisepreds=2, polydegree=12,
                        hrfShifts=20, hrfBasis=NA, crossvalidationgroups=4, selectionthresh=0.1,
                        multievents=FALSE, whichcols=NA, baseshift=0, verbose=FALSE )
{
designmatrix<-ashift( designmatrixIn, c(baseshift,0) )
allruns<-unique( blockNumb )
if ( all(is.na( whichcols )) ) whichcols<-1:ncol(designmatrix)
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
neventstot<-0
for ( runs in allruns ) 
  {
  kkt<-which( blockNumb == runs )
  neventstot<-neventstot+sum(designmatrix[kkt,whichcols])
  }
designnames<-colnames(designmatrix)[whichcols]
if ( verbose ) print(designnames)
if ( ! all(is.na(hrfBasis) ) ) bl<-length(hrfBasis) else bl<-hrfShifts
runhrfs<-data.frame(matrix( rep(0,length(allruns)*bl), ncol=bl))
eventhrfs<-data.frame(matrix( rep(0,neventstot*bl), ncol=bl))
eventrows<-rep(0,neventstot)
eventbetas<-data.frame(matrix( rep(0,neventstot*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
ct<-1
for ( runs in allruns ) 
  {
  if ( verbose ) print(paste("run%:",rct/length(allruns)*100,"event%:",(ct-1)/neventstot*100,"..."))
  kkt<-which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  submat<-boldmatrix[kkt,]
  oneeventmat<-matrix( rowMeans(denoisedes), ncol=1 )
  dd<-glmDenoiseR( submat, denoisedes, hrfBasis=hrfBasis, selectionthresh=selectionthresh,
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, hrfShifts=hrfShifts,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree, baseshift=0 )
  if ( verbose ) plot( ts(dd$hrf) )
  glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
  glmdf<-data.frame( dd$hrfdesignmat, glmdfnuis )
  nevents<-sum(denoisedes[,whichcols])
  if ( nevents > 0 & TRUE ) {
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,whichcols] ) > 0 )
          {
          col<-which( denoisedes[row,whichcols] > 0 )
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
          twoeventmat<-cbind( denoisematmod1[,col] ,
                              rowSums(denoisematmod2) )
          twoeventmat[,1]<-conv(twoeventmat[,1],dd$hrf)[1:nrow(twoeventmat)]
          twoeventmat[,2]<-conv(twoeventmat[,2],dd$hrf)[1:nrow(twoeventmat)]
          glmdf<-data.frame( twoeventmat, glmdfnuis )
          mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
          mylm<-bigLMStats( mylm , 0.001 )
          eventbetas[ct,]<-mylm$beta.t[1,]
          eventrows[ct]<-rownames( denoisedes )[row]
          eventhrfs[ct,]<-dd$hrf
          rownames(eventbetas)[ct]<-paste(designnames[col],eventrows[ct],sep='.')
          if ( verbose ) print(paste(rownames(eventbetas)[ct],"ct:",ct,'...',ct/neventstot*100,"%...Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          ct<-ct+1
      }  
    }
  } # nevents > 0 
  runhrfs[rct,] <- dd$hrf
  rct<-rct+1
  }# runs
return(list(runhrfs=runhrfs,eventbetas=eventbetas,eventhrfs=eventhrfs,
            eventrows=eventrows,hrf=colMeans(runhrfs)))
}
