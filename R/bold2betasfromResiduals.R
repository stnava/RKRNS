bold2betasfromResiduals <- function( boldmatrix, designmatrix,
  blockNumb, hrfBasis, multievents=FALSE, whichcols=NA, 
  collapsedesign=TRUE, verbose=FALSE )
{
if ( all(is.na( whichcols )) ) whichcols<-1:ncol(designmatrix)
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
neventstot<-0
allruns<-unique( blockNumb )
for ( runs in allruns )
  {
  kkt<-which( blockNumb == runs )
  neventstot<-neventstot+sum(designmatrix[kkt,whichcols])
  }
designnames<-colnames(designmatrix)[whichcols]
if ( verbose ) print(designnames)
bl<-length(hrfBasis)
eventrows<-rep(0,neventstot)
eventbetas<-data.frame(matrix( rep(0,neventstot*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
ct<-1
  {
  kkt<-1:nrow(boldmatrix)
  denoisedes<-designmatrix[kkt,]
  submat<-boldmatrix[kkt,]
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
          twoeventmat[,1]<-conv(twoeventmat[,1],hrfBasis)[1:nrow(twoeventmat)]
          twoeventmat[,2]<-conv(twoeventmat[,2],hrfBasis)[1:nrow(twoeventmat)]
          glmdf<-data.frame( twoeventmat )
          mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
          mylm<-bigLMStats( mylm , 0.001 )
          eventbetas[ct,]<-mylm$beta.t[1,]
          eventrows[ct]<-rownames( denoisedes )[row]
          rownames(eventbetas)[ct]<-paste(designnames[col],eventrows[ct],sep='.')
          if ( verbose ) print(paste(rownames(eventbetas)[ct],"ct:",ct,'...',ct/neventstot*100,"%...Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          ct<-ct+1
      }
    }
  } # nevents > 0
  rct<-rct+1
  }# runs
return( list( eventbetas=eventbetas, eventrows=eventrows) )
}
