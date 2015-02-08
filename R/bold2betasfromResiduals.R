#' Convert a bold time series and design matrix to event-wise betas
#' 
#' Inspired by discussion with Kendrick Kay regarding his glm denoise tool
#' http://journal.frontiersin.org/Journal/10.3389/fnins.2013.00247/abstract.
#' Here, we apply glm denoise run by run to get event-wise betas.
#' 
#' 
#' @param boldmatrix input raw bold data in time by space matrix
#' @param designmatrix input design matrix - binary/impulse entries for event
#' related design, blocks otherwise
#' @param blockNumb numbers for the rows that should be treated together as
#' runs
#' @param maxnoisepreds number of noise predictors to explore e.g. a list of
#' values 1 to 10 or a few values to try
#' @param polydegree number of polynomial predictors
#' @param bl basis length for hrf estimation
#' @param crossvalidationgroups number of data splits in glm cross-validation
#' @return returns a list with relevant output
#' @author Avants BB
#' @examples
#' 
#' # get example image
#' fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="")
#' eximg<-antsImageRead(fn,3)
#' fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="")
#' mask<-antsImageRead(fn,3)
#' bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
#' boldImage<-bb$simbold
#' mat<-timeseries2matrix( bb$simbold, bb$mask )
#' runs<-bb$desmat$Run;
#' # produce assumed hrf - bold2betasfromResiduals
#' # expects that you used glmdenoiser to denoise all data
#' # then pass the residualized matrix to this function
#' hrf<-hemodynamicRF( 12, onsets=1, durations=1, rt=1,
#'     cc=0.1,a1=3,a2=8,b1=0.5, b2=0.75 )
#' dd<-glmDenoiseR( mat, bb$desmat[,1:4], hrfBasis=hrf, hrfShifts = 4 ,
#'   crossvalidationgroups=runs, maxnoisepreds=2 , selectionthresh=0.1 ,
#'   collapsedesign=T, polydegree=4, myintercept=0 )
#' # residualize bold matrix against all nuisance variables with a single model
#' rmat<-residuals( lm( mat ~ 0 + dd$polys + dd$noiseu ) )
#' # rmat<-dd$denoisedBold # denoised run by run data - should test both in practice
#' btsc<-bold2betasfromResiduals( boldmatrix=rmat,  verbose=T,
#'       designmatrix=bb$desmat[,1:4],
#'       blockNumb=runs,  hrfBasis=dd$hrf )
#' # now some decoding
#' mylabs<-rep("",nrow(btsc$eventbetas))
#' for ( i in 1:nrow(btsc$eventbetas) ) mylabs[i]<-substr( rownames(btsc$eventbetas)[i],1,2)
#' mylabs<-as.factor(mylabs)
#' # sample for training
#' inds<-sample(1:nrow(btsc$eventbetas),size=round(nrow(btsc$eventbetas)*3./4.))
#' # basic voxel selection & classification
#' zz<-apply(btsc$eventbetas,FUN=mean,MARGIN=2)
#' zze<-zz/apply(btsc$eventbetas,FUN=sd,MARGIN=2)
#' th<-0.6
#' ff<-which( abs(zze) > th )
#' mydf<-data.frame( lab=mylabs,  vox=data.matrix(btsc$eventbetas)[,ff]) #, runs=runs[btsc$eventrows] )
#' mdl<-svm( lab ~., data=mydf[inds,])
#' err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
#' print(paste("NPredVox",length(ff),"Correct",err*100))
#' 
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
