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
#' # produce assumed hrf - bold2betas will use glmdenoiser
#' # internally to refine this assumed HRF across hrfShifts
#' # number of baseline shifts
#' b1<-hemodynamicRF( 12, onsets=1, durations=1, rt=1,
#'     cc=0.1,a1=3,a2=8,b1=0.5, b2=0.75 )
#' # looks ok so do the same w/in bold2betas
#' btsc<-bold2betas( boldmatrix=data.matrix(mat),  verbose=T,
#'       designmatrix=bb$desmat[,1:4], baseshift=0,
#'       blockNumb=runs, maxnoisepreds=2, hrfBasis=b1,
#'       hrfShifts=4, polydegree=4, selectionthresh=0.1 )
#' mylabs<-rep("",nrow(btsc$eventbetas))
#' for ( i in 1:nrow(btsc$eventbetas) ) mylabs[i]<-substr( rownames(btsc$eventbetas)[i],1,2)
#' mylabs<-as.factor(mylabs)
#' # sample for training
#' inds<-sample(1:nrow(btsc$eventbetas),size=round(nrow(btsc$eventbetas)*3./4.))
#' # basic voxel selection & classification
#' zz<-apply(btsc$eventbetas,FUN=mean,MARGIN=2)
#' zze<-zz/apply(btsc$eventbetas,FUN=sd,MARGIN=2)
#' th<-0.5
#' ff<-which( abs(zze) > th )
#' mydf<-data.frame( lab=mylabs,  vox=data.matrix(btsc$eventbetas)[,ff]) #, runs=runs[btsc$eventrows] )
#' mdl<-svm( lab ~., data=mydf[inds,])
#' err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
#' print(paste("NPredVox",length(ff),"Correct",err*100))
#' # cca voxel selection & classification
#' ccamats<-list( data.matrix(btsc$eventbetas)[inds,] ,
#'                data.matrix(bb$desmat[btsc$eventrows,1:4])[inds,] )
#' ##
#' ## this particular cca configuration selects orthogonal information
#' ## with initialization similar to that above.  however, the voxel
#' ## information is regularized and clustered first. each component
#' ## is also associated with particular "contrasts" (loosely speaking).
#' ## the result is 8 predictors associated with connected anatomical
#' ## regions where each predictor is an averaged beta map.
#' ## .... 2 init ideas below.
#' zze1<-zze2<-zze3<-zze4<-zze
#' zze1[ zze < th ]<-0
#' zze2[ zze < th*0.5 ]<-0
#' zze3[ zze > (th*(-1)) ]<-0
#' zze4[ zze > (th*(-0.5)) ]<-0
#' initcca<-t(cbind(zze1,zze2,zze3,zze4))
#' initcca<-t( svd( btsc$eventbetas, nu=0, nv=5 )$v )
#' initcca<-initializeEigenanatomy( initcca, mask=mask, nreps=2 )$initlist
#' nv<-length(initcca)
#' mycca<-sparseDecom2( inmatrix=ccamats, initializationList=initcca,
#'   sparseness=c( -0.001, -0.95 ), nvecs=nv, its=10, cthresh=c(250,0),
#'   uselong=0, smooth=0.0, mycoption=1, inmask=c(mask,NA) )
#' ccaout<-(data.matrix(imageListToMatrix( mycca$eig1, mask )))
#' ff<-which(colSums(abs(ccaout))>1.e-4)
#' mydf<-data.frame( lab=mylabs,
#'   vox=data.matrix(btsc$eventbetas) %*% t(ccaout) )
#' mdl<-svm( lab ~., data=mydf[inds,])
#' err<-sum(mydf[-inds,]$lab==predict( mdl, newdata=mydf[-inds,]))/nrow(mydf[-inds,])
#' print(paste("CCA:NPredVox",length(ff),"Correct",err*100))
#' vox<-antsImageClone(mask)
#' vox[mask==1]<-colSums(abs(ccaout))
#' antsImageWrite(vox,'temp.nii.gz')
#' 
bold2betas <- function( boldmatrix, designmatrixIn, blockNumb, maxnoisepreds=2, polydegree=4,
                        hrfShifts=20, hrfBasis=NA, crossvalidationgroups=4, selectionthresh=0.1,
                        multievents=FALSE, whichcols=NA, baseshift=0, tr=1, collapsedesign=TRUE,
                        verbose=FALSE )
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
# runhrfs<-data.frame(matrix( rep(0,length(allruns)*bl), ncol=bl))
eventhrfs<-data.frame(matrix( rep(0,neventstot*bl), ncol=bl))
eventrows<-rep(0,neventstot)
eventbetas<-data.frame(matrix( rep(0,neventstot*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
ct<-1
# for ( runs in allruns ) 
  {
#  if ( verbose ) print(paste("run%:",rct/length(allruns)*100,"event%:",(ct-1)/neventstot*100,"..."))
  kkt<-1:nrow(boldmatrix) # which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  submat<-boldmatrix[kkt,]
  dd<-glmDenoiseR( submat, denoisedes, hrfBasis=hrfBasis, selectionthresh=selectionthresh, tr=tr,
    crossvalidationgroups=blockNumb, maxnoisepreds=maxnoisepreds, hrfShifts=hrfShifts,
    collapsedesign=collapsedesign, reestimatenoisepool=F, polydegree = polydegree, baseshift=0 )
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
#  runhrfs[rct,] <- dd$hrf
  rct<-rct+1
  }# runs
return(list(eventbetas=eventbetas,eventhrfs=eventhrfs,
            eventrows=eventrows,hrf=colMeans(eventhrfs)))
}
