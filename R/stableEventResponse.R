#' Estimate effect sizes per-voxel, per stimulus, cross-validated
#' 
#' Use leave-a-run-out to summarize response variability across runs for a
#' stimulus class.  Returns a matrix of beta effect sizes per stimulus class
#' where the beta effect size are defined as mean of the beta response across
#' left out runs divided by standard deviation of the same.
#' 
#' 
#' @param boldmatrix input raw bold data in time by space matrix
#' @param designmatrix input design matrix - binary/impulse entries for event
#' related design, blocks otherwise
#' @param runIDs numbers for the rows that should be treated together as runs
#' @param polydegree number of polynomial predictors
#' @param hrf input hrf to use in regression
#' @param baseshift basis shift for design matrix
#' @return returns a list of cross-validated effect sizes for different stimuli
#' classes
#' @author Avants BB
#' @examples
#' 
#' fn<-paste(path.package("RKRNS"),"/extdata/111157_mocoref_masked.nii.gz",sep="") 
#' eximg<-antsImageRead(fn,3)
#' fn<-paste(path.package("RKRNS"),"/extdata/subaal.nii.gz",sep="") 
#' mask<-antsImageRead(fn,3)
#' bb<-simulateBOLD(option="henson",eximg=eximg,mask=mask)
#' boldImage<-bb$simbold
#' mat<-timeseries2matrix( bb$simbold, bb$mask )
#' runs<-bb$desmat$Run; 
#' hrf<-hemodynamicRF( 20, onsets=1, durations=1, rt=1,cc=0.1 )
#' stb<-stableEventResponse(mat,  bb$desmat[,1:4], runs,  hrf=hrf )
#' var1i<-antsImageClone( mask ) 
#' var1i[mask==1]<-stb[1,] # then write this out ...
#' # or run eigseg 
#' stb2<-stb
#' stb2[ stb < 6 ]<-0
#' ee<-eigSeg(mask, matrixToImages( stb2, mask) )
#' ImageMath(3,ee,'ClusterThresholdVariate',ee,mask,5)
#' # antsImageWrite ...
#' 
stableEventResponse <- function( boldmatrix, designmatrixIn, runIDs, hrf,
                        verbose=F, polydegree=4, baseshift=0, timevals=NA )
{
# for each event class, estimate betas for each fold based
# on leave one run out cross-validation
# return statistics on the stability of beta responses
# i.e. find voxels that have good, stable beta values
if ( var( runIDs ) == 0 )
  {
  print("Need distinct runIDs")
  return(NA)
  }
designmatrix<-ashift( designmatrixIn, c(baseshift,0) )
whichcols<-1:ncol(designmatrix)
allruns<-unique( runIDs )
neventstot<-0
for ( runs in allruns ) 
  {
  kkt<-which( runIDs == runs )
  neventstot<-neventstot+sum(designmatrix[kkt,whichcols])
  }
designnames<-colnames(designmatrix)[whichcols]
if ( verbose ) print(designnames)
if ( all(is.na(timevals)) ) timevals<-1:nrow(designmatrix)
p<-stats::poly( timevals ,degree=polydegree )
runf<-as.factor( runIDs )
betastab<-matrix( rep(0,length(whichcols)*ncol(boldmatrix)),
                      ncol=ncol(boldmatrix) )
colct<-1
for ( mycol in whichcols )
  {
  rct<-1
  eventbetas<-matrix( rep(0,length(allruns)*ncol(boldmatrix)),
                      ncol=ncol(boldmatrix) )
  for ( runs in allruns ) 
    {
    if ( verbose ) print(paste("mycol",mycol,"run%:",rct/length(allruns)*100))
    kkt<-runIDs != runs 
    twoeventmat<-cbind( designmatrix[kkt, mycol] ,
                 rowSums(designmatrix[kkt,-mycol]) )
    twoeventmat[,1]<-conv(twoeventmat[,1],hrf)[1:nrow(twoeventmat)]
    twoeventmat[,2]<-conv(twoeventmat[,2],hrf)[1:nrow(twoeventmat)]
    glmdf<-data.frame( twoeventmat, p=p[kkt,], runf=runf[kkt] )
    mylm<-lm(  data.matrix(boldmatrix[kkt,]) ~ . , data=glmdf )
    mylm<-bigLMStats( mylm , 0.001 )
    eventbetas[rct,]<-mylm$beta.t[1,]
    rct<-rct+1
    }
#  betalist<-lappend( betalist , eventbetas )
  betastab[colct,]<-apply( eventbetas, FUN=mean, MARGIN=2)/
                    apply( eventbetas, FUN=sd, MARGIN=2)
  colct<-colct+1
  }# runs
return(betastab)
}
