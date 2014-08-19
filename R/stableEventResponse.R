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
