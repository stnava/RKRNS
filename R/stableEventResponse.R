stableEventResponse <- function( boldmatrix, designmatrix, blockNumb, eventselect='coffee',
      polydegree=10, bl=50, crossvalidationgroups=4, selectionthresh=0.1, maxnoisepreds=4 )
{
eselect<-eventselect
if ( all( is.character(eventselect) ) ) eselect<-grep( eventselect, colnames(designmatrix) )
# otherwise treat as numeric
if (length(eselect)==1) stim<-designmatrix[,eselect] else {
  stim<-rowSums(designmatrix[,eselect])
}
notstim<-rowSums(designmatrix[,-eselect])
# now subset stim by bl
whrows<-which( stim == 1 )
debug<-F
if ( debug ) print(whrows)
# cross val over whrows - leaving out a row ...
# est betas at each, then 
whtimes<-c()
maxrow<-nrow(boldmatrix)
for ( i in 1:length(whrows) )
  {
  lo<-whrows[i]-bl
  hi<-whrows[i]+bl
  if ( lo < 1 ) lo<-1
  if ( hi > maxrow ) hi<-maxrow
  if ( i == 1 ) whtimes<-c( lo:hi ) else whtimes<-c(whtimes,lo:hi)
  }
whtimes<-unique(whtimes)
if ( debug ) print(whtimes)
if ( debug ) print(maxrow)
stim<-stim[whtimes]
notstim<-notstim[whtimes]
des<-cbind(stim=stim,notstim=notstim)
colnames(des)[1]<-eventselect
colnames(des)[2]<-paste("Not.",eventselect,sep='')
# estimate hrf --- add run as a factor ....
runfactor=blockNumb[whtimes]
runf<-as.factor(runfactor)
if (var(as.numeric(runf))==0) return(NA)
dd<-glmDenoiseR( boldmatrix[whtimes,], des, whichbase=NA, selectionthresh=selectionthresh,
    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, hrfbasislength=bl,
    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree, runfactor=runf )

# get hrf & noise vecs - estimate beta
des[,1]<-conv(des[,1],dd$hrf)[1:nrow(des)]
des[,2]<-conv(des[,2],dd$hrf)[1:nrow(des)]
if ( debug ) plot( ts(dd$hrf) )
glmdfnuis<-data.frame( noiseu=dd$noiseu, polys=dd$polys )
allruns<-unique( blockNumb[whtimes] )
nevents<-length(allruns)
eventbetas<-data.frame(matrix( rep(0,nevents*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
rct<-1
for ( runs in allruns )
  {
  print(paste("ct",rct,"run",runs,":",rct/length(allruns)*100,"%"))
  kkt<-which( blockNumb[whtimes] != runs )
  print(kkt)
  submat<-boldmatrix[kkt,]
  glmdf<-data.frame( des[kkt,], glmdfnuis[kkt,] )
  mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
  mylm<-bigLMStats( mylm , 0.001 )
  eventbetas[rct,]<-mylm$beta.t[1,]
  rownames(eventbetas)[rct]<-paste(runs,sep='.')
  rct<-rct+1
  }
return( list(cveventbetas=eventbetas, glmdenoiz=dd ) )
#betas<-bold2betas( boldmatrix[whtimes,], des, blockNumb[whtimes], 
#                  maxnoisepreds=12, polydegree=4 )
# return(betas)
# for ( run in unique( blockNumb[whtimes] ) ) {
#  dd<-glmDenoiseR( boldmatrix[whtimes,], des2, whichbase=NA, selectionthresh=selectionthresh,
#    crossvalidationgroups=crossvalidationgroups , maxnoisepreds=maxnoisepreds, hrfbasislength=bl,
#    collapsedesign=T, reestimatenoisepool=F, polydegree = polydegree )
# }
}
