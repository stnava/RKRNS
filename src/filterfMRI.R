print("################## filter the fmri ##################")
imatf<-imat
# mycompcor<-compcor( data.matrix(imat), 2 )
if ( removeSentLengthEffects & !removeEventOverlap ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+nchar   ) ) )
if (!removeSentLengthEffects &  removeEventOverlap ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+eventss ) ) )
if ( removeSentLengthEffects &  removeEventOverlap ) imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ 0+eventss + nchar ) ) )
colnames(imatf)<-colnames(imat)
fh<-filterhighfrequency
fl<-filterlowfrequency
imatb<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" )) # , opt="trig" ))
colnames(imatf)<-colnames(imat)
for ( i in 1:ncol(imatf) ) {
  imatf[,i]<-data.frame(stl( ts( winsor(imatb[,i],winsorval),frequency=trendfrequency) , "per" )$time.series)$trend
}
globsig<-rowMeans( imat )
globsigb<-rowMeans( imatb )
globsigf<-rowMeans( imatf )
par(mfrow=c(2,1))
plotinds<-1:400
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
plotinds<-1:2000
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
# imatf<-imatb
