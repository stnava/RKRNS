print("################## filter the fmri ##################")
fl<-0.01
fh<-1  # because of expected bold response < 25secs, > 5 seconds
imatf<-imat
# mycompcor<-compcor( data.matrix(imat), 2 )
# imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ eventss + nchar  ) ) )
colnames(imatf)<-colnames(imat)
imatb<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" )) # , opt="trig" ))
colnames(imatf)<-colnames(imat)
for ( i in 1:ncol(imatf) ) {
  imatf[,i]<-data.frame(stl( ts( winsor(imatb[,i],0.01),frequency=2) , "per" )$time.series)$trend
}
globsigb<-rowMeans( imatb )
globsigf<-rowMeans( imatf )
par(mfrow=c(2,1))
plotinds<-2000:2100
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
plotinds<-2000:4000
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
