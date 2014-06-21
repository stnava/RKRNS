print("################## filter the fmri ##################")
fl<-0.02
fh<-0.2 # because of expected bold response < 25secs, > 5 seconds
imatf<-imat
# mycompcor<-compcor( data.matrix(imat), 2 )
imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ eventss + nchar  ) ) )
colnames(imatf)<-colnames(imat)
imatb<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" )) # , opt="trig" ))
colnames(imatf)<-colnames(imat)
for ( i in 1:ncol(imatf) ) {
  imatf[,i]<-data.frame(stl( ts( winsor(imatb[,i],0.02),frequency=12) , "per" )$time.series)$trend
}
globsig<-rowMeans( imat )
globsigb<-rowMeans( imatb )
globsigf<-rowMeans( imatf )
par(mfrow=c(2,1))
plotinds<-2000:2400
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
plotinds<-2000:11000
plot( scale(globsig[plotinds]), type='l')
points( scale(globsigb[plotinds]), type='l', col='red')
points( scale(globsigf[plotinds]), type='l', col='green')
