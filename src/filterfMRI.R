print("################## filter the fmri ##################")
fl<-0.05
fh<-1
globsig<-rowMeans( imat )
mycompcor<-compcor( data.matrix(imat), 2 )
imatf<-data.frame(residuals( lm(  data.matrix(imat) ~ eventss + nchar  ) ) ) #+ mycompcor)
colnames(imatf)<-colnames(imat)
imatf<-data.frame(frequencyFilterfMRI( imatf,  tr=tr, freqLo=fl, freqHi=fh, opt="butt" ))
colnames(imatf)<-colnames(imat)
dmatsblock<-dmats
blocksize<-20
for ( i in 1:blocksize ) dmatsblock<-dmatsblock+ashift(dmats,c(i,0))
for ( i in 1:ncol(imatf) ) {
  imatf[,i]<-data.frame(stl( ts( winsor(imatf[,i],0.01),frequency=100) , "per" )$time.series)$trend
}
