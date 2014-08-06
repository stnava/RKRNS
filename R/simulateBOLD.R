simulateBOLD<-function(ntime=2000,nstim=30,signalscale=0.5)
{
  rs<-2
  while ( max(rs) > 1 ) {
  simbold<-replicate(10000, rnorm(ntime)) # this is noise, 10,000 vox, ntime time points
  # add linear and quadratic noise 
#  simbold<-simbold+rnorm(ntime)*c(1:ntime)/ntime*0.001
#  simbold<-simbold+rnorm(ntime)*sample(c(1:ntime)^2)/ntime*0.001
  # create 3 stimuli, a, b & c
  s1<-sample( c( rep(1,nstim), rep(0,ntime-nstim) ) ) # 10 repeats
  s2<-shift( s1, 150 ) # 10 repeats
  s3<-shift( s2, 150 ) # 10 repeats
  hrf1<-shift(hemodynamicRF( nrow(simbold), onsets=which(s1>0), durations=1, rt=1 ), 0 )
  hrf2<-shift(hemodynamicRF( nrow(simbold), onsets=which(s2>0), durations=1, rt=1 ), 0 )
  hrf3<-shift(hemodynamicRF( nrow(simbold), onsets=which(s3>0), durations=1, rt=1 ), 0 )
  simbold[,100:200]<-simbold[,100:200]+as.numeric(hrf1)*0.5*signalscale
  simbold[,320:400]<-simbold[,320:400]+as.numeric(hrf2)*0.4*signalscale
  simbold[,380:600]<-simbold[,380:600]+as.numeric(hrf3)*0.8*signalscale
  desmat<-data.frame( a=s1, b=s2, c=s3 ) 
  rs<-rowSums(desmat)
  }
  ts1<-ts(rowMeans(simbold[,380:600]))	
  plot( ts(data.frame(s3=s3,hrf3=hrf3,bold=ts1,boldraw=simbold[,390]) ))
  return(list(simbold=simbold,desmat=desmat))
}
