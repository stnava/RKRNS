simulateBOLD<-function(ntime=2000,nstim=30,signalscale=0.5,TR=0.5, lowfnoise=0.01 , physnoise=0.01, temporalnoise=0.01, option=c("basic","henson"), eximg=NA, mask=NA )
{
  if ( option[1] == "henson" ) {
    # from http://www.jstatsoft.org/v44/i10/paper
    dim <- c(53, 63, 46)
    nscan <- 384
    TR <- 2
    total.time <- nscan * TR
    n1on <- c(6.75, 15.75, 18, 27, 29.25, 31.5, 36, 42.75, 65.25,
      74.25, 92.25, 112.5, 119.25, 123.75, 126, 137.25, 141.75,
      144, 146.25, 155.25, 159.75, 162, 164.25, 204.75, 238.5)
    onsets.N1 <- n1on * TR
    n2on <- c(13.5, 40.5, 47.25, 56.25, 90, 94.5, 96.75, 135,
      148.5, 184.5, 191.25, 202.5, 216, 234, 236.25, 256.5, 261,
      281.25, 290.25, 303.75, 310.5, 319.5, 339.75, 342)
    onsets.N2 <- n2on * TR
    f1on <- c(0, 2.25, 9, 11.25, 22.5, 45, 51.75, 60.75, 63,
      76.5, 78.75, 85.5, 99, 101.25, 103.5, 117, 130.5, 150.75,
      171, 189, 227.25, 265.5, 283.5, 285.75, 288, 344.25)
    onsets.F1 <- f1on * TR
    f2on <- c(33.75, 49.5, 105.75, 153, 157.5, 168.75, 177.75,
      180, 182.25, 198, 222.75, 240.75, 254.25, 267.75, 270, 274.4,
      294.75, 299.25, 301.5, 315, 317.25, 326.25, 333, 335.25,
      337.5, 346.5)
    onsets.F2 <- f2on * TR

    designmat<-data.frame( N1=rep(0,nscan), N2=rep(0,nscan),
                           F1=rep(0,nscan), F2=rep(0,nscan),
                           Run=rep(1,nscan) )
    designmat$N1[ round(n1on)  ]<-1
    designmat$N2[ round(n2on)  ]<-1
    designmat$F1[ round(f1on) ]<-1
    designmat$F2[ round(f2on) ]<-1
    designmat$Run[ (round(nscan/2)+1):nscan ]<-2
    
    region.1A.center <- c(13, 13, 11)
    region.1A.radius <- 4
    region.1B.center <- c(40, 18, 9)
    region.1B.radius <- 6
    region.1C.center <- c(10, 45, 24)
    region.1C.radius <- 3
    region.2.center <- c(15, 16, 31)
    region.2.radius <- 5
    region.3.center <- c(12, 16, 13)
    region.3.radius <- 5
    onsets <- list(onsets.N1, onsets.N2, onsets.F1, onsets.F2)
    onsets.regions <- list(onsets, onsets, onsets, onsets, onsets)
    dur <- list(0, 0, 0, 0)
    dur.regions <- list(dur, dur, dur, dur, dur)
    region.1a.d <- list(160.46, 140.19, 200.16, 160.69)
    region.1b.d <- list(140.51, 120.71, 160.55, 120.44)
    region.1c.d <- list(120.53, 120.74, 140.02, 100.48)
    region.2.d <- list(-0.24, 10.29, 80.18, 160.24)
    region.3.d <- list(200.81, 50.04, 240.6, 50.83)
    effect <- list( region.1a.d, region.1b.d, region.1c.d, region.2.d,region.3.d )
# now use the mask
    design <- simprepTemporal(regions = 5, onsets = onsets.regions,
      durations = dur.regions, hrf = "double-gamma", TR = TR,
      totaltime = total.time, effectsize = effect)
    spatial <- simprepSpatial(regions = 5, coord = list(region.1A.center,
      region.1B.center, region.1C.center, region.2.center, region.3.center),
      radius = c(region.1A.radius, region.1B.radius, region.1C.radius,
      region.2.radius, region.3.radius), form = "sphere", fading = 0.01)
    if ( is.na(mask) ) mask<-getMask(eximg,cleanup=TRUE)
    sim.data <- simVOLfmri(design = design, image = spatial, base = as.array(eximg),
      SNR = 3.87, noise = "mixture", type = "rician", rho.temp = c(0.142,
      0.108, 0.084), rho.spat = 0.4, w = c(0.05, 0.1, 0.01,
      0.09, 0.05, 0.7), dim = dim(eximg), nscan = nscan, vee = 0,
      template = as.array(mask), spat = "gaussRF")
  antsimg<-as.antsImage(sim.data)
  antsimg<-antsImageClone( antsimg, 'float' )
  spc<-as.numeric(rep(TR,4))
  spc[1:3]<-antsGetSpacing(mask)
  antsSetSpacing( antsimg, spc )
  antsSetDirection(antsimg,diag(4))
  return( list(simbold=antsimg,desmat=designmat,mask=mask) )
  }
  rs<-2
  while ( max(rs) > 1 ) {
  simbold<-replicate(10000, rnorm(ntime)) # this is noise, 10,000 vox, ntime time points
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
  # nuisance
  if ( physnoise > 0 )
   {
   n.phys <-t(physnoise(dim = ncol(simbold), sigma = 15, nscan = nrow(simbold), TR = TR ) )
   simbold <- simbold + n.phys*temporalnoise
   }
  if ( lowfnoise > 0 )
   {
   n.low <- t(lowfreqdrift(dim = ncol(simbold), nscan = nrow(simbold), TR = TR, freq = 120))
   simbold <- simbold + n.low*lowfnoise
   }
  if ( temporalnoise > 0 )
   {
   n.temp <-t(temporalnoise(dim = ncol(simbold), sigma = 15, nscan = nrow(simbold), rho = c(0.4, -0.2)))
   simbold <- simbold + n.temp*temporalnoise
   }
#  print("n.task")
#  return( tasknoise(act.image = as.array(desmat), sigma = 15) )
  ts1<-ts(rowMeans(simbold[,380:600]))	
  plot( ts(data.frame(s3=s3,hrf3=hrf3,bold=ts1,boldraw=simbold[,390]) ))
  return(list(simbold=(simbold),desmat=desmat))
}
