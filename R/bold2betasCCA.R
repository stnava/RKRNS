bold2betasCCA <- function( boldmatrix, designmatrix, blockNumb, bl=12, baseshift=5, mask=NA, sparseness=c(0,0),
                          multievents=FALSE, polydegree=10, bestvoxnum=50, nvecs=5, whichcols=NA,
                          mycoption=1, its=10, onlyhrf=FALSE )
{
    # first use FIR and lm to estimate hrf and betas for polynomials and stimuli
    # second estimate noise regressors from the voxels that relate highly to the poly regressors
    # third estimate a new hrf from CCA variables along with voxels/responses
par(mfrow=c(1,2))
rct<-1
if ( all(is.na( whichcols )) ) whichcols<-1:ncol(designmatrix)
allruns<-unique( blockNumb )
print(allruns)
neventstot<-0
for ( runs in allruns ) 
  {
  kkt<-which( blockNumb == runs )
  neventstot<-neventstot+sum(designmatrix[kkt,whichcols])
  }
print( neventstot )
designnames<-colnames(designmatrix)[whichcols]
runhrfs<-data.frame(matrix( rep(0,length(allruns)*bl), ncol=bl))
eventhrfs<-data.frame(matrix( rep(0,neventstot*bl), ncol=bl))
eventrows<-rep(0,neventstot)
eventbetas<-data.frame(matrix( rep(0,neventstot*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
ct<-1
for ( runs in allruns ) 
  {
  print(paste("run%:",rct/length(allruns)*100,"event%:",(ct-1)/neventstot*100,"..."))
  kkt<-which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  submat<-boldmatrix[kkt,]
  nevents<-sum(denoisedes[,whichcols])
  if ( nevents > 0 )
  {
  allevents<-matrix( rowSums( denoisedes ), ncol=1 )
  alleventsfir<-finiteImpuleResponseDesignMatrix( allevents, n=bl, baseshift=baseshift )
  p<-stats::poly( 1:nrow(alleventsfir) ,degree=polydegree )
  locsp<-c(bestvoxnum/ncol(submat), -0.9 )
#  compcornuis<-compcor( submat, 2 )
  mycca<-sparseDecom2( inmatrix=list( data.matrix(residuals(lm( submat~0+p))), alleventsfir ),
                       sparseness=locsp, nvecs=nvecs, its=its, cthresh=c(bestvoxnum,0),
                       uselong=0, smooth=0, mycoption=mycoption, inmask=c(mask,NA) )
  mytype<-typeof( mycca$eig1[[1]] )
  if ( mytype == "double" ) ccamat<-t(mycca$eig1) else ccamat<-imageListToMatrix( mycca$eig1, mask )
  if ( mean( data.matrix(mycca$eig2) ) < 0 ) mycca$eig2<-mycca$eig2*(-1)
  myhrf2<-as.numeric(mycca$eig2[1:bl,1])
  meanhrfval<-mean(myhrf2)
  mxdf<-abs(max(myhrf2)-meanhrfval)
  mndf<-abs(min(myhrf2)-meanhrfval)
  if ( mndf > mxdf ) myhrf2<-myhrf2*(-1)
  runhrf<-myhrf2/max(myhrf2)
  runhrf<-data.frame(stl(ts(runhrf, frequency = 2),"per")$time.series)$trend
  plot(ts(runhrf))
  if ( ! onlyhrf ) {      
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,whichcols] ) > 0 )
          {
          myhrf<-runhrf
          col<-which( denoisedes[row,whichcols] > 0 )
          denoisematmod1<-matrix( rep(0,nrow(denoisedes)), ncol=1)
          denoisematmod2<-denoisedes
          denoisematmod1[row,1]<-1
          denoisematmod2[row,col]<-0
          # now get the hrf estimating part 
          twoeventmat<-cbind( denoisematmod1, rowSums(denoisematmod2))
#          twoeventmat<-finiteImpuleResponseDesignMatrix( twoeventmat, n=bl, baseshift=baseshift )
#          twoeventmat<-finiteImpuleResponseDesignMatrix( denoisematmod2, n=bl, baseshift=baseshift )
          p<-stats::poly( 1:nrow(twoeventmat) ,degree=polydegree )
          doconvolution<-TRUE
          if ( doconvolution ) 
              for ( tk in 1:ncol(twoeventmat) )
                  twoeventmat[,tk]<-conv( twoeventmat[,tk], runhrf )[1:length(twoeventmat[,tk])]
          glmdf<-data.frame( twoeventmat , p=p )
          mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
          mylm<-bigLMStats( mylm , 0.001 )
          initmat<-mylm$beta.t[1:(ncol(twoeventmat)),]
          pblock<-mylm$beta.t[(ncol(twoeventmat)+1):nrow(mylm$beta.t),]
          pblock<-apply( (pblock) , FUN=sum, MARGIN=2)
#          initmat<-rbind(initmat,pblock)
          initdf<-initializeEigenanatomy( initmat , mask=mask, 2 )
          initmat<-rbind( runhrf, runhrf )
          initdf2<-initializeEigenanatomy( initmat , mask=NA, 2 )
          mask=initdf$mask
          nvecs<-length( initdf$initlist )
          locdes<-data.matrix(glmdf)
          noisevox<-abs(pblock)>mean(abs(pblock)) # & abs(temp)<mean(abs(temp))
#          ccamat<-residuals( lm(  data.matrix(submat) ~ p + 1 ) )
          ccamatin<-residuals( lm(  data.matrix(submat) ~ p + twoeventmat[,2] ) )
#          noisesvd<-svd( ccamat , nu=8, nv=0 )$u
#          locdes<-cbind( twoeventmat ) #, noisesvd )
          locdes<-matrix(twoeventmat[,1],ncol=1) 
          locdes<-finiteImpuleResponseDesignMatrix( denoisematmod1, n=bl, baseshift=baseshift )
          locdes<-cbind( locdes ) # , twoeventmat[,2] )
          mask2<-initdf2$mask
          mycca<-sparseDecom2( inmatrix=list(  ccamatin, locdes ),
                               initializationList =initdf$initlist,
                               initializationList2=initdf2$initlist, sparseness=sparseness,
                               nvecs=length(initdf$initlist), its=its, cthresh=c(bestvoxnum,0),
                               uselong=0, smooth=0, mycoption=mycoption, inmask=c(mask,mask2) )
          mytype<-typeof( mycca$eig1[[1]] )
          if ( mytype == "double" ) ccamat<-t(mycca$eig1) else ccamat<-imageListToMatrix( mycca$eig1, mask )
          mytype<-typeof( mycca$eig2[[1]] )
          if ( mytype == "double" ) ccamat2<-(mycca$eig2) else ccamat2<-t(imageListToMatrix( mycca$eig2, initdf2$mask ))
          if ( mean( data.matrix(ccamat2) ) < 0 ) ccamat2<-ccamat2*(-1)
#          bestv<-0
#          bestval<-0
#          for ( nv in 1:nvecs ) {
#              normval<-sum(abs(ccamat2[,nv]))
#              betablock<-sum((ccamat2[1:ncol(locdes),nv]))
#              if ( abs(betablock/normval) > bestval ) {
#                bestval <- abs(betablock/normval)
#                bestv<-nv
#              }
#          }
           bestv<-1
#          p1<-( data.matrix(submat) %*% t(ccamat) )[,bestv]
#          p2<-( data.matrix(locdes)  %*% data.matrix(ccamat2) )[,bestv]
#          eps<-1.e-4
#          ww<-( abs(p1) >= eps & abs(p2) >= eps )
#          locor<-cor.test( p1[ww], p2[ww] )$est
#          plot( p1[ww], p2[ww] )
#          print(paste("best",bestval,bestv,'cor',locor))
          hrf<-as.numeric(ccamat2[1:bl,bestv])
          meanhrfval<-mean(hrf)
          mxdf<-abs(max(hrf)-meanhrfval)
          mndf<-abs(min(hrf)-meanhrfval)
          if ( mndf > mxdf ) hrf<-hrf*(-1)
          myhrf<-hrf
#          myhrf<-data.frame(stl(ts(myhrf, frequency = 2),"per")$time.series)$trend
#          eventbetas[ct,]<-ccamat[bestv,]/max(abs(ccamat[bestv,]))
          if (  max(ccamat[1,]) <= 0  ) ccamat[1,]<-ccamat[1,]*(-1)
          y<-myhrf  # ccamatin %*% t(ccamat)
          loctime<-(row+baseshift):(row+baseshift+bl-1)
          if ( max(loctime) > nrow(ccamatin) ) loctime<-loctime[1:which(loctime==nrow(ccamatin))]
          temp<-as.numeric(cor( data.matrix(ccamatin[loctime,]), y[1:length(loctime)] ))
          eventbetas[ct,]<-temp
          eventhrfs[ct,]<-myhrf
          eventrows[ct]<-kkt[row]
          rownames(eventbetas)[ct]<-paste(designnames[col],eventrows[ct],sep='.')
          plot(ts( myhrf ) )
          print(paste(rownames(eventbetas)[ct],"ct:",ct,'...',ct/neventstot*100,"% Mx",max(temp),"Me",mean(temp) ))
          ct<-ct+1
      }
  } # row in ...
  } # cca-hrf
  runhrfs[rct,] <- runhrf
  }# runs
  rct<-rct+1
}
return(list(runhrfs=runhrfs ,eventbetas=eventbetas,eventhrfs=eventhrfs,
            eventrows=eventrows,hrf=colMeans(runhrfs)))
}
