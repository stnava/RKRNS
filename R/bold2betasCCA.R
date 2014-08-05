bold2betasCCA <- function( boldmatrix, designmatrix, blockNumb, bl=12, baseshift=5, mask=NA, sparseness=c(0,0), multievents=FALSE, polydegree=10, bestvoxnum=50, uselm=FALSE , nvecs=5 )
{
par(mfrow=c(1,2))
fulleventbetas<-matrix(c(NA,NA),nrow=2)
rct<-1
allruns<-unique( blockNumb )
for ( runs in allruns[1] ) 
  {
  print(paste("ct",rct,"run",runs,":",rct/length(allruns)*100,"%"))
  kkt<-which( blockNumb == runs )
  denoisedes<-designmatrix[kkt,]
  submat<-boldmatrix[kkt,]
  nevents<-sum(denoisedes)
  eventbetas<-data.frame(matrix( rep(0,nevents*ncol(boldmatrix)), ncol=ncol(boldmatrix)))
  ct<-1
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,] ) > 0 )
          {
          col<-which( denoisedes[row,] > 0 )
          denoisematmod1<-denoisedes*0
          denoisematmod2<-denoisedes
          denoisematmod1[row,col]<-1
          denoisematmod2[row,col]<-0
          if ( multievents ) {
            rows<-which(denoisedes[,col]==1)
            rowsx<-which( (rows-row) > 0 )
            rows<-c( row, rows[rowsx[1]] )
            denoisematmod1[rows,col]<-1
            denoisematmod2[-rows,col]<-0
          }
          # now get the hrf estimating part 
          twoeventmat<-cbind( denoisematmod1[,col], rowSums(denoisematmod2))
          twoeventmat<-finiteImpuleResponseDesignMatrix( twoeventmat, n=bl, baseshift=baseshift )
          p<-stats::poly( 1:nrow(twoeventmat) ,degree=polydegree )
          glmdf<-data.frame( twoeventmat , z=p )
          if ( uselm ) {
            mylm<-lm(  data.matrix(submat) ~ . , data=glmdf )
            mylm<-bigLMStats( mylm , 0.001 )
            betablock<-mylm$beta.t[1:bl,]
            if ( sparseness[2] > 0 ) { betablock[betablock<0]<-0 }
            eventbetas[ct,]<-apply( abs(betablock) , FUN=sum, MARGIN=2)
            temp<-(as.numeric(eventbetas[ct,]))
            tempord<-sort(temp,decreasing=TRUE)
            bestvoxels<-which( temp > tempord[bestvoxnum]  )
            myhrf<-rowSums( (betablock[,bestvoxels] ) )
          } else { # cca
            mycca<-sparseDecom2( inmatrix=list(  data.matrix(submat) , data.matrix(glmdf)  ),
                                sparseness=sparseness, nvecs=nvecs, its=11, cthresh=c(bestvoxnum,0),
                                uselong=0, smooth=0, mycoption=1, inmask=c(mask,NA) )
            mytype<-typeof( mycca$eig1[[1]] )
            if ( mytype == "double" ) ccamat<-t(mycca$eig1) else ccamat<-imageListToMatrix( mycca$eig1, mask )
            if ( mean( data.matrix(mycca$eig2) ) < 0 ) mycca$eig2<-mycca$eig2*(-1)
            # find predictor most related to betablock1
            bestv<-0
            bestval<-0
            for ( nv in 1:nvecs ) {
              normval<-sum(abs(mycca$eig2[,nv]))
              betablock<-sum(abs(mycca$eig2[1:bl,nv]))
              if ( (betablock/normval) > bestval ) {
                bestval <- (betablock/normval)
                bestv<-nv
              }
            }
            print(paste("best",bestval,bestv))
            betablock<-mycca$eig2[1:bl,bestv]
            myhrf<-( (betablock ) )
            if ( mean(myhrf) < 0 ) myhrf<-myhrf*(-1)
            eventbetas[ct,]<-ccamat[bestv,]/max(abs(ccamat[bestv,]))
            plot( mycca$projections[,bestv], mycca$projections2[,bestv] )
            temp<-(as.numeric(eventbetas[ct,]))
          }
          plot(ts( myhrf ) )
#          rownames(eventbetas)[ct]<-paste(colnames(denoisedes)[col],".",rownames(boldmatrix)[row],sep='')
#          print(paste(rownames(eventbetas)[ct],"Mx",max(abs(mylm$beta.t[1,])),"Me",mean(abs(mylm$beta.t[1,])) ))
          print(paste(ct/nevents*100,"% Mx",max(temp),"Me",mean(temp) ))
          ct<-ct+1
      }  
    }
  if ( is.na( fulleventbetas[1,1] ) ) fulleventbetas<-eventbetas else fulleventbetas<-rbind(fulleventbetas,eventbetas)
  rct<-rct+1
  }# runs
return(fulleventbetas)
}
