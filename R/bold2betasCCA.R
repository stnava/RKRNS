bold2betasCCA <- function( boldmatrix, designmatrix, blockNumb, bl=12, baseshift=5, mask=NA, sparseness=c(0,0), multievents=FALSE, polydegree=10, bestvoxnum=50, uselm=FALSE , nvecs=5, whichcols=NA,
  mycoption=1, its=10 )
{
par(mfrow=c(1,2))
rct<-1
if ( all(is.na( whichcols )) ) whichcols<-1:ncol(designmatrix)
allruns<-unique( blockNumb )
neventstot<-0
for ( runs in allruns ) 
  {
  kkt<-which( blockNumb == runs )
  neventstot<-neventstot+sum(designmatrix[kkt,whichcols])
  }
designnames<-colnames(designmatrix)[whichcols]
print(designnames)
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
  for ( row in 1:nrow(denoisedes) ) {
      if ( sum( denoisedes[row,whichcols] ) > 0 )
          {
          col<-which( denoisedes[row,whichcols] > 0 )
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
            locdes<-data.matrix(glmdf)
            mycca<-sparseDecom2( inmatrix=list(  data.matrix(submat) , locdes  ),
                                sparseness=sparseness, nvecs=nvecs, its=its, cthresh=c(bestvoxnum,0),
                                uselong=0, smooth=0, mycoption=mycoption, inmask=c(mask,NA) )
#            locdes[,1:bl]<-locdes[sample(1:nrow(submat)),1:bl] 
#            permcca<-sparseDecom2( inmatrix=list(  data.matrix(submat) , locdes ),
#                                sparseness=sparseness, nvecs=nvecs, its=5, cthresh=c(bestvoxnum,0),
#                                uselong=0, smooth=0, mycoption=0, inmask=c(mask,NA) )
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
            p1<-( data.matrix(submat) %*% t(ccamat) )[,bestv]
            p2<-( data.matrix(glmdf)  %*% data.matrix(mycca$eig2) )[,bestv]
            eps<-1.e-4
            ww<-( abs(p1) >= eps & abs(p2) >= eps )
            locor<-cor.test( p1[ww], p2[ww] )$est
            plot( p1[ww], p2[ww] )
            print(paste("best",bestval,bestv,'cor',locor))
            betablock<-mycca$eig2[1:bl,bestv]
            myhrf<-( (betablock ) )
            if ( mean(myhrf) < 0 ) myhrf<-myhrf*(-1)
            eventbetas[ct,]<-ccamat[bestv,]/max(abs(ccamat[bestv,]))
            temp<-(as.numeric(eventbetas[ct,]))
          }
          eventhrfs[ct,]<-myhrf
          eventrows[ct]<-kkt[row]
          rownames(eventbetas)[ct]<-paste(designnames[col],eventrows[ct],sep='.')
          plot(ts( myhrf ) )
          print(paste(rownames(eventbetas)[ct],"ct:",ct,'...',ct/neventstot*100,"% Mx",max(temp),"Me",mean(temp) ))
          ct<-ct+1
      }  
    }
  }
  rct<-rct+1
  }# runs
return(list(eventbetas=eventbetas,eventhrfs=eventhrfs,eventrows=eventrows))
}
