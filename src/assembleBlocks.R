#######################assemble#######################
aalimg<-antsImageRead( "aal/111157_aal2.nii.gz" , 3 )
bmask<-antsImageRead( "ref/111157_mask.nii.gz" , 3 )
aalmask<-antsImageClone( aalimg )
aalmask[ aalmask > 1 ]<-1
subaal<-antsImageClone( aalimg )
subaal[ aalimg > 0 ]<-0
for ( lab in labs ) subaal[ aalimg == lab ]<-lab
isplits<-paste("aalsplit/",whichblocks,"_moco",lab,".mha",sep='')
isplits<-paste("nii/",whichblocks,".nii.gz",sep='')
isplits<-paste("moco/",whichblocks,"_moco.nii.gz",sep='')
print(paste("assemble",afn))
throwaway<-8
ncompcor<-8
if ( ! file.exists( afn ) | ! exists("imat")  ) {
  imat<-timeseries2matrix( antsImageRead(isplits[1],4), subaal )
  for( j in 1:throwaway ) imat[j,]<-apply( imat[(throwaway+1):nrow(imat),], FUN=mean, MARGIN=2 )
  if ( docompcor ) {
    imatfull<-timeseries2matrix( antsImageRead(isplits[1],4), bmask )
    mycompcorv<-compcor( imatfull, ncompcor, variance_extreme = 0.95, returnv=TRUE )
    highvarmat<-compcor( imatfull, ncompcor, variance_extreme = 0.95, returnhighvarmat=TRUE )
    mycompcor<- scale(highvarmat %*% mycompcorv)
    mycompcor<-compcor( imatfull, ncompcor, variance_extreme = 0.95 )
    imat<-residuals(lm(imat~0+mycompcor))
    }
  imat<-imat/mean(imat)
  usedesignrow[1:nrow(imat)]<-TRUE
  nisplits<-length(isplits)
  for ( i in 2:nisplits ) {
    img<-antsImageRead(isplits[i],4)
    if ( dim(img)[1] == dim(aalimg)[1] &  dim(img)[2] == dim(aalimg)[2] &  dim(img)[3] == dim(aalimg)[3] )
      {
      imat2<-timeseries2matrix( img, subaal )
      for( j in 1:throwaway ) imat2[j,]<-apply( imat2[(throwaway+1):nrow(imat2),], FUN=mean, MARGIN=2 )
      if ( docompcor ) {
         imatfull<-timeseries2matrix( img, bmask )
#         highvarmat<-compcor( imatfull, ncompcor, variance_extreme = 0.95, returnhighvarmat=TRUE )
#         mycompcor<- scale(highvarmat %*% mycompcorv)
         mycompcor<-compcor( imatfull, ncompcor, variance_extreme = 0.95 )
         imat2<-residuals(lm(imat2~0+mycompcor))
       }
      imat2<-imat2/mean(imat2)
      oldrow<-nrow(imat)
      imat<-rbind( imat, imat2)
      newrow<-nrow(imat)
      usedesignrow[(oldrow+1):newrow]<-TRUE
      } 
  print(paste(i/length(isplits)*100,"%"))
  }
  imat<-data.frame(imat)
  colnames(imat)<-aal$label_name[labs]
  write.csv(imat,afn,row.names=F)
}
imat<-data.frame(imat)
colnames(imat)<-aal$label_name[labs]
dmat<-dmat[usedesignrow,]
if ( nrow(dmat) != nrow(imat) ) print("CHECK DIMENSIONS MATCH!")
print("assembly done")
