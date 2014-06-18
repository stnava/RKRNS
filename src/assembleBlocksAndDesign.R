#######################assemble#######################
aalimg<-antsImageRead( "aal/111157_aal.nii.gz" , 3 )
labs<-as.numeric( c(13,79,89) )
labs<-as.numeric(1:90)
subaal<-antsImageClone( aalimg )
subaal[ aalimg > 0 ]<-0
for ( lab in labs ) subaal[ aalimg == lab ]<-lab
whichblocks<-Sys.glob("design/*csv")
whichblocks2<-unlist( strsplit( whichblocks, "/" ) )[c(1:length(whichblocks))*2]
whichblocks<-gsub("_design.csv","",whichblocks2)
isplits<-paste("aalsplit/",whichblocks,"_moco",lab,".mha",sep='')
isplits<-paste("nii/",whichblocks,".nii.gz",sep='')
isplits<-paste("moco/",whichblocks,"_moco.nii.gz",sep='')
dsplits<-paste("design/",whichblocks,"_design.csv",sep='')
# imat<-as.matrix(antsImageRead(isplits[1],2))
afn<-paste("assembly/assembled_aal_",labs[1],"_",labs[length(labs)],".csv",sep='')
dfn<-paste("assembly/assembled_design_",labs[1],"_",labs[length(labs)],".csv",sep='')
print(paste("assemble",afn))
if ( file.exists( afn ) & ! exists("imat")  ) {
  print(paste("read",afn))
  imat<-data.frame( read.csv( afn ) )
  dmat<-data.frame( read.csv( dfn ) )
}
throwaway<-8
if ( ! exists("imat") |  ! file.exists(afn)) {
imat<-timeseries2matrix( antsImageRead(isplits[1],4), subaal )
for( j in 1:throwaway ) imat[j,]<-apply( imat[(throwaway+1):nrow(imat),], FUN=mean, MARGIN=2 )
imat<-imat/mean(imat)
dmat<-data.frame( read.csv( dsplits[1] ) )
for ( i in 2:length(isplits) ) {
#  imat2<-as.matrix(antsImageRead(isplits[i],2))
  img<-antsImageRead(isplits[i],4)
  if ( dim(img)[1] == dim(aalimg)[1] &  dim(img)[2] == dim(aalimg)[2] &  dim(img)[3] == dim(aalimg)[3] )
    {
    imat2<-timeseries2matrix( img, subaal )
    for( j in 1:throwaway ) imat2[j,]<-apply( imat2[(throwaway+1):nrow(imat2),], FUN=mean, MARGIN=2 )
    imat2<-imat2/mean(imat2)
    imat<-rbind( imat, imat2)
    dmat2<-data.frame( read.csv( dsplits[i] ) )
    dmat<-rbind( dmat, dmat2 )
    }
  print(paste(i/length(isplits)*100,"%"))
  }
  imat<-data.frame(imat)
  colnames(imat)<-aal$label_name[labs]
  write.csv(imat,afn,row.names=F)
  write.csv(dmat,dfn,row.names=F)
}
imat<-data.frame(imat)
colnames(imat)<-aal$label_name[labs]
if ( nrow(dmat) != nrow(imat) ) print("CHECK DIMENSIONS MATCH!")
print("assembly done")


