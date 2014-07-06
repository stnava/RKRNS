assembleDesign <- function( blocksCSVlist , datadir, assembledDesignOutPrefix,  assembledImageOutPrefix )
{
# blocksCSVlist<-Sys.glob(paste(datadir,"design/*csv",sep=''))   # INPUT csv list
# dfn<-paste(datadir,"assembly/assembled_design_",labs[1],"_",labs[length(labs)],"test.csv",sep='') # INPUT out csv name
# afn<-paste(datadir,"assembly/assembled_aal_",labs[1],"_",labs[length(labs)],"test.mha",sep='') # INPUT out img name
whichblocks2<-unlist( strsplit( blocksCSVlist, "/" ) )[c(1:length(blocksCSVlist))*2]
whichblocks<-gsub("_design.csv","",whichblocks2)
dsplits<-paste(datadir,"design/",whichblocks,"_design.csv",sep='')
imat<-NA
if ( file.exists( afn ) & file.exists( dfn ) ) {
  imat<-as.matrix( antsImageRead( afn , 2 ) )
  dmat<-data.frame( read.csv( dfn ) )
} else {
  dmat<-data.frame( read.csv( dsplits[1] ) )
  for ( i in 2:length(dsplits) ) {
    dmat2<-data.frame( read.csv( dsplits[i] ) )
    dmat<-rbind( dmat, dmat2 )
  }
}
usedesignrow<-rep(FALSE,nrow(dmat))
# return imat and dmat
return( list( dmat=dmat, imat=imat, usedesignrow=usedesignrow ) )
}
