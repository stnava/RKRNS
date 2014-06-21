print("######### assemble overarching data-organization & design files #########")
whichblocks<-Sys.glob("design/*csv")
whichblocks2<-unlist( strsplit( whichblocks, "/" ) )[c(1:length(whichblocks))*2]
whichblocks<-gsub("_design.csv","",whichblocks2)
dsplits<-paste("design/",whichblocks,"_design.csv",sep='')
dfn<-paste("assembly/assembled_design_",labs[1],"_",labs[length(labs)],"test.csv",sep='')
afn<-paste("assembly/assembled_aal_",labs[1],"_",labs[length(labs)],"test.csv",sep='')
print(paste("assemble",dfn))
if ( ! exists("dmat") | ! exists("imat") ) {
  if ( file.exists( afn ) & file.exists( dfn ) ) {
    print(paste("read",afn))
    imat<-data.frame( read.csv( afn ) )
    dmat<-data.frame( read.csv( dfn ) )
  } else {
    dmat<-data.frame( read.csv( dsplits[1] ) )
    for ( i in 2:length(dsplits) ) {
      dmat2<-data.frame( read.csv( dsplits[i] ) )
      dmat<-rbind( dmat, dmat2 )
    }
  }
}
print("design assembly done")
usedesignrow<-rep(FALSE,nrow(dmat))
