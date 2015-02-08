#' Simple assembleDesigning function.
#' 
#' Assembles the study design from raw block data
#' 
#' 
#' @param blocksCSVlist a list of block-wise csv files detailing the design
#' @param datadir root data directory containing subject design, bold, etc
#' @param assembledDesignOutPrefix output file name for assembled design file -
#' will be read if it already exists (rather than reassembled)
#' @param assembledImageOutPrefix output file name for assembled bold file -
#' will be read if it already exists (rather than reassembled)
#' @return list is output containing organizational variables for design matrix
#' @author Avants BB, Phillips JS
#' @examples
#' 
#' \dontrun{
#' assembly<-assembleDesign( blocksCSVlist, datadir, dfn, afn )
#' }
#' 
assembleDesign <- function( blocksCSVlist , datadir, assembledDesignOutPrefix,  assembledImageOutPrefix )
{
# blocksCSVlist<-Sys.glob(paste(datadir,"design/*csv",sep=''))   # INPUT csv list
# dfn<-paste(datadir,"assembly/assembled_design_",labs[1],"_",labs[length(labs)],"test.csv",sep='') # INPUT out csv name
# afn<-paste(datadir,"assembly/assembled_aal_",labs[1],"_",labs[length(labs)],"test.mha",sep='') # INPUT out img name
testblock<-unlist( strsplit( blocksCSVlist[1], "/" ) )
whichblocks2<-unlist( strsplit( blocksCSVlist, "/" ) )[ seq( length(testblock), length(testblock)*length(blocksCSVlist), by=length(testblock) ) ]
whichblocks<-gsub("_design.csv","",whichblocks2)
dsplits<-paste(datadir,"design/",whichblocks,"_design.csv",sep='')
imat<-NA
if ( file.exists( afn ) & file.exists( dfn ) ) {
  imat<-as.matrix( as.array( antsImageRead( afn , 2 ) ) )
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
