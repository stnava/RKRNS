print("#####organize design#####")
if ( ! exists("wordinds" ) ) {
wordinds<-39:280
sentinds<-281:ncol(dmat)
dmatw<-dmat[,wordinds]
dmats<-dmat[,sentinds]
eventsw<-apply( dmatw, FUN=sum, MARGIN=1 )
eventswr<-apply( dmatw, FUN=sum, MARGIN=2 )
eventss<-apply( dmats, FUN=sum, MARGIN=1 )
if ( !is.na( removeSentLengthEffects ) ) removeSentLengthEffects<-as.numeric( dmat$nchar )
if ( !is.na(removeEventOverlap) ) removeEventOverlap<-eventss
eventssr<-apply( dmats, FUN=sum, MARGIN=2 )
dmatsblock<-dmats
blocksize<-responselength
nchar<-as.numeric( dmat$nchar )
for ( i in 1:blocksize ) {
    dmatsblock<-dmatsblock+ashift(dmats,c(i,0))
    nchar<-nchar+shift(nchar,1)
}
events1<-apply( dmats     , FUN=sum, MARGIN=1 )
eventtimes<-which( events1 > 0 )
eventss<-apply( dmatsblock, FUN=sum, MARGIN=1 )
sentences<-colnames(dmats)
words<-colnames(dmatw)
}
