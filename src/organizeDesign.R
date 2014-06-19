print("#####organize design#####")
wordinds<-39:280
sentinds<-281:ncol(dmat)
dmatw<-dmat[,wordinds]
dmats<-dmat[,sentinds]
eventsw<-apply( dmatw, FUN=sum, MARGIN=1 )
eventswr<-apply( dmatw, FUN=sum, MARGIN=2 )
eventss<-apply( dmats, FUN=sum, MARGIN=1 )
nchar<-as.numeric( dmat$nchar )
eventssr<-apply( dmats, FUN=sum, MARGIN=2 )
tr<-0.5
dmatsblock<-dmats
blocksize<-20
nchar<-as.numeric( dmat$nchar )
for ( i in 1:blocksize ) {
    dmatsblock<-dmatsblock+ashift(dmats,c(i,0))
    nchar<-nchar+shift(nchar,1)
}
events1<-apply( dmats     , FUN=sum, MARGIN=1 )
eventss<-apply( dmatsblock, FUN=sum, MARGIN=1 )


