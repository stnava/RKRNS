#####organize design#####
wordinds<-14:280
sentinds<-281:ncol(dmat)
dmatw<-dmat[,wordinds]
dmats<-dmat[,sentinds]
eventsw<-apply( dmatw, FUN=sum, MARGIN=1 )
eventswr<-apply( dmatw, FUN=sum, MARGIN=2 )
eventss<-apply( dmats, FUN=sum, MARGIN=1 )
nchar<-as.numeric( dmat$nchar )
eventssr<-apply( dmats, FUN=sum, MARGIN=2 )
tr<-0.5


