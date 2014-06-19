print("#################correlation of sentences based on bold#################")
# 0. create a vector denoting times of events w/ child vs not child at sentence level
eventtimes<-which( events1 > 0 )
sentnames<-colnames( dmats )
sentnames<-colnames(dmatw)
# 1. for each event, extract submatrix of bold, then vectorize that matrix
responselength<-20/tr # e.g. 15 seconds div by 0.5 tr => 30 volumes
eventinds<-matrix(rep(NA,length(eventtimes)*responselength),nrow=length(eventtimes))
eventinds[,1]<-eventtimes
for (  i in 2:ncol(eventinds) ) eventinds[,i]<-eventinds[,i-1]+1
eventinds[ eventinds > nrow(imatf) ]<-nrow(imatf)
featspace<-matrix(rep(NA, ncol(imatf)*length(eventtimes)*(responselength) ),nrow=length(eventtimes))
fspacenames<-rownames(featspace)
lablspace<-data.frame( dmats[eventtimes,] ) #JHU
dmatsnames<-colnames(dmats)
for ( i in 1:nrow( eventinds ) ) {
  kk<-(imatf[eventinds[i,],])
#  resSpec <- spec.mtm(imatf[eventinds[i,],10], k=10, nw=5.0, nFFT = "default",
#                    centreWithSlepians = TRUE, Ftest = TRUE,
#                    jackknife = FALSE, maxAdaptiveIterations = 100,
#                    plot = TRUE, na.action = na.fail) 
#  featspace[i,]<-c( apply(kk,FUN=median,MARGIN=2) )# /c( apply(kk,FUN=sd,MARGIN=2) ) # , apply(kk,FUN=sd,MARGIN=2) )
  featspace[i,]<-unlist(kk)
  fspacenames[i]<-dmatsnames[ which( dmats[eventtimes[i],  ] == 1  ) ]
#  if ( childVnochild[i] == 0 ) lablspace[i,]<-c(1,0) else lablspace[i,]<-c(0,1)
  if ( i %% 500 == 0 ) print(paste(i/length( eventtimes )*100,"%"))
}
# for ( i in 1:ncol(featspace) ) featspace[,i]<-rank(featspace[,i])
# colnames(featspace)<-colnames(imat)
rownames(featspace)<-(fspacenames)
agg<-aggregate( featspace , list(Sent=fspacenames), mean )
aggcor<-cor(t(data.matrix(agg[,2:ncol(agg)])))
colnames(aggcor)<-rownames(aggcor)<-dmatsnames
pdf("fspace_corr.pdf",width=256,height=256)
pheatmap(aggcor,fontsize = 80)
dev.off()
#pdf("fspace_corr.pdf",width=4096,height=4096)
#pheatmap(cor(featspace))
#dev.off()
