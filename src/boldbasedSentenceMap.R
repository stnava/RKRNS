print("#################correlation of sentences based on bold#################")
# 0. create a vector denoting times of events w/ child vs not child at sentence level
eventtimes<-which( events1 > 0 )
sentnames<-colnames( dmats )
sentnames<-colnames(dmatw)
# 1. for each event, extract submatrix of bold, then vectorize that matrix
responselength<-12/tr # e.g. 15 seconds div by 0.5 tr => 30 volumes
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
sentspace<-matrix(rep(NA, length(eventtimes)*(nsentences) ),nrow=length(eventtimes))
for ( i in 1:nrow( eventinds ) )   sentspace[i,]<-sentsimilarity[  which(colnames(sentsimilarity) == fspacenames[i]  ), ]    
sentspace<-matrix(rep(NA, length(eventtimes)*ncol(eigsent) ),nrow=length(eventtimes))
for ( i in 1:nrow( eventinds ) )   sentspace[i,]<-eigsent[  which(colnames(sentsimilarity) == fspacenames[i]  ), ]    
# for ( i in 1:ncol(featspace) ) featspace[,i]<-rank(featspace[,i])
# colnames(featspace)<-colnames(imat)
rownames(featspace)<-(fspacenames)
agg<-aggregate( featspace , list(Sent=fspacenames), mean )
aggcor<-cor(t(data.matrix(agg[,2:ncol(agg)])))
colnames(aggcor)<-rownames(aggcor)<-dmatsnames
pdf("fspace_corr.pdf",width=32,height=32)
pheatmap(aggcor,fontsize = 10)
dev.off()
#pdf("fspace_corr.pdf",width=4096,height=4096)
#pheatmap(cor(featspace))
#dev.off()

stop()
# dictionaries of all bold responses 
eanat1<-sparseDecom( featspace[ grep("red",rownames(featspace)),]    , sparseness=-0.9, nvecs=10, its=5, mycoption=1 )
eanat2<-sparseDecom( featspace[ grep("artist",rownames(featspace)),] , sparseness=-0.99, nvecs=10, its=5, mycoption=1 )
max(abs(cor(sentspace[grep("red",rownames(featspace)),],eanat1$projections[,])))
max(abs(cor(sentspace[grep("artist",rownames(featspace)),],eanat2$projections[,])))
artcor<-cor(sentspace[grep("artist",rownames(featspace)),],eanat2$projections[,])
rownames(artcor)<-paste("ESent",1:nwv,sep='')
pheatmap(artcor)
eanatmat<- featspace[grep("artist",rownames(featspace)),] %*% as.matrix(eanat2$eig)
sentmat<-sentspace[grep("artist",rownames(featspace)),]
summary(lm(sentmat~eanatmat))
nperm<-50
mysparse<-c( 0.1 , -0.5 )
myrob<-0
simccamats<-list( eanatmat, sentmat )
simspacecca<-sparseDecom2( inmatrix=simccamats, nvecs=10, sparseness=mysparse, its=50, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), ell1=-11 , perms=nperm, robust=0 )  # subaal

stop()
par(mfrow=c(1,1))
for ( i in 1:10 ) 
    for ( j in 4:4 ) {
        mm1<-matrix(eanat1$eig[,j],nrow=responselength)
        plot(mm1[,1],type='l');
        }
for ( i in 1:10 ) 
    for ( j in 4:6 ) {
        mm2<-matrix(eanat2$eig[,j],nrow=responselength)
        plot(mm2[,i],type='l',col='red');
        }


# dictionaries of sentence specific bold responses 
# these go into cca?
stop()
##### quick cca on fspace sentspace
nperm<-0
mysparse<-c( 0.2 , 0.5 )
myrob<-0
simccamats<-list( featspace, sentspace )
simspacecca<-sparseDecom2( inmatrix=simccamats, nvecs=4, sparseness=mysparse, its=5, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), ell1= 1 , perms=nperm, robust=myrob )  # subaal
stop()
