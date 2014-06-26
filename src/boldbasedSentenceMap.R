print("#################correlation of sentences based on bold#################")
# 0. create a vector denoting times of events w/ child vs not child at sentence level
eventtimes<-which( events1 > 0 )
sentnames<-colnames( dmats )
sentnames<-colnames(dmatw)
# 1. for each event, extract submatrix of bold, then vectorize that matrix
fspacenames<-rep("", length(eventtimes) )
dmatsnames<-colnames(dmats)
for ( i in 1:nrow( eventinds ) ) {
#  resSpec <- spec.mtm(imatf[eventinds[i,],10], k=10, nw=5.0, nFFT = "default",
#                    centreWithSlepians = TRUE, Ftest = TRUE,
#                    jackknife = FALSE, maxAdaptiveIterations = 100,
#                    plot = TRUE, na.action = na.fail) 
#  featspace[i,]<-c( apply(kk,FUN=median,MARGIN=2) )# /c( apply(kk,FUN=sd,MARGIN=2) ) # , apply(kk,FUN=sd,MARGIN=2) )
  featspace[i,]<-unlist(kk)
  fspacenames[i]<-dmatsnames[ which( dmats[eventtimes[i],  ] == 1  ) ]
}
featspaceOrg<-timeserieswindow2matrix( data.matrix( imatf ), subaal, eventtimes, responselength, 3, c(1,1,1,0.25) )
featspace<-featspaceOrg$eventmatrix
mask4d<-featspaceOrg$mask4d
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

if ( FALSE ) {
# dictionaries of all bold responses 
eanat1<-sparseDecom( featspace[ grep("red",rownames(featspace)),]    , sparseness=-0.9, nvecs=10, its=5, mycoption=1 )
eanat2<-sparseDecom( featspace[ grep("artist",rownames(featspace)),] , sparseness=-0.99, nvecs=10, its=5, mycoption=1 )
max(abs(cor(sentspace[grep("red",rownames(featspace)),],eanat1$projections[,])))
max(abs(cor(sentspace[grep("artist",rownames(featspace)),],eanat2$projections[,])))
artcor<-cor(sentspace[grep("artist",rownames(featspace)),],eanat2$projections[,])
rownames(artcor)<-paste("ESent",1:eigsentbasislength,sep='')
pheatmap(artcor)
eanatmat<- featspace[grep("artist",rownames(featspace)),] %*% as.matrix(eanat2$eig)
sentmat<-sentspace[grep("artist",rownames(featspace)),]
summary(lm(sentmat~eanatmat))
nperm<-50
mysparse<-c( 0.1 , -0.5 )
myrob<-0
simccamats<-list( eanatmat, sentmat )
simspacecca<-sparseDecom2( inmatrix=simccamats, nvecs=10, sparseness=mysparse, its=50, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), ell1=-11 , perms=nperm, robust=0 )  # subaal

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
##### quick cca on fspace sentspace
nperm<-0
mysparse<-c( 0.2 , 0.5 )
myrob<-0
simccamats<-list( featspace, sentspace )
simspacecca<-sparseDecom2( inmatrix=simccamats, nvecs=4, sparseness=mysparse, its=5, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), ell1= 1 , perms=nperm, robust=myrob )  # subaal
}
