print("########basic validation########")
nl<-nrow(  featspace )
# featspace<-scale(featspace,center=FALSE)
print(paste("NL rows:",nl))
inds1<-1:(nl/2)+5
inds2<-(max(inds1)+1):nl
meanfeat1<-apply(featspace[ inds1  , ], FUN=mean,MARGIN=2)
meanfeat2<-apply(featspace[ inds1  , ], FUN=mean,MARGIN=2)
ww1<-apply(featspace[ (lablspace$The.small.boy.feared.the.storm==1)[inds1]  , ], FUN=mean,MARGIN=2)
ww2<-apply(featspace[ (lablspace$The.small.boy.feared.the.storm==1)[inds2]  , ], FUN=mean,MARGIN=2)
an1<-apply(featspace[ (lablspace$The.red.plane.flew.through.the.cloud.==1)[inds1]  , ], FUN=mean,MARGIN=2)
an2<-apply(featspace[ (lablspace$The.red.plane.flew.through.the.cloud.==1)[inds2]  , ], FUN=mean,MARGIN=2)
ch1<-apply(featspace[ (lablspace$The.window.was.dusty.==1)[inds1]  , ], FUN=mean,MARGIN=2)
ch2<-apply(featspace[ (lablspace$The.window.was.dusty.==1)[inds2]  , ], FUN=mean,MARGIN=2)
dd1<-apply(featspace[ (lablspace$The.coffee.was.hot.==1)[inds1]  , ], FUN=mean,MARGIN=2)
dd2<-apply(featspace[ (lablspace$The.coffee.was.hot.==1)[inds2]  , ], FUN=mean,MARGIN=2)
meantfeat1<-meanfeat2<-1

sum(abs(ch2-ww1))
sum(abs(an2-ww1))
sum(abs(ww2-ww1))
sum(abs(dd2-ww1))
print("an1")
sum(abs(ch2-an1))
sum(abs(an2-an1))
sum(abs(ww2-an1))
sum(abs(dd2-an1))
print("dd1")
sum(abs(ch2-dd1))
sum(abs(an2-dd1))
sum(abs(ww2-dd1))
sum(abs(dd2-dd1))


docca<-T 
if ( docca == TRUE ) {
  longc<-1
  nperm<-0
  nv<-1; its<-25
  mysparse<-c( 0.1, 0.1 )
  myrob<-1
  redlist<-grep("red", fspacenames )
  l1<-1:(length(redlist)/2)
  l2<-((max(l1)+1):length(redlist))
  ccamats1<-list( featspace[ redlist[l1]  , ], featspace[ redlist[l2]  , ] )
#  ccamats1<-list( featspace[ redlist  , ], sentspace[ redlist , ] )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=longc, ell1= 1 , perms=nperm, robust=myrob )  # subaal
  mm<-matrix(fcca1$eig1[,1],nrow=responselength)
#  for ( i in 1:90 ) { print(i); plot(mm[,i],type='l'); if ( mean(abs(mm[,i])>0)) Sys.sleep(4) ; }

  grnlist<-( grep("artist", fspacenames ) )
  l1<-(1:(length(grnlist)/2))
  l2<-((max(l1)+1):(length(grnlist)))
  ccamats2<-list( featspace[ grnlist[l1]  , ], featspace[ grnlist[l2]  , ] )
  fcca2<-sparseDecom2( inmatrix=ccamats2, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=longc, ell1= 1 , perms=nperm, robust=myrob )  # subaal
#  mm<-matrix(fcca1$eig1[,1],nrow=responselength)
#  for ( i in 1:90 ) { print(i); plot(mm[,i],type='l'); if ( mean(abs(mm[,i])>0)) Sys.sleep(4) ; }

  
  ccamats3<-list( ccamats1[[1]] , ccamats2[[2]] )
  fcca3<-sparseDecom2( inmatrix=ccamats3, nvecs=nv, sparseness=mysparse, its=its, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=longc, ell1= 1 , perms=nperm, robust=myrob )  # subaal


  print( fcca1$ccasummary[1,] ) 
  print( fcca2$ccasummary[1,] ) 
  print( fcca3$ccasummary[1,] ) 
}

# print( summary(glm( lablspace$The.window.was.dusty ~ featspace, family="binomial" ) ) )


decode<-FALSE
if ( decode == TRUE ) {
  botr<-( scale( featspace[inds1,] )  )
  bote<-( scale( featspace[inds2,] )  )
  trdf<-data.frame( ch=factor(lablspace[inds1,colnames(lablspace)=="The.small.boy.feared.the.storm."]), bold=botr )
  tedf<-data.frame( ch=factor(lablspace[inds2,colnames(lablspace)=="The.small.boy.feared.the.storm."]), bold=bote )
  myrf<-randomForest( ch ~ . , data=trdf )
  pred<-predict( myrf , tedf )
  cor.test(pred,tedf$ch)
}



ccadecode<-FALSE
if ( ccadecode ) {
  nperm<-2
  nvecs<-2
  mysparse<-rep( 0.05, nvecs )
  myrob<-0
  matvec<-matrix(  lablspace$The.window.was.dusty.[  ] , ncol=1 ) 
  ccamats1<-list(  matvec , featspace[, ]  )
  fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=nvecs, sparseness=mysparse, its=15, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=0, ell1= 1 , perms=nperm, robust=myrob, z=0.1 )  # subaal
  mm<-matrix(fcca1$eig2[,1],nrow=30)
  for ( i in 1:90 ) { print(paste(colnames(imat)[i],mean(abs(mm[,i])))); plot(mm[,i],type='l'); if ( mean(abs(mm[,i])>0)) Sys.sleep(4) ; }
}
