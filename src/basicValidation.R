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


l1<-inds1[221:2000]
l2<-inds2[221:1000]
ccamats1<-list( featspace[ (lablspace$The.window.was.dusty.==1)[ l1 ]  , ], featspace[ (lablspace$The.window.was.dusty.==1)[ l2 ]  , ] )
fcca1<-sparseDecom2( inmatrix=ccamats1, nvecs=2, sparseness=c( -0.2, -0.2 ), its=33, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=1, ell1= 1 , perms=500, robust=1 )  # subaal
l1<-inds1[1000:2000]
l2<-inds2[221:1000]
ccamats2<-list( featspace[ (lablspace$The.coffee.was.hot.==1)[ l1 ]  , ], featspace[ (lablspace$The.coffee.was.hot.==1)[ l2 ]  , ] )
fcca2<-sparseDecom2( inmatrix=ccamats2, nvecs=2, sparseness=c( -0.2, -0.2 ), its=33, mycoption=1, inmask=c(NA,NA ), cthresh=c(0,10), uselong=1, ell1= 1 , perms=500, robust=1 )  # subaal

print( summary(glm( lablspace$The.window.was.dusty ~ featspace, family="binomial" ) ) )
library(randomForest)
library(e1071)
botr<-( scale( featspace[inds1,] )  )
bote<-( scale( featspace[inds2,] )  )
trdf<-data.frame( ch=factor(lablspace[inds1,colnames(lablspace)=="The.small.boy.feared.the.storm."]), bold=botr )
tedf<-data.frame( ch=factor(lablspace[inds2,colnames(lablspace)=="The.small.boy.feared.the.storm."]), bold=bote )
myrf<-randomForest( ch ~ . , data=trdf )
pred<-predict( myrf , tedf )
cor.test(pred,tedf$ch)
