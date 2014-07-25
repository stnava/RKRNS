glmDenoiseR <- function( boldmatrix, designmatrixIn , hrfbasislength=50, kfolds=4, whichbase = NA, pvalthresh=0.25, maxnoisepreds=12 )
{
designmatrix<-designmatrixIn[,colMeans(designmatrixIn)>0 ]
mygamma <- function(x, a1 = 6.,   a2 = 12., b1 = 0.9, b2 = 0.99, cc = 0.35) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- (x/d1)^a1
    c2 <- cc * (x/d2)^a2
    res <- c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)
    res
  }
# start with a reasonable default set of bases 
b1<-mygamma(c(1:hrfbasislength),  1,  5 )# ,0.9,0.9,0.05) 
b2<-mygamma(c(1:hrfbasislength),  5, 10 )#,0.9,0.9,0.05) 
b3<-mygamma(c(1:hrfbasislength), 10, 15 )#,0.9,0.9,0.05) 
b4<-mygamma(c(1:hrfbasislength), 15, 20 )#,0.9,0.9,0.05) 
b5<-mygamma(c(1:hrfbasislength), 20, 25 )#,0.9,0.9,0.05) 
b6<-mygamma(c(1:hrfbasislength), 25, 30 )#,0.9,0.9,0.05)
basismat<-cbind( b1, b2, b3, b4, b5, b6 )
if ( ! all(is.na(whichbase)) ) basismat<-basismat[,whichbase[ whichbase < ncol(basismat)]]
basismat<-as.matrix( basismat )
#################################################
####### convolutions on basis image matrix ######
nbase<-ncol(basismat) # number of basis functions
designmatrixext<-interleaveMatrixWithItself( designmatrix, nbase )
print("####### convolutions on basis image matrix ######")
for ( i in 1:ncol(designmatrixext) )
  {
  k<-( i %% nbase )
  if ( k == 0 ) k<-nbase
  designmatrixext[,i]<-conv( designmatrixext[,i]  , basismat[,k] )[1:nrow(designmatrixext)]
  }
# overall description of the method
# 1. regressors include - design + trends + noise-pool
# 2. need to find noise-pool by looking at high p-value voxels
# 3. cross-validate predictions using different numbers of noise regressors
# 4. select best n for predictors from noise pool
# 5. return the noise mask and the value for n
# 6. should this be done separtely on each session or run ?
# 7. note that we dont know which events are which ...

# make regressors
p<-stats::poly( 1:nrow(designmatrixext) ,degree=4)
polyplusdesign<-cbind(designmatrixext,p)
mylm<-lm( data.matrix(boldmatrix)  ~  polyplusdesign )
biglm1<-bigLMStats( mylm, 0.01 )
noisevox<-biglm1$pval.model > pvalthresh
noiseu<-svd( data.matrix(boldmatrix)[,noisevox], nv=0, nu=maxnoisepreds )$u
groups<-c()
grouplength<-round(nrow(boldmatrix)/10)-1
for ( k in 1:kfolds ) groups<-c(groups,rep(k,grouplength))
bestnoiselev<-rep(0,maxnoisepreds)
for ( i in 1:maxnoisepreds )
  {
  noisepredcorrs<-rep(0,kfolds)
  for ( k in 1:kfolds )
    {
    selector <- groups!=k
    mydf<-data.frame( polyplusdesign[selector,] , noiseu[selector,1:i] )
    mylm1<-lm( data.matrix( boldmatrix[selector,] )  ~  . , data=mydf )
    selector <- groups==k
    mydf<-data.frame( polyplusdesign[selector,] , noiseu[selector,1:i] )
    predmat<-predict(mylm1,newdata=mydf)
    realmat<-boldmatrix[selector,]
    abscor<-0
    for ( v in 1:ncol(predmat) ) abscor<-abscor+abs( cor( predmat[,v], realmat[,v] ) )
    abscor<-abscor/ncol(predmat)
    noisepredcorrs[k]<-abscor
    }
  print(paste("NoiseU:",i,"MeanCorr",mean(noisepredcorrs)))
  bestnoiselev[i]<-mean(noisepredcorrs)
  }
print(bestnoiselev)
mylm<-lm( data.matrix(boldmatrix)  ~  polyplusdesign + noiseu[selector,1:which.min(bestnoiselev)]  )
biglm2<-bigLMStats( mylm, 0.01 )
return( list( biglm1=biglm1, biglm2=biglm2 ) )
# v<-residuals(m)+coefficients(m)[1]
}
