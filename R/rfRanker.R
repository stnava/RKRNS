rfRanker <- function( mydfin , trainingRows, targetSuccessRate=0.15, usesparse=FALSE, verbose=FALSE ) {
  library(randomForest)
  library(RRF)
  if ( !usesparse ) mdl<-randomForest( lab ~., data=mydfin[trainingRows,] )
  if (  usesparse ) mdl<-RRF( lab ~., data=mydfin[trainingRows,] )
  err<-sum(mydfin[-trainingRows,]$lab==predict( mdl, newdata=mydfin[-trainingRows,]))/nrow(mydfin[-trainingRows,])
  truth<-mydfin[-trainingRows,]$lab
  pred<-predict( mdl, newdata=mydfin[-trainingRows,], norm.votes=TRUE)
  myvotes<-predict( mdl, newdata=mydfin[-trainingRows,] , type='vote', norm.votes=TRUE )
  xtab <- table(pred, truth)
  cm<-confusionMatrix(xtab)
  successthresh<-targetSuccessRate*length(levels(mydfin$lab)) # in top 15%
  succct<-0
  for ( j in 1:length(pred) ) {
    eventid<-as.numeric( names(pred)[j] )
    myvoteval<-myvotes[j,  mydfin$lab[eventid] == colnames(myvotes) ]
    ranker<-sum(myvotes[j,  ]>=myvoteval)
    if ( ranker <= successthresh ) {
      succct<-succct+1
      if ( verbose ) print(paste(mydfin$lab[eventid],myvoteval,ranker))
    }
  }
  if ( verbose ) print(paste( succct/length(pred) *100 , "%"))
  xtab <- table(pred, truth)
  cm<-confusionMatrix(xtab)
  dfout<-data.frame( errRate=err, successPercent=succct/length(pred) *100, confusionMatrix=cm )
  return( dfout )
}

