rkrnsRanker <- function( mydfin , trainingRows, targetSuccessRate=0.15, whichmodel='randomForest', verbose=FALSE ) {
  library(randomForest)
  library(RRF)
  if ( whichmodel == 'randomForest' ) mdl<-randomForest( lab ~., data=mydfin[trainingRows,] )
  if ( whichmodel == 'RRF' ) mdl<-RRF( lab ~., data=mydfin[trainingRows,] )
  if ( whichmodel == 'svm' ) mdl<-svm( lab ~., data=mydfin[trainingRows,],  kernel='linear',type='C-classification', probability = TRUE )
  err<-sum(mydfin[-trainingRows,]$lab==predict( mdl, newdata=mydfin[-trainingRows,]))/nrow(mydfin[-trainingRows,])
  truth<-mydfin[-trainingRows,]$lab
  if ( whichmodel == 'randomForest' | whichmodel == 'RRF' ) {
    pred<-predict( mdl, newdata=mydfin[-trainingRows,], norm.votes=TRUE)
    myvotes<-predict( mdl, newdata=mydfin[-trainingRows,] , type='vote', norm.votes=TRUE )
  }
  if ( whichmodel == 'svm' ) {
    pred<-predict( mdl, newdata=mydfin[-trainingRows,] )
    myvotes<-attr(predict( mdl, newdata=mydfin[-trainingRows,] , probability=TRUE ),"probabilities")
  }
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
  nodedf<-data.frame( nodename=mydfin$lab[trainingRows], nodeid=mydfin$lab[trainingRows] )
  ww<-  classificationNetwork( nodesIn=nodedf, truth, pred ,outfile='classificationNetwork.html', mycharge=-2066,zoom=T)
  dfout<-list( errRate=err, successPercent=succct/length(pred) *100, confusionMatrix=cm )
  return( dfout )
}

