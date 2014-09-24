eigenSentences <- function(  wordembed, functiontoapply = sum , normalize=F, eigsentbasislength = NA, sentencesIn=NA, eventdata=NA ) {
######## build the average eigenword for each sentence #########
if ( is.na(sentencesIn) )
  {
  data(sentences,package="RKRNS")
  sentencesIn<-sentences
  }
words<-wordembed$WhichWord
nwords<-length(words)
nsentencesIn<-length(sentencesIn$Sentence)
if ( is.na(eigsentbasislength) ) eigsentbasislength<-ncol( wordembed ) - 1
if ( eigsentbasislength > (ncol( wordembed ) - 1) ) eigsentbasislength<-ncol( wordembed ) - 1
sclfactor<-1
if ( typeof(functiontoapply) == "character" ) if ( functiontoapply  == "svd" ) sclfactor<-2
if ( typeof(functiontoapply) == "character" ) if ( functiontoapply  == "nvn" ) sclfactor<-3
data(eigen_sentences, package = "RKRNS")
eigsent<-matrix( rep( 0, (nsentencesIn) * eigsentbasislength * sclfactor ) ,  nrow=nsentencesIn )
rownames( eigsent )<-sentencesIn$Sentence
sentct<-0 # should = nsentencesIn
for ( i in 1:nsentencesIn )
  {
  locwords<-unlist(strsplit(as.character(sentencesIn$Sentence[i]),"\\."))
  locwords<-words[ words %in% locwords ]
  if ( length(locwords) > 0 )
    {
    sentct<-sentct+1
    sentmat<-matrix( rep( 0, length(locwords) * eigsentbasislength ) ,  nrow=length(locwords) )
    for ( j in 1:length(locwords) )
      {    
      whichword2<-wordembed$WhichWord == locwords[j]
      wordvec <- wordembed[whichword2,2:(2+eigsentbasislength-1)]
      wvec<-as.numeric( wordvec )
      sentmat[j,]<-wvec
      }
    sentmat<-antsrimpute( sentmat )
    if ( typeof(functiontoapply) == "character" )
      {
      if ( functiontoapply  == "svd" )
        {
        pj<-svd(sentmat)$u %*% sentmat
        if ( dim(pj)[1] == 1 ) pj<-rbind(pj,pj)
        eigsent[sentct,]<-as.numeric(pj[1:2,])
        }
      if ( functiontoapply  == "joao" )
        {
        whichesent<-which( eigen_sentences$Sentence == sentencesIn$Sentence[i])
        esent<-eigen_sentences[ whichesent, 4:(4+eigsentbasislength-1) ]
        eigsent[sentct,]<-as.numeric(esent)
        }
      if ( functiontoapply  == "nvn" )
        {
        ww<-which( eventdata$sentences == sentencesIn$Sentence[i] )[1]
        whichword2<-as.character(wordembed$WhichWord)== as.character(eventdata$enouns1[ww])
        if ( !is.na(whichword2) ) wordvec1 <- wordembed[whichword2,2:(2+eigsentbasislength-1)] else wordvec1<-rep(NA,eigsentbasislength)
        whichword2<-as.character(wordembed$WhichWord)== as.character(eventdata$everbs1[ww])
        if ( (sum(whichword2) > 0 ) & !is.na(whichword2) )
           wordvec2 <- wordembed[whichword2,2:(2+eigsentbasislength-1)] else wordvec2<-rep(NA,eigsentbasislength)
        whichword2<-as.character(wordembed$WhichWord)== as.character(eventdata$enouns2[ww])
        if ( !is.na(whichword2) ) wordvec3 <- wordembed[whichword2,2:(2+eigsentbasislength-1)] else wordvec3<-rep(NA,eigsentbasislength)
        eigsent[sentct,]<-c(as.numeric(wordvec1),as.numeric(wordvec2),as.numeric(wordvec3))
        }
      }
    else
      {
      eigsent[sentct,]<-apply(sentmat,FUN=functiontoapply,MARGIN=2)
      }
    }
  }
  eigsent<-antsrimpute(eigsent)
if ( normalize )
  {
  eigsentmag<-sqrt( rowSums(eigsent * eigsent) )
  eigsent<-eigsent/eigsentmag
  }
return( eigsent )
}
