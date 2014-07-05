eigenSentences <- function( functiontoapply = mean , normalize=T, eigsentbasislength = NA, option="wiki", sentencesIn=NA ) {
######## build the average eigenword for each sentence #########
if ( option == "wiki" )
  {
  data(wiki_words,package="RKRNS")
  wordembed<-wiki_words
  }
if ( option == "reuters" )
  {
  data(reuters_words,package="RKRNS")
  wordembed<-reuters_words
  }
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
eigsent<-matrix( rep( 0, (nsentencesIn) * eigsentbasislength ) ,  nrow=nsentencesIn )
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
    eigsent[sentct,]<-apply(sentmat,FUN=functiontoapply,MARGIN=2)
    }
  }
if ( normalize )
  {
  eigsentmag<-sqrt( rowSums(eigsent * eigsent) )
  eigsent<-eigsent/eigsentmag
  }
return( eigsent )
}
