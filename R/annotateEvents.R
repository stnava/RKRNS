#' Will match eigensentences to events and given unique labels to words,
#' sentences etc.
#'
#' Returns a dataframe for every event with annotated event times, word
#' categories (V, N,... ), sentence descriptors (labels, other variables) etc.
#'
#'
#' @param sentencesIn the list of possible sentences
#' @param words the list of possible words
#' @param eventtimes the index to the temporal volumes at which events occur
#' @param eventnames a name for every event so of equal length to eventtimes
#' @return returns a data frame
#' @author Avants BB
#' @examples
#'
#'   data(reuters_words,package="RKRNS")
#'   data(sentences,package="RKRNS")
#'   eventtimes<-c(10,25,99)
#'   eventnames<-list(sentences$Sentence[1],sentences$Sentence[55],sentences$Sentence[19])
#'   eventdf<-annotateEvents( sentences$Sentence, reuters_words$WhichWord, eventtimes, unlist(eventnames) )
#'
#' @export annotateEvents
annotateEvents <- function( sentencesIn, words, eventtimes, eventnames )
{
  sentenceids<-1:length(unique(sentencesIn))
  sentencedf<-data.frame(  sentences=sentencesIn, sentenceids=sentenceids )
  wordids<-1:length(unique(words))
  worddf<-data.frame( words=words, wordids=wordids )
  eventdata<-data.frame(eventtimes=eventtimes,sentences=eventnames)
  enouns1<-rep(NA,nrow(eventdata))
  enouns1lab<-rep(NA,nrow(eventdata))
  everbs1<-rep(NA,nrow(eventdata))
  everbs1lab<-rep(NA,nrow(eventdata))
  enouns2<-rep(NA,nrow(eventdata))
  enouns2lab<-rep(NA,nrow(eventdata))
  sentlab<-rep(NA,nrow(eventdata))
#  sent_token_annotator <- Maxent_Sent_Token_Annotator()
#  word_token_annotator <- Maxent_Word_Token_Annotator()
#  pos_tag_annotator <- Maxent_POS_Tag_Annotator()
  for ( i in 1:nrow(eventdata) ) {
    s<-eventdata$sentences[i]
    s<-gsub('[.]',' ',s)
    realwords<-unlist( strsplit(s,' ') )
    rpl<-as.character(sentencedf$sentences) == eventdata$sentences[i]
    if ( sum(rpl) > 0 ) sentlab[i]<-sentencedf$sentenceids[ rpl ]
    if ( i == 1 ) print("THIS IS A HACK THAT SHOULD BE FIXED")
    wordsmod<-paste(words,"ed",sep='')
    wordsmod2<-paste(words,"d",sep='')
    locwords<-realwords[ ( realwords  %in% words) |
                         ( realwords  %in% wordsmod )|
                         ( realwords  %in% wordsmod2 ) ]

    # find position of locwords in realwords
    locwordpos<-rep(0,length(locwords))
    ct<-1
    for ( k in locwords ) {
      pos<-1
      for ( j in realwords ) {
        if( length( grep(k,j) ) > 0 ) locwordpos[ct]<-pos;
        pos<-pos+1
        }
      ct<-ct+1
    }
    s<-paste(locwords,collapse=' ')
  }
eventdata<-cbind( eventdata, sentlab=sentlab )

return( eventdata )

}
