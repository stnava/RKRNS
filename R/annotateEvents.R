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
  sent_token_annotator <- Maxent_Sent_Token_Annotator()
  word_token_annotator <- Maxent_Word_Token_Annotator()
  pos_tag_annotator <- Maxent_POS_Tag_Annotator()
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
    a2 <- annotate(s, list(sent_token_annotator, word_token_annotator))
    a3 <- annotate(s, pos_tag_annotator, a2)
    a3w <- subset(a3, type == "word")
    tags <- sapply(a3w$features, `[[`, "POS")
#    grepverb<-grep("V",tags[locwordpos])
    grepverb<-grep("V",tags)
    if (length(grepverb)>0){
    everbs1[i]<-locwords[grepverb]
    rpl<-as.character(worddf$words) == everbs1[i]
    if ( sum(rpl) > 0 ) everbs1lab[i]<-worddf$wordids[ rpl ]
    }
#    grepnouns<-grep("NN",tags[locwordpos])
    grepnouns<-grep("NN",tags)
    enouns1[i]<-locwords[grepnouns[1]]
    if ( !is.na(enouns1[i]) )
      {
      rpl<-as.character(worddf$words) == enouns1[i]
      if ( sum(rpl) > 0 ) enouns1lab[i]<-worddf$wordids[ rpl ]
      nextnoun<-grepnouns[1]
      if ( length(grepnouns) > 1 ) nextnoun<-grepnouns[2]
      enouns2[i]<-locwords[nextnoun]
      rpl<-as.character(worddf$words) == enouns2[i]
      if ( sum(rpl) > 0 ) enouns2lab[i]<-worddf$wordids[ rpl ]
      }
  }
eventdata<-cbind( eventdata, enouns1=enouns1, enouns1lab=enouns1lab, everbs1=everbs1, everbs1lab=everbs1lab, enouns2=enouns2, enouns2lab=enouns2lab, sentlab=sentlab )

return( eventdata )
  
}

