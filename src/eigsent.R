print("######## build the average eigenword for each sentence #########")
wordembedr<-read.csv(paste(srcdir,"../data/reuters_words.csv",sep=''))
wordembedw<-read.csv(paste(srcdir,"../data/wiki_words.csv",sep=''))
wordembed<-wordembedw
words<-colnames(dmatw)
nwords<-length(words)
sentences<-colnames(dmats)
nsentences<-length(sentences)
eigsent<-matrix( rep( 0, (nsentences) * eigsentbasislength ) ,  nrow=nsentences )
rownames( eigsent )<-paste("sent",1:nsentences)
sentct<-0 # should = nsentences
for ( i in 1:nsentences ) {
    locwords<-unlist(strsplit(as.character(sentences[i]),"\\."))
    locwords<-words[ words %in% locwords ]
    if ( length(locwords) > 0 )
      {
      sentct<-sentct+1
      for ( j in 1:length(locwords) ) {    
        whichword2<-wordembed$WhichWord == locwords[j]
        wordvec <- wordembed[whichword2,2:(2+eigsentbasislength-1)]
        wvec<-as.numeric( wordvec )
#        if ( j == 1 )
            eigsent[sentct,]<-eigsent[sentct,]+wvec # else eigsent[sentct,]<-eigsent[sentct,]*wvec*(1/j)# /length(locwords)
      }
    }
}
eigsentmag<-sqrt( rowSums(eigsent * eigsent) )
eigsent<-eigsent/eigsentmag
rownames(eigsent)<-sentences
pdf("eigsentcor.pdf",width=32,height=32)
pheatmap(cor(t(eigsent)))
dev.off()
if ( sentct != nsentences ) stop("nsent")
sentsimilarity<-cor(t(eigsent))
diag(sentsimilarity)<-mean(sentsimilarity)
colnames(sentsimilarity)<-rownames(sentsimilarity)<-sentences
