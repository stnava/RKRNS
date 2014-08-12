#!/usr/bin/env Rscript
library(RKRNS)
RKRNSsrcdir<-"./"
render( paste(RKRNSsrcdir,"vignettes/RKRNS.Rmd",sep='') ,'pdf_document')
render( paste(RKRNSsrcdir,"vignettes/RKRNS.Rmd",sep='') ,'revealjs_presentation')
q()
render( "RKRNS/demo/RKRNS2.Rmd" ,'pdf_document')
