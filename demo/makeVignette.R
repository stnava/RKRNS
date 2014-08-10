#!/usr/bin/env Rscript
library(RKRNS)
RKRNSsrcdir<-"./"
render( paste(RKRNSsrcdir,"vignettes/RKRNS.Rmd",sep='') ,'pdf_document')
