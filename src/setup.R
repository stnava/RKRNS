print("#########setup#########TODO: need to incorporate motion parameters!")
library(ANTsR)
istest<-FALSE
subject<-"111157"
tr<-as.numeric(0.5)
responselength<-12/tr # e.g. 15 seconds div by 0.5 tr => 30 volumes
labs<-as.numeric(1:90) # label numbers to use ... need to know which label set is at hand
labs<-as.numeric( c(13,79,89) ) # lang network
throwaway<-8
ncompcor<-5
compcorvarval<-0.95
filterlowfrequency<-0.1
filterhighfrequency<-1 # because of expected bold response < 25secs, > 5 seconds
trendfrequency<-3
winsorval<-0.01
removeSentLengthEffects<-TRUE
removeEventOverlap<-FALSE
eigsentbasislength<-100
aalimg<-antsImageRead( paste("aal/",subject,"_aal2.nii.gz",sep='') , 3 )
bmask<-antsImageRead( paste("ref/",subject,"_mask.nii.gz",sep='') , 3 )
imagedir<-"moco/"
imagepostfix<-"_moco.nii.gz"
data("aal",package="ANTsR")
########################## that's the important stuff, above ##########################
mygamma <- function(x, a1 = 6.,   a2 = 12., b1 = 0.9, b2 = 0.9, cc = 0.35) {
    d1 <- a1 * b1
    d2 <- a2 * b2
    c1 <- (x/d1)^a1
    c2 <- cc * (x/d2)^a2
    res <- c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)
    res
  }
print("# define deconvolution bases")
basislength<-50
b1<-mygamma(c(1:basislength),  1,  5 ,0.9,0.9,0.05) 
b2<-mygamma(c(1:basislength),  5, 10 ,0.9,0.9,0.05) 
b3<-mygamma(c(1:basislength), 10, 15 ,0.9,0.9,0.05) 
b4<-mygamma(c(1:basislength), 15, 20 ,0.9,0.9,0.05) 
b5<-mygamma(c(1:basislength), 20, 25 ,0.9,0.9,0.05) 
b6<-mygamma(c(1:basislength), 25, 30 ,0.9,0.9,0.05)
basismat<-cbind( b1, b2, b3, b4, b5, b6 )
##########################################
library(ggplot2)
library(fpc)
library(waveslim)
library(FitARMA)
library(forecast)
require(multitaper)
library(pheatmap)
library(psych)
library(visreg)
library(randomForest)
library(e1071)
library(kernlab)
library(RRF)
library(vbmp)
library(tgp)
#library(rpud)

