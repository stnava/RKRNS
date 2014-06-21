docompcor<-FALSE
eigsentbasislength<-25
library(ANTsR)
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
data("aal",package="ANTsR")
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
