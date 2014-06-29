print("#########setup#########TODO: need to incorporate motion parameters!")
library(ANTsR)
istest<-FALSE
subject<-"111157"
tr<-as.numeric(0.5)
responselength<-8/tr # e.g. 15 seconds div by 0.5 tr => 30 volumes
labs<-as.numeric(1:90) # label numbers to use ... need to know which label set is at hand
labs<-as.numeric( c(13,79,89) ) # lang network
throwaway<-8
ncompcor<-5
compcorvarval<-0.95
filterlowfrequency<-0.1 # 0.05 if you multiply by TR
filterhighfrequency<-1.0 # 0.4 # because of expected bold response < 25secs, > 5 seconds
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
misclassnetwork <- function( nodesIn, realClass, predClass ,  whichviz="force", outfile='temp.html', zoom=F , opac=0.8, mycharge=-400 ) {
  if (nargs() == 0) {
    return(NA)
  }
  if ( length(realClass) != length( predClass ) ) return(NA)
  pckg <- try(require(d3Network))
  if (!pckg) {
    getPckg("d3Network")
  }
  library(d3Network)
  nodeoff<-nrow(nodesIn)
  nodedata<-nodesIn
  nodedata<-cbind( nodedata, zId=0:(nrow(nodedata)-1 ) )
  nodedata2<-nodedata
  nodedata2$nodename<-paste("P",nodedata$nodename,sep='_')
  nodedata2$zId<-nodedata2$zId+nrow(nodedata)
  nodedata<-rbind( nodedata, nodedata2 )
  row.names(nodedata)<-as.numeric(c(1:nrow(nodedata)))
  nedges<-rep(NA,length(realClass))
  nedges<-length(realClass)
  edgedata<-data.frame( from=rep(nodedata$nodename[1],nedges), to=rep(nodedata$nodename[1],nedges), weight=rep(0,nedges), from2=rep(NA,nedges) , to2=rep(NA,nedges) )
  weightmatrix<-matrix( rep(0,nrow(nodesIn)^2) , nrow=nrow(nodesIn) )
  rownames(weightmatrix)<-nodedata$nodename[1:nodeoff]
  colnames(weightmatrix)<-nodedata$nodename[(nodeoff+1):nrow(nodedata)]
  for ( i in 1:nedges )
    {
    whichTnode<-nodedata$nodeid == realClass[i]
    whichPnode<-nodedata$nodeid == predClass[i]
    whichPnode[1:nodeoff]<-FALSE
    whichTnode[(nodeoff+1):nrow(nodedata)]<-FALSE
    if (  any(whichTnode ) ) {
      kk<-(nodedata$nodename[ whichTnode ])
      edgedata$from[i]<-kk
      edgedata$from2[i]<-nodedata$zId[ whichTnode ]
      edgedata$to[i]<-kk
      edgedata$to2[i]<-nodedata$zId[ whichTnode ]
      matTentry<-as.numeric(row.names(nodedata)[whichTnode])
      if (   any(whichPnode ) ) {
        kk<-(nodedata$nodename[ whichPnode ])
        edgedata$to[i]<-kk
        edgedata$to2[i]<-nodedata$zId[ whichPnode ]
        matPentry<-as.numeric(row.names(nodedata)[whichPnode])
        weightmatrix[matTentry,matPentry-nodeoff]<-weightmatrix[matTentry,matPentry-nodeoff]+1
        }
      } else {
      edgedata$from[i]<-0
      edgedata$from2[i]<-0
      edgedata$to[i]<-0
      edgedata$to2[i]<-0
      }
    }
  for ( i in 1:nedges )
    {
    whichTnode<-nodedata$nodeid == realClass[i]
    whichPnode<-nodedata$nodeid == predClass[i]
    whichPnode[1:nodeoff]<-FALSE
    whichTnode[(nodeoff+1):nrow(nodedata)]<-FALSE
    if ( any(whichTnode) & any(whichPnode) )
      {
      matTentry<-as.numeric(row.names(nodedata)[whichTnode])
      matPentry<-as.numeric(row.names(nodedata)[whichPnode])
      edgedata$weight[i]<-weightmatrix[matTentry,matPentry-nodeoff]
      }
    }
    if ( whichviz == "Sankey" ) {
    d3Sankey(Links = edgedata, Nodes = nodedata, Source = "from2",
           Target = "to2", Value = "weight", NodeID = "nodeid",
           fontsize = 12, nodeWidth = 30, width = 700,file=outfile)
  } else if ( whichviz == "Simple" ) {
      d3SimpleNetwork( edgedata,  fontsize = 12, linkDistance = 200, opacity = opac,
                      file = outfile , width = 1200, height = 1200, charge = mycharge)
  } else {
  if (!is.na(outfile) ) d3ForceNetwork(Links = edgedata, Nodes = nodedata, Source = "from2", Target = "to2",
Value = "weight", NodeID = "nodename", Group = "nodeid", width = 1200, height = 1200, opacity = opac,
   file = outfile, linkDistance = 200, fontsize = 12, charge = mycharge )
  else  d3ForceNetwork(Links = edgedata, Nodes = nodedata, Source = "from2", Target = "to2",
Value = "weight", NodeID = "nodename", Group = "nodeid", width = 1200, height = 1200, opacity = opac, linkDistance = 200, fontsize = 12, charge = mycharge )
  }
  return( list( nodes=nodedata,  edges=edgedata, weightmatrix=weightmatrix ) )
}
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
library(seewave)
library(rCharts)
library(rmarkdown)
#library(rpud)

