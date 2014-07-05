classificationNetwork <- function( nodesIn, realClass, predClass ,  whichviz="force", outfile='temp.html', zoom=F , opac=0.8, mycharge=-400 ) {
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
Value = "weight", NodeID = "nodename", Group = "nodeid", width = 1200, height = 1200, opacity = opac,zoom=zoom ,
   file = outfile, linkDistance = 200, fontsize = 12, charge = mycharge )
  else  d3ForceNetwork(Links = edgedata, Nodes = nodedata, Source = "from2", Target = "to2",zoom=zoom ,
Value = "weight", NodeID = "nodename", Group = "nodeid", width = 1200, height = 1200, opacity = opac, linkDistance = 200, fontsize = 12, charge = mycharge )
  }
  return( list( nodes=nodedata,  edges=edgedata, weightmatrix=weightmatrix ) )
}
