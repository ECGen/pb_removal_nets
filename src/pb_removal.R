###Analysis of Art's Pemphigus removal data
###Current file started 18 Oct 2013

###Import data
pbr.08 <- read.csv('../data/keith_pb_removal_2008.csv')
pbr.09 <- read.csv('../data/keith_pb_removal_2009.csv')
pbr.08[,2] <- as.character(pbr.08[,2])
pbr.09[,2] <- as.character(pbr.09[,2])
                                        #fix typos
pbr.08[pbr.08[,2]=="SSG1-1c",2] <- "SSG1-1 c" 
pbr.08[pbr.08[,2]=="S7-3x",2] <- "S7-3 x"
pbr.08[pbr.08[,2]=="S7-3c",2] <- "S7-3 c"
pbr.08[pbr.08[,2]=="S4-9c",2] <- "S4-9 c"
pbr.08[pbr.08[,2]=="S3-21 x ",2] <- "S3-21 x"
pbr.08[pbr.08[,2]=="S3-21 c ",2] <- "S3-21 c"
pbr.09[pbr.09[,2]=="SSG1-1c",2] <- "SSG1-1 c" 
pbr.09[pbr.09[,2]=="S7-3x",2] <- "S7-3 x"
pbr.09[pbr.09[,2]=="S7-3c",2] <- "S7-3 c"
pbr.09[pbr.09[,2]=="S4-9c",2] <- "S4-9 c"
pbr.09[pbr.09[,2]=="S3-21 x ",2] <- "S3-21 x"
pbr.09[pbr.09[,2]=="S3-21 c ",2] <- "S3-21 c"
                                        #
pbr.08[is.na(pbr.08)] <- 0
pbr.09[is.na(pbr.09)] <- 0
geno.08 <- pbr.08[,1]
geno.09 <- pbr.09[,1]
tree.08 <- as.character(sapply(as.character(pbr.08[,2]),function(x) strsplit(x,split=' ')[[1]][1]))
tree.09 <- as.character(sapply(as.character(pbr.09[,2]),function(x) strsplit(x,split=' ')[[1]][1]))
trt.08 <- as.character(sapply(as.character(pbr.08[,2]),function(x) strsplit(x,split=' ')[[1]][2]))
trt.09 <- as.character(sapply(as.character(pbr.09[,2]),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #remove env data
pbr.08 <- pbr.08[,-1:-2]
pbr.09 <- pbr.09[,-1:-2]
                                        #separate by treatment
pbr.08 <- split(pbr.08,trt.08)
pbr.09 <- split(pbr.09,trt.09)
                                        #make binary
as.binary <- function(x){x[x!=0] <- 1;return(x)}
co.08 <- lapply(pbr.08,as.binary)
co.09 <- lapply(pbr.09,as.binary)

###Networks
library(sna)
source('~/projects/dissertation/projects/lichen_coo/src/seenetR.R')

my.gplot <- function(scn){
    e.col <- sign(scn)
    e.col[e.col==1] <- 'grey'
    e.col[e.col==-1] <- 'red'
    coord <- gplot(abs(scn),displaylabels=TRUE,gmode='graph',pad=1.5,
                   edge.col=e.col,edge.lwd=log(abs(scn)),vertex.col='lightblue',
                   mode='circle')
}

net.08 <- lapply(co.08,dep.net)
net.09 <- lapply(co.09,dep.net)
                                        #reduced nets
                                        #reduce speceis names to numbers
rownames(net.08[[1]]) <- colnames(net.08[[1]]) <- as.character(1:ncol(net.08[[1]]))
rownames(net.08[[2]]) <- colnames(net.08[[2]]) <- as.character(1:ncol(net.08[[2]]))
rownames(net.09[[1]]) <- colnames(net.09[[1]]) <- as.character(1:ncol(net.09[[1]]))
rownames(net.09[[2]]) <- colnames(net.09[[2]]) <- as.character(1:ncol(net.09[[2]]))
                                        #gplots
                                        #5.1 4.1 4.1 2.1
                                        #0.8466 0.6806 0.6806 0.3486
## par(mfrow=c(2,2),mar=c(0.1,0.1,0.1,2.1),mai=c(0.8466,0.6806,0.6806,0.3486)*0.2)
## my.gplot(net.08[[1]])
## mtext('2008',2,font=2)
## title('Control')
## my.gplot(net.08[[2]])
## title('Removed')
## my.gplot(net.09[[1]])
## mtext('2009',2,font=2)
## my.gplot(net.09[[2]])
                                        #edge weight distribution
## par(mfrow=c(2,2))
## hist(net.08[[1]])
## hist(net.08[[2]])
## hist(net.09[[1]])
## hist(net.09[[2]])
                                        #qap test
g.08 <- array(0,dim=c(nrow(net.08[[1]]),ncol(net.08[[1]]),2))
g.08[,,1] <- net.08[[1]]
g.08[,,2] <- net.08[[2]]
g.09 <- array(0,dim=c(nrow(net.09[[1]]),ncol(net.09[[1]]),2))
g.09[,,1] <- net.09[[1]]
g.09[,,2] <- net.09[[2]]
                                        #test
qap.08 <- qaptest(g.08,gcor,g1=1,g2=2)
dput(qap.08,file='../results/qap08.Rdata')
qap.09 <- qaptest(g.09,gcor,g1=1,g2=2)
dput(qap.09,file='../results/qap09.Rdata')
