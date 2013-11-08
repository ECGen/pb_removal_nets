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
                                        #remove species occurring less than 10 times
pbr.08 <- pbr.08[,apply(pbr.08,2,sum)>=5]
pbr.09 <- pbr.09[,apply(pbr.09,2,sum)>=5]
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
source('~/projects/dissertation/projects/art_coo/src/helper_func.R')

print('Modeling networks')
net.08 <- lapply(co.08,dep.net)
net.09 <- lapply(co.09,dep.net)
                                        #reduced nets
                                        #reduce speceis names to numbers
rownames(net.08[[1]]) <- colnames(net.08[[1]]) <- as.character(1:ncol(net.08[[1]]))
rownames(net.08[[2]]) <- colnames(net.08[[2]]) <- as.character(1:ncol(net.08[[2]]))
rownames(net.09[[1]]) <- colnames(net.09[[1]]) <- as.character(1:ncol(net.09[[1]]))
rownames(net.09[[2]]) <- colnames(net.09[[2]]) <- as.character(1:ncol(net.09[[2]]))
                                        #correlations
edge.08 <- list(net.08[[1]][net.08[[1]]!=0|net.08[[2]]!=0],net.08[[2]][net.08[[1]]!=0|net.08[[2]]!=0])
names(edge.08) <- names(net.08)
plot(edge.08[[2]]~edge.08[[1]])
edge.09 <- list(net.09[[1]][net.09[[1]]!=0|net.09[[2]]!=0],net.09[[2]][net.09[[1]]!=0|net.09[[2]]!=0])
names(edge.09) <- names(net.09)
plot(edge.09[[2]]~edge.09[[1]])
                                        #gplots
par(mfrow=c(2,2))
coord.08=mgp2(net.08[[1]],v.cex=(apply(pbr.08[[1]],2,sum)/max(apply(pbr.08[[1]],2,sum))),scalar=1)
mgp2(net.08[[2]],v.cex=(apply(pbr.08[[2]],2,sum)/max(apply(pbr.08[[2]],2,sum))),scalar=1,my.coord=coord.08)
coord.09=mgp2(net.09[[1]],v.cex=(apply(pbr.09[[1]],2,sum)/max(apply(pbr.09[[1]],2,sum))),scalar=1)
mgp2(net.09[[2]],v.cex=(apply(pbr.09[[2]],2,sum)/max(apply(pbr.09[[2]],2,sum))),scalar=1,my.coord=coord.09)
                                        #qap test
g.08 <- array(0,dim=c(nrow(net.08[[1]]),ncol(net.08[[1]]),2))
g.08[,,1] <- net.08[[1]]
g.08[,,2] <- net.08[[2]]
g.09 <- array(0,dim=c(nrow(net.09[[1]]),ncol(net.09[[1]]),2))
g.09[,,1] <- net.09[[1]]
g.09[,,2] <- net.09[[2]]
                                        #test
print('running qap 2008')
qap.08 <- qaptest(g.08,gcor,g1=1,g2=2)
dput(qap.08,file='../results/qap08.Rdata')
print('running qap 2009')
qap.09 <- qaptest(g.09,gcor,g1=1,g2=2)
dput(qap.09,file='../results/qap09.Rdata')
print('gcor matrices')
dput(gcor(g.08),file='../data/gcor08.Rdata')
dput(gcor(g.09),file='../data/gcor09.Rdata')
print('Done!')

library(sna)
qap.08 <- dget('../results/qap08.Rdata')
qap.09 <- dget('../results/qap09.Rdata')
summary(qap.08)
summary(qap.09)
