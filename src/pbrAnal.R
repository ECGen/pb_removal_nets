### PBR Network Analyses
### PB network modularity
### 12 June 2015 

library(bipartite)
library(methods)
source('global.R')
source('pbrDataLoader.R')

### Setup data
pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
pbr.avg <- list(c.avg=meanMat(pbr$c08,pbr$c09),x.avg=meanMat(pbr$x08,pbr$x09))
pbr.avg <- lapply(pbr.avg,function(x) apply(x,2,floor))
pbr <- c(pbr,pbr.avg)
pbr. <- lapply(pbr,function(x) x[,colnames(x) != 'pb'])
pbr <- c(pbr,pbr.)
null.dir <- paste("../data/npb",names(pbr),sep="/")

### Null model based modularity analysis
obs <- dir('../results/',full.names=TRUE)
x <- lapply(obs[grepl('rdata',obs)],dget)
obs <- unlist(lapply(x,slot,name='likelihood'))
names(obs) <- c("c.avg.npb","c.avg","c08.npb","c08","c09.npb","c09","x.avg.npb","x.avg","x08.npb","x08","x09.npb","x09")

null <- dir('../results/',full.names=TRUE)
x <- lapply(null[grepl('.csv',null)],read.csv)
names(x) <- c("c.avg","c.avg.npb","c09","c09.npb","x.avg","x.avg.npb","x09","x09.npb")
null <- x

### 2009 c and x
obs.c09 <- obs[names(obs)=='c09']
obs.x09 <- obs[names(obs)=='x09']
null.c09 <- unlist(null[names(null)=='c09'])
null.x09 <- unlist(null[names(null)=='x09'])
zc09 <- (obs.c09 - mean(null.c09)) / sd(null.c09)
zx09 <- (obs.x09 - mean(null.x09)) / sd(null.x09)
pc09 <- length(null.c09[null.c09 >= obs.c09]) / length(null.c09)
px09 <- length(null.x09[null.x09 >= obs.x09]) / length(null.x09)

### average c and x
obs.cAVG <- obs[names(obs) == 'c.avg']
obs.xAVG <- obs[names(obs) == 'x.avg']
null.cAVG <- unlist(null[names(null)=='c.avg'])
null.xAVG <- unlist(null[names(null)=='x.avg'])
zcAVG <- (obs.cAVG - mean(null.cAVG)) / sd(null.cAVG)
zxAVG <- (obs.xAVG - mean(null.xAVG)) / sd(null.xAVG)
pcAVG <- length(null.cAVG[null.cAVG >= obs.cAVG]) / length(null.cAVG)
pxAVG <- length(null.xAVG[null.xAVG >= obs.xAVG]) / length(null.xAVG)


z <- c(zc09,zx09,zcAVG,zxAVG)
p <- c(pc09,px09,pcAVG,pxAVG)


