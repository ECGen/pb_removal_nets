### PBR Network Analyses
### PB network modularity
##18Apr2015

rm(list=ls())
nitr <- 1000
mat <- as.numeric(commandArgs(TRUE)[1]) ### which matrix do you want to use?
## [1] "c08"       "x08"       "c09"       "x09"       "c.avg"     "x.avg"    
## [7] "c08.npb"   "x08.npb"   "c09.npb"   "x09.npb"   "c.avg.npb" "x.avg.npb"

print(mat)

library(bipartite)
library(methods)
## source('global.R')
## source('pbrDataLoader.R')

### Setup data
## pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
## pbr.avg <- list(c.avg=meanMat(pbr$c08,pbr$c09),x.avg=meanMat(pbr$x08,pbr$x09))
## pbr.avg <- lapply(pbr.avg,function(x) apply(x,2,floor))
## pbr <- c(pbr,pbr.avg)
## pbr. <- lapply(pbr,function(x) x[,colnames(x)!='pb'])
## names(pbr.) <- paste(names(pbr.),'.npb',sep='')
## pbr <- c(pbr,pbr.)
pbr <- c("c08","x08","c09","x09","c.avg","x.avg","c08.npb","x08.npb","c09.npb","x09.npb","c.avg.npb","x.avg.npb")
names(pbr) <- pbr
null.dir <- paste("../../data/npb",names(pbr),sep="/")

### Null modules
x <- lapply(paste(null.dir[[mat]],1:nitr,sep='/'),read.csv,header=FALSE)
mod <- numeric()
for (j in 1:length(x)){
    mod[j] <- slot(computeModules(x[[j]]),name='likelihood')
}
write.csv(mod,file=paste('../../results/nullMods',names(pbr)[mat],'.csv',sep=''),row.names=FALSE)
