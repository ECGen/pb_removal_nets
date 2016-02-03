### Contribution to Modularity
### MKLau 18Jun2015

restart <- 49

mat <- 'c.avg'
nits <- 100
null.dir <- '/N/u/mklau/Mason/mklau/pb_removal_nets/data/npb'

library(methods)
library(bipartite)
source('../src/global.R')
source('../src/pbrDataLoader.R')

### Load data
pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
pbr.avg <- list(c.avg=meanMat(pbr$c08,pbr$c09),x.avg=meanMat(pbr$x08,pbr$x09))
pbr.avg <- lapply(pbr.avg,function(x) apply(x,2,floor))
pbr.avg.npb <- lapply(pbr.avg,function(x) x[,-1])
names(pbr.avg.npb) <- paste(names(pbr.avg),'.npb',sep='')
pbr.avg <- c(pbr.avg,pbr.avg.npb)
x <- pbr.avg[names(pbr.avg) == mat][[1]]

### Load null data
null.dir <- paste(null.dir,mat,sep='/')
null <- lapply(1:nits,function(x,y) read.csv(paste(y,x,sep='/'),header=FALSE),y=null.dir)

print('Species Contribution')

start <- Sys.time()
print(start)

x <- t(x)

for (i in (restart):nrow(x)){
    out <- numeric()
    for (j in 1:length(null)){
        y <- x
        y[i,] <- as.numeric(t(null[[j]])[i,])
        mod <- computeModules(y)
        out[j] <- slot(mod,'likelihood')
    }
    write.csv(out,paste('results/',rownames(x)[i],sep=''))
    print(Sys.time())
    print(abs(Sys.time()-start))
}

