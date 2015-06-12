### PBR Network Analyses
### PB network modularity
##18Apr2015

n.null <- 5000
mknull <- FALSE
getobs <- FALSE
getnull <- TRUE

library(bipartite)
library(methods)
source('global.R')
source('pbrDataLoader.R')

##1. Compute observed modularity
##2. Generate null set
##3. Calculate modularity for null set
##4. Calculate contribution to modularity for each tree
##5. Caclulate contribution to modularity for each tree
##6. Repeat for each year
##7. Repeat using the average across years

### Setup data
pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
pbr.avg <- list(c.avg=meanMat(pbr$c08,pbr$c09),x.avg=meanMat(pbr$x08,pbr$x09))
pbr.avg <- lapply(pbr.avg,function(x) apply(x,2,floor))
pbr <- c(pbr,pbr.avg)
pbr. <- lapply(pbr,function(x) x[,colnames(x)!='pb'])
names(pbr.) <- paste(names(pbr.),'.npb',sep='')
pbr <- c(pbr,pbr.)
null.dir <- paste("../data/npb",names(pbr),sep="/")

### Null communities
if (mknull){
    sapply(null.dir,dir.create,recursive=TRUE)

    for (i in 1:length(pbr)){
        mkNull(pbr[[i]],n.null,null.dir[i])
    }
}

### Observed Modularity
if (getobs){
    obs <- lapply(pbr,computeModules)
    for (i in 1:length(obs)){
        dput(obs[[i]],file=paste(null.dir[[i]],'/',names(obs)[i],'.rdata',sep=''))
    }
}

### Null modularity
### Stupidly parallelized

### Contribution to modularity



### Stats
