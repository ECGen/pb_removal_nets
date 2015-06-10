### PBR Network Analyses
### PB network modularity
##18Apr2015

mknull <- TRUE
getobs <- TRUE
getnull <- TRUE
n.null <- 5000

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

### Null communities
if (mknull){
    null.dir <- paste("../data/npb",names(pbr),sep="/")
    sapply(null.dir,dir.create,recursive=TRUE)

    for (i in 1:length(pbr)){
        mkNull(pbr[[i]],n.null,null.dir[i])
    }
}

### Observed Modularity
if (getobs){
    obs <- lapply(pbr,computeModules)
    for (i in 1:length(obs)){
        out <- computeModules(pbr[[i]])
        dput(out,file=paste(null.dir,'/',names(obs)[[i]],'.rdata',sep=''))
    }
}

### Null modularity
if (getnull){
    out <- list()
    mod <- numeric()
    for (i in 1:length(null.dir)){
        for (j in 1:length(dir(null.dir[[i]]))){
            x <- read.csv(paste(null.dir[[i]],j,sep='/'))
            mod[j] <- slot(computeModules(x),name='likelihood')
        }
        out[[i]] <- mod
    }
    out <- do.call(rbind,out)
    write.csv(out,file='../results/nullMods.csv',row.names=FALSE)
}

### Contribution to modularity



### Stats
