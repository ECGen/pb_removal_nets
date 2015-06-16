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
obs.dir <- dir('../results/',full.names=TRUE)
obs.names <- sub('.rdata','',sub('../results//','',obs.dir[grepl('rdata',obs.dir)]))
x <- lapply(obs.dir[grepl('rdata',obs.dir)],dget)
obs <- unlist(lapply(x,slot,name='likelihood'))
names(obs) <- obs.names


null <- dir('../results/',full.names=TRUE)
null.names <- sub('.csv','',sub('../results//','',null[grepl('.csv',null)]))
null.names <- sub('nullMods','',null.names)
null <- lapply(null[grepl('.csv',null)],read.csv)
names(null) <- null.names

obs <- obs[names(obs) %in% names(null)]
null <- null[names(null) %in% names(obs)]
obs <- obs[order(names(obs))]
null <- null[order(names(null))]
if (all(names(null) == names(obs))){'Good!'}else{'Not so good...'}

### 
coa <- list()
for (i in 1:length(obs)){coa[[i]] <- coStats(obs[i],null[[i]][,1])}
coa <- do.call(rbind,coa)
rownames(coa) <- names(null)
colnames(coa) <- c('z','pval','obs','mean','sd')
coa
