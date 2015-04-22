### PBR Network Analyses
### PB network modularity
##18Apr2015

library(bipartite)
library(methods)
source('pbrDataLoader.R')

### compute modules without 

##1. Compute observed modularity
##2. Generate null set
##3. Calculate modularity for null set
##4. Calculate contribution to modularity for each tree
##5. Cacluulate contribution to modularity for each tree
##6. Repeat for each year
##7. Repeat using the average across years

###Remove PB

library(bipartite)

pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
pbr <- lapply(pbr,function(x) x[,colnames(x)!='pb'])

### Observed Modularity
obs <- lapply(pbr,computeModularity)
out <- write.table(unlist(lapply(obs,slot,name='likelihood')),file='../results/pbrObsMods.txt',col.names=FALSE,row.names=FALSE)

### Null Models
n.null <- 5000
dir.create('./tmp')
for (l in 1:length(pbr)){
    for (i in 1:n.null){
        write.table(nullmodel(pbr[[l]],method='r2d',N=1)[[1]],
                    file=paste('./tmp/pbrNul',names(pbr)[l],i,sep=''),
                    row.names=FALSE,col.names=FALSE,sep=',')
    }
}

### Null modularity
out <- numeric(length=length(dir('./tmp')))
names(out) <- dir('./tmp')
for (i in 1:length(dir(tmp))){
    x <- read.table(dir(tmp)[i])
    mod <- computeModules(x)
    out[i] <- slot(x,'likelihood')
}

write.table(out,file='../results/pbrNulMods.txt',col.names=FALSE)

### Calculate z values
null <- out
trt <- substr(names(null),8,10)
null <- split(null,trt)
obs <- read.table(file='../results/pbrObsMods.txt')
z <- numeric(length(null))
p <- numeric(length(null))
for (i in 1:length(null)){
    z[i] <- (obs[i] - mean(null[[i]])) / sd(null[[i]])
    p[i] <- length(null[[i]][null[[i]] >= obs[i]]) / length(null[[i]])
}

out <- rbind(obs,z,p)
rownames(out) <- c('obs','z','p')
colnames(out) <- names(null)
write.table(out,file='../results/pbrModOutput.csv',col.names=TRUE,row.names=TRUE,sep=',')

### Calculate contribution to modularity (trees)

### Calculate contribution to modularity (arthropods)


