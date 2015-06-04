### PBR Network Analyses
### PB network modularity
##18Apr2015

library(bipartite)
library(methods)
source('global.R')
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

pbr <- list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x)
pbr <- lapply(pbr,function(x) x[,colnames(x)!='pb'])

### Observed Modularity
obs <- lapply(pbr,computeModules)

### Null modules
null.dir <- paste("npb",names(pbr),sep="/")
for (i in 1:length(pbr)){
    mkNull(pbr[[i]],5000,null.dir[i])
}
### Contribution to modularity

### Stats
