###Modularity analysis for the pbr study
###MKLau 11Jul2014

source('../src/pbr_load_data.R')
source('../src/pbr_funcs.R')

###Manage data
g08 <- geno.08[[1]]
g09 <- geno.09[[1]]
n08c <- apply(pbr.08$c[,apply(sign(pbr.08$c),2,sum)>3],2,function(x) x/max(x))
n08x <- apply(pbr.08$x[,apply(sign(pbr.08$x),2,sum)>3],2,function(x) x/max(x))
n09c <- apply(pbr.09$c[,apply(sign(pbr.09$c),2,sum)>3],2,function(x) x/max(x))
n09x <- apply(pbr.09$x[,apply(sign(pbr.09$x),2,sum)>3],2,function(x) x/max(x))
###Test modularity
dput(list(sigMods(n09c),sigMods(n08x)),file='../results/sigMods09.rda')
###Calculate contribution to modularity
dput(list(cMods(n09c),cMods(n08x)),file='../results/cMods09.rda')
