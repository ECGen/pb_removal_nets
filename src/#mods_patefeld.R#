###Testing for modularity in Art's removal data using
###the Patefield method Bluthgen 2008
###MKLau 29Jul2014

source('~/projects/packages/ComGenR/R/mod_test.R')
source('../src/pbr_load_data.R')
n09c <- pbr.09$c;n09x <- pbr.09$x #re-name
n08c <- pbr.08$c;n08x <- pbr.08$x #re-name

nits <- 100

###2009
start <- Sys.time()
##control with pb
mt09 <- list(testMods(n09c,nits=nits))
##control w/o pb
mt09[[2]] <- list(testMods(n09c[,-1],nits=nits))
##exclusion with pb
mt09[[3]] <- list(testMods(n09x,nits=nits))
##exclusion w/o pb
mt09[[4]] <- list(testMods(n09x[,-1],nits=nits))
###Output
end <- Sys.time()
time09 <- end-start
print(time09)
dput(mt09,file='../results/mt09.rdata')

###2008
start <- Sys.time()
##control with pb
mt08 <- list(testMods(n08c,nits=nits))
##control w/o pb
mt08[[2]] <- list(testMods(n08c[,-1],nits=nits))
##exclusion with pb
mt08[[3]] <- list(testMods(n08x,nits=nits))
##exclusion w/o pb
mt08[[4]] <- list(testMods(n08x[,-1],nits=nits))
###Output
end <- Sys.time()
time08 <- end-start
print(time08)
dput(mt08,file='../results/mt08.rdata')
