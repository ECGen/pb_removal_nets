###Testing for the modularity of the network
##18 Feb 2014

library(bipartite)
library(methods)
                                        #
pbr.08 <- read.csv('../data/keith_pb_removal_2008.csv')
pbr.09 <- read.csv('../data/keith_pb_removal_2009.csv')
pbr.08[,2] <- as.character(pbr.08[,2])
pbr.09[,2] <- as.character(pbr.09[,2])
                                        #fix typos
pbr.08[pbr.08[,2]=="SSG1-1c",2] <- "SSG1-1 c" 
pbr.08[pbr.08[,2]=="S7-3x",2] <- "S7-3 x"
pbr.08[pbr.08[,2]=="S7-3c",2] <- "S7-3 c"
pbr.08[pbr.08[,2]=="S4-9c",2] <- "S4-9 c"
pbr.08[pbr.08[,2]=="S3-21 x ",2] <- "S3-21 x"
pbr.08[pbr.08[,2]=="S3-21 c ",2] <- "S3-21 c"
pbr.09[pbr.09[,2]=="SSG1-1c",2] <- "SSG1-1 c" 
pbr.09[pbr.09[,2]=="S7-3x",2] <- "S7-3 x"
pbr.09[pbr.09[,2]=="S7-3c",2] <- "S7-3 c"
pbr.09[pbr.09[,2]=="S4-9c",2] <- "S4-9 c"
pbr.09[pbr.09[,2]=="S3-21 x ",2] <- "S3-21 x"
pbr.09[pbr.09[,2]=="S3-21 c ",2] <- "S3-21 c"
                                        #
pbr.08[is.na(pbr.08)] <- 0
pbr.09[is.na(pbr.09)] <- 0
geno.08 <- pbr.08[,1]
geno.09 <- pbr.09[,1]
tree.08 <- as.character(sapply(as.character(pbr.08[,2]),function(x) strsplit(x,split=' ')[[1]][1]))
tree.09 <- as.character(sapply(as.character(pbr.09[,2]),function(x) strsplit(x,split=' ')[[1]][1]))
trt.08 <- as.character(sapply(as.character(pbr.08[,2]),function(x) strsplit(x,split=' ')[[1]][2]))
trt.09 <- as.character(sapply(as.character(pbr.09[,2]),function(x) strsplit(x,split=' ')[[1]][2]))
                                        #remove env data
pbr.08 <- pbr.08[,-1:-2]
pbr.09 <- pbr.09[,-1:-2]
                                        #separate by treatment
pbr.08 <- split(pbr.08,trt.08)
pbr.09 <- split(pbr.09,trt.09)
geno.08 <- split(geno.08,trt.08)
geno.09 <- split(geno.09,trt.09)
                                        #add pb?
pb08 <- read.csv('../data/P betae excl treatment effect 2008.csv')
colnames(pb08) <- as.character(pb08[1,])
pb08 <- data.frame(genotype=pb08[,1],tree=pb08[,2],x=as.numeric(pb08[,3]),c=as.numeric(pb08[,4]))
pb08.tree <- c(as.character(pb08[,2]),as.character(pb08[,2]))
pb08.trt <- c(rep('x',nrow(pb08)),rep('c',nrow(pb08)))
pb <- c(pb08[,3],pb08[,4])
trt.tree.pbr <- paste(tree.08,trt.08)
trt.tree.pb <- paste(pb08.tree,pb08.trt)
pb <- pb[match(trt.tree.pbr,trt.tree.pb)]
check <- trt.tree.pb[match(trt.tree.pbr,trt.tree.pb)]
if (all(check==trt.tree.pbr)){print('Good to go!')}else{warning('Oh no!!!')}
pb <- split(pb,trt.08)
pbr.08[[1]] <- cbind(pb=pb[[1]],pbr.08[[1]])
pbr.08[[2]] <- cbind(pb=pb[[2]],pbr.08[[2]])
                                        #
pb09 <- read.csv('../data/2009Pbetaetreatmenteffect.csv')
pb09 <- data.frame(genotype=pb09[,1],tree=pb09[,2],x=as.numeric(pb09[,3]),c=as.numeric(pb09[,4]))
pb09[pb09[,3]==1,3] <- 'x'
pb09[pb09[,3]==2,3] <- 'c'
pb09 <- pb09[match(paste(tree.09,trt.09),paste(pb09[,2],pb09[,3])),]
pb <- split(pb09[,4],pb09[,3])
if (all(split(pb09[,2],pb09[,3])[[1]]==split(tree.09,trt.09)[[1]])){print('Good to go!')}else{print('Holy crap!')}
pbr.09[[1]] <- cbind(pb=pb[[1]],pbr.09[[1]])
pbr.09[[2]] <- cbind(pb=pb[[2]],pbr.09[[2]])
                                        #
m <- sign(pbr.09$c)[,apply(pbr.09$c,2,sum)!=0]
n <- 1000
print('Null Modules')
null <- list()
mod.null <- list()
for (i in 1:n){null[[i]] <- commsimulator(m,method='r1',thin=100);print(i)}
for (i in 1:length(null)){mod.null[[i]] <- slot(computeModules(null[[i]]),name='likelihood');print(i)}
print('Observed Modules')
mod.null <- unlist(mod.null)
mod.obs <- computeModules(m)
z <- (mod.obs-mean(mod.null))/sd(mod.null)
p <- c(lower.p=length(mod.null[mod.null<=mod.obs]),upper.p=length(mod.null[mod.null>=mod.obs]))/n
out <- list(null=mod.null,results=c(obs=mod.obs,z=z,p))
print('saving...')
dput(out,file='../results/pbr_mod_out.R')

