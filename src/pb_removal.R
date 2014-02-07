###Analysis of Art's Pemphigus removal data
###Current file started 18 Oct 2013

###Import data
add.pb <- TRUE
rm.resistant <- TRUE
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
                                        #remove species occurring less than 5 times
pbr.08 <- pbr.08[,apply(pbr.08,2,sum)>=5]
pbr.09 <- pbr.09[,apply(pbr.09,2,sum)>=5]
                                        #separate by treatment
pbr.08 <- split(pbr.08,trt.08)
pbr.09 <- split(pbr.09,trt.09)
geno.08 <- split(geno.08,trt.08)
geno.09 <- split(geno.09,trt.09)
                                        #add pb?
if (add.pb){
  library(xlsx)
  pb08 <- as.matrix(read.xlsx('../data/P betae excl treatment effect 2008.xlsx',sheetIndex=1))
  colnames(pb08) <- as.character(pb08[1,])
  pb08 <- pb08[-1,]
  pb08 <- data.frame(genotype=pb08[,1],tree=pb08[,2],x=as.numeric(pb08[,3]),c=as.numeric(pb08[,4]))
  pb08.tree <- c(as.character(pb08[,2]),as.character(pb08[,2]))
  pb08.trt <- c(rep('x',nrow(pb08)),rep('c',nrow(pb08)))
  pb <- c(pb08[,3],pb08[,4])
  trt.tree.pbr <- paste(tree.08,trt.08)
  trt.tree.pb <- paste(pb08.tree,pb08.trt)
  pb <- pb[match(trt.tree.pbr,trt.tree.pb)]
  check <- trt.tree.pb[match(trt.tree.pbr,trt.tree.pb)]
  if (all(check==trt.tree.pbr)){print('You are good to go!')}else{warning('Oh no!!!')}
  pb <- split(pb,trt.08)
  pbr.08[[1]] <- cbind(pb=pb[[1]],pbr.08[[1]])
  pbr.08[[2]] <- cbind(pb=pb[[2]],pbr.08[[2]])
                                        #
  pb09 <- as.matrix(read.xlsx('../data/2009Pbetaetreatmenteffect.xls',sheetIndex=1))
  pb09 <- data.frame(genotype=pb09[,1],tree=pb09[,2],x=as.numeric(pb09[,3]),c=as.numeric(pb09[,4]))
  pb09[pb09[,3]==1,3] <- 'x'
  pb09[pb09[,3]==2,3] <- 'c'
  pb09 <- pb09[match(paste(tree.09,trt.09),paste(pb09[,2],pb09[,3])),]
  pb <- split(pb09[,4],pb09[,3])
  if (all(split(pb09[,2],pb09[,3])[[1]]==split(tree.09,trt.09)[[1]])){print('Good to go!')}else{print('Holy crap!')}
  pbr.09[[1]] <- cbind(pb=pb[[1]],pbr.09[[1]])
  pbr.09[[2]] <- cbind(pb=pb[[2]],pbr.09[[2]])
}else{}
                                        #make binary
as.binary <- function(x){x[x!=0] <- 1;return(x)}
co.08 <- lapply(pbr.08,as.binary)
co.09 <- lapply(pbr.09,as.binary)
                                        #remove resistants
if (rm.resistant){
  co.08[[1]] <- co.08[[1]][geno.08[[1]]!='1008'&geno.08[[1]]!='1020',]
  co.08[[2]] <- co.08[[2]][geno.08[[2]]!='1008'&geno.08[[2]]!='1020',]
  co.09[[1]] <- co.09[[1]][geno.09[[1]]!='1008'&geno.09[[1]]!='1020',]
  co.09[[2]] <- co.09[[2]][geno.09[[2]]!='1008'&geno.09[[2]]!='1020',]
}else{}

###Networks
library(sna)
library(ComGenR)
source('~/projects/pb_removal_nets/src/helper_funcs.R')
                                        #
print('Modeling networks')
net.08 <- lapply(co.08,co.net)
net.09 <- lapply(co.09,co.net)
                                        #reduced nets
                                        #reduce speceis names to numbers
rownames(net.08[[1]]) <- colnames(net.08[[1]]) <- as.character(1:ncol(net.08[[1]]))
rownames(net.08[[2]]) <- colnames(net.08[[2]]) <- as.character(1:ncol(net.08[[2]]))
rownames(net.09[[1]]) <- colnames(net.09[[1]]) <- as.character(1:ncol(net.09[[1]]))
rownames(net.09[[2]]) <- colnames(net.09[[2]]) <- as.character(1:ncol(net.09[[2]]))
                                        #correlations
edge.08 <- list(net.08[[1]][net.08[[1]]!=0|net.08[[2]]!=0],net.08[[2]][net.08[[1]]!=0|net.08[[2]]!=0])
names(edge.08) <- names(net.08)
plot(edge.08[[2]]~edge.08[[1]])
edge.09 <- list(net.09[[1]][net.09[[1]]!=0|net.09[[2]]!=0],net.09[[2]][net.09[[1]]!=0|net.09[[2]]!=0])
names(edge.09) <- names(net.09)
plot(edge.09[[2]]~edge.09[[1]])
                                        #gplots
par(mfcol=c(2,2))
coord.08=mgp2(net.08[[1]],v.cex=(apply(pbr.08[[1]],2,sum)/max(apply(pbr.08[[1]],2,sum))),scalar=1,displaylabels=TRUE)
title(main='2008')
mtext(side=2,text='Control',font=2)
mgp2(net.08[[2]],v.cex=(apply(pbr.08[[2]],2,sum)/max(apply(pbr.08[[2]],2,sum))),scalar=1,my.coord=coord.08,displaylabels=TRUE)
mtext(side=2,text='Removed',font=2)
coord.09=mgp2(net.09[[1]],v.cex=(apply(pbr.09[[1]],2,sum)/max(apply(pbr.09[[1]],2,sum))),scalar=1,displaylabels=TRUE)
title(main='2009')
mgp2(net.09[[2]],v.cex=(apply(pbr.09[[2]],2,sum)/max(apply(pbr.09[[2]],2,sum))),scalar=1,my.coord=coord.09,displaylabels=TRUE)
                                        #
par(mfrow=c(1,2))
hist(nodeDist(net.08),xlim=c(0,max(nodeDist(net.08))),main='2008')
abline(v=mean(nodeDist(net.08)),lty=2)
hist(nodeDist(net.09),xlim=c(0,max(nodeDist(net.09))),main='2009')
abline(v=mean(nodeDist(net.09)),lty=2)
t.test(nodeDist(net.08))
t.test(nodeDist(net.09))

##                                         #qap test

## g.08 <- array(0,dim=c(nrow(net.08[[1]]),ncol(net.08[[1]]),2))
## g.08[,,1] <- net.08[[1]]
## g.08[,,2] <- net.08[[2]]
## g.09 <- array(0,dim=c(nrow(net.09[[1]]),ncol(net.09[[1]]),2))
## g.09[,,1] <- net.09[[1]]
## g.09[,,2] <- net.09[[2]]
##                                         #test
## print('running qap 2008')
## qap.08 <- qaptest(g.08,gcor,g1=1,g2=2)
## dput(qap.08,file='../results/qap08.Rdata')
## print('running qap 2009')
## qap.09 <- qaptest(g.09,gcor,g1=1,g2=2)
## dput(qap.09,file='../results/qap09.Rdata')
## print('gcor matrices')
## dput(gcor(g.08),file='../data/gcor08.Rdata')
## dput(gcor(g.09),file='../data/gcor09.Rdata')
## print('Done!')
                                        #qap
library(sna)
                                        #see qap_test.R for how analysis was run
readLines('../results/qap_results_08.txt')
readLines('../results/qap_results_09.txt')
qap.d.08 <- read.csv('../results/qap_dist_08.csv')
qap.d.09 <- read.csv('../results/qap_dist_09.csv')
hist(qap.d.08[,1])
hist(qap.d.09[,1])
                                        #nestedness plots and looking at node weights
                                        # percent max convert
library(gplots)
smc.09c <- apply(pbr.09$c,2,function(x) if (any(x!=0)){x/max(x)}else{x})
tot.09c <- apply(smc.09c,1,sum)
mu <- tapply(tot.09c,geno.09$c,mean)
se <- tapply(tot.09c,geno.09$c,function(x) sd(x)/sqrt(length(x)))
1000      1008      1017      1020       996    coal 3     HE-10      Rm-2      T-15      WC-5 
1           2       3          4          5       6          7         8          9        10
deg.09c <- c(3,6,7,1,10,8,5,2,4,9)
se <- se[deg.09c]
mu <- mu[deg.09c]
barplot2(mu,plot.ci=TRUE,ci.u=mu+se,ci.l=mu-se,ylab='Toal Percent Species Maximums')
                                        #
source('~/projects/packages/ComGenR_development/src/cgREML.R')
round(cgREML(tot.09c,geno.09$c),5)
