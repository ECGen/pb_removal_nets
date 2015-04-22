### text mining with Aaron's papers
### 2Apr2015
library(wordcloud)
source('./tmFuncs.R')

add.pb <- TRUE
rm.resistant <- FALSE
binary <- FALSE
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
pbr.08 <- pbr.08[,apply(pbr.08,2,sum)>=10]
pbr.09 <- pbr.09[,apply(pbr.09,2,sum)>=10]
                                        #separate by treatment
pbr.08 <- split(pbr.08,trt.08)
pbr.09 <- split(pbr.09,trt.09)
geno.08 <- split(geno.08,trt.08)
geno.09 <- split(geno.09,trt.09)
                                        #add pb?
if (add.pb){
  library(gdata)
  pb08 <- as.matrix(read.xls('../data/P betae excl treatment effect 2008.xlsx',sheet=1))
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
  pb09 <- as.matrix(read.xls('../data/2009Pbetaetreatmenteffect.xls',sheet=1))
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
if (binary){
    co.08 <- lapply(pbr.08,as.binary)
    co.09 <- lapply(pbr.09,as.binary)
}else{
    co.08 <- pbr.08;co.09 <- pbr.09
}
                                        #remove resistants
if (rm.resistant){
  co.08[[1]] <- co.08[[1]][geno.08[[1]]!='1008'&geno.08[[1]]!='1020',]
  co.08[[2]] <- co.08[[2]][geno.08[[2]]!='1008'&geno.08[[2]]!='1020',]
  co.09[[1]] <- co.09[[1]][geno.09[[1]]!='1008'&geno.09[[1]]!='1020',]
  co.09[[2]] <- co.09[[2]][geno.09[[2]]!='1008'&geno.09[[2]]!='1020',]
}else{}

tmt.08 <- lapply(co.08,function(x) apply(x,2,sum))
tmt.09 <- lapply(co.09,function(x) apply(x,2,sum))

pdf('../results/pbrWC08c.pdf')
wordcloud(names(tmt.08$c),as.numeric(tmt.08$c), 
          scale=c(7,0.15), random.order=FALSE, 
          rot.per=0.35, use.r.layout=FALSE, 
          colors=brewer.pal(8,"Dark2"),min.freq=1)
dev.off()
system('open ../results/pbrWC08c.pdf')

pdf('../results/pbrWC08x.pdf')
wordcloud(names(tmt.08$x),as.numeric(tmt.08$x), 
          scale=c(7,0.15), random.order=FALSE, 
          rot.per=0.35, use.r.layout=FALSE, 
          colors=brewer.pal(8,"Dark2"),min.freq=1)
dev.off()
system('open ../results/pbrWC08x.pdf')

pdf('../results/pbrWC09c.pdf')
wordcloud(names(tmt.09$c),as.numeric(tmt.09$c), 
          scale=c(7,0.15), random.order=FALSE, 
          rot.per=0.35, use.r.layout=FALSE, 
          colors=brewer.pal(8,"Dark2"),min.freq=1)
dev.off()
system('open ../results/pbrWC09c.pdf')

pdf('../results/pbrWC09x.pdf')
wordcloud(names(tmt.09$x),as.numeric(tmt.09$x), 
          scale=c(7,0.15), random.order=FALSE, 
          rot.per=0.35, use.r.layout=FALSE, 
          colors=brewer.pal(8,"Dark2"),min.freq=1)
dev.off()
system('open ../results/pbrWC09x.pdf')

