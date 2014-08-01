###PBR: modularity analysis on null communities
library(bipartite)
library(methods)

##2008 c p
x <- list()
for (i in 1:length(dir(paste('./data/null08c/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null08c/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods08cp.rdata',append=TRUE)
}

##2008 c n
x <- list()
for (i in 1:length(dir(paste('./data/null08cnpb/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null08cnpb/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods08cn.rdata',append=TRUE)
}

##2008 x p
x <- list()
for (i in 1:length(dir(paste('./data/null08x/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null08x/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods08xp.rdata',append=TRUE)
}

##2008 x n
x <- list()
for (i in 1:length(dir(paste('./data/null08xnpb/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null08xnpb/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods08xn.rdata',append=TRUE)
}

##2009 c p
x <- list()
for (i in 1:length(dir(paste('./data/null09c/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null09c/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods09cp.rdata',append=TRUE)
}

##2009 c n
x <- list()
for (i in 1:length(dir(paste('./data/null09cnpb/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null09cnpb/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods09cn.rdata',append=TRUE)
}

##2009 x p
x <- list()
for (i in 1:length(dir(paste('./data/null09x/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null09x/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods09xp.rdata',append=TRUE)
}

##2009 x n
x <- list()
for (i in 1:length(dir(paste('./data/null09xnpb/',sep='')))){
  x[[i]] <- read.csv(paste('./data/null09xnpb/',i,sep=''))
}

for (i in 1:length(x)){
  y <- computeModules(x[[i]])
  cat(paste('',slot(y,name='likelihood')),file='./data/mods09xn.rdata',append=TRUE)
}
