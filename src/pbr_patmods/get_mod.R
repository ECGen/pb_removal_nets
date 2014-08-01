###Running modularity on hoth

library(bipartite)
library(methods)
args <- commandArgs()
x <- read.csv(args[length(args)-1])
y <- computeModules(x,steps=100000)
cat(paste('',slot(y,name='likelihood')),file=args[length(args)],append=TRUE)

