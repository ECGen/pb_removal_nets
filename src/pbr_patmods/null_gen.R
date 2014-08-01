###Null community generator for pbr mods
nits <- 100
n08c <- dget('./data/pbr08c')
n08e <- dget('./data/pbr08e')
n09c <- dget('./data/pbr09c')
n09e <- dget('./data/pbr09e')

###2008
##control with pb
out <- r2dtable(nits,apply(n08c,1,sum),apply(n08c,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08c/',i,sep=''),row.names=FALSE)
}

##control w/o pb
out <- r2dtable(nits,apply(n08c[,-1],1,sum),apply(n08c[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08cnpb/',i,sep=''),row.names=FALSE)
}

##exclusion with pb
out <- r2dtable(nits,apply(n08e,1,sum),apply(n08e,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08e/',i,sep=''),row.names=FALSE)
}

##exclusion w/o pb
out <- r2dtable(nits,apply(n08e[,-1],1,sum),apply(n08e[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08enpb/',i,sep=''),row.names=FALSE)
}

###2009
##control with pb
out <- r2dtable(nits,apply(n09c,1,sum),apply(n09c,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09c/',i,sep=''),row.names=FALSE)
}

##control w/o pb
out <- r2dtable(nits,apply(n09c[,-1],1,sum),apply(n09c[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09cnpb/',i,sep=''),row.names=FALSE)
}

##exclusion with pb
out <- r2dtable(nits,apply(n09e,1,sum),apply(n09e,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09e/',i,sep=''),row.names=FALSE)
}

##exclusion w/o pb
out <- r2dtable(nits,apply(n09e[,-1],1,sum),apply(n09e[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09enpb/',i,sep=''),row.names=FALSE)
}
