###Null community generator for pbr mods
nits <- 100
n08c <- dget('./data/pbr08c.rdata')
n08x <- dget('./data/pbr08x.rdata')
n09c <- dget('./data/pbr09c.rdata')
n09x <- dget('./data/pbr09x.rdata')

###Create directories
dir.create('./data/null08c') 
dir.create('./data/null08cnpb') 
dir.create('./data/null08x') 
dir.create('./data/null08xnpb') 
dir.create('./data/null09c') 
dir.create('./data/null09cnpb') 
dir.create('./data/null09x') 
dir.create('./data/null09xnpb')

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
out <- r2dtable(nits,apply(n08x,1,sum),apply(n08x,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08x/',i,sep=''),row.names=FALSE)
}

##exclusion w/o pb
out <- r2dtable(nits,apply(n08x[,-1],1,sum),apply(n08x[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null08xnpb/',i,sep=''),row.names=FALSE)
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
out <- r2dtable(nits,apply(n09x,1,sum),apply(n09x,2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09x/',i,sep=''),row.names=FALSE)
}

##exclusion w/o pb
out <- r2dtable(nits,apply(n09x[,-1],1,sum),apply(n09x[,-1],2,sum))
for (i in 1:length(out)){
  write.csv(out[[i]],file=paste('./data/null09xnpb/',i,sep=''),row.names=FALSE)
}
