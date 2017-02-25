### PBR one off plots
### MKLau
### 18 Jun 2015

library(gplots)
library(enaR)
source('global.R')
source('pbrDataLoader.R')
library('TeachingDemos')
library('RColorBrewer')
library('ggplot2')
library('bipartite')

coa <- read.csv('../data/pbr_coa.csv')
coa <- coa[c(4,5,8,9),]
x <- unlist(lapply(strsplit(as.character(coa$X),split='\\.'),function(x) x[length(x)]))
x[x == 'avg'] <- 1;x[x == 'npb'] <- 0
trace <- substr(coa$X,1,1)
z <- coa$z
par(cex.lab=1.25,cex.axis=1.25)
interaction.plot(x,trace,z,type='b',pch=c('c','x'),legend=FALSE,xlab=expression(italic('P. betae')),ylab='z-score (Standardized Modularity)',lty=c(1,2))

null <- rnorm(1000,mean=coa[4,5],sd=coa[4,6])
par(cex.lab=1.25,cex.axis=1.25)
hist(null,xlim=c(0,coa[4,4]+0.05),xlab='modularity',ylab='frequency',main='')
abline(v=coa[4,4],lty=2,col='black',lwd=2)

### see the office window
cmSppCavg <- lapply(dir('../data/conModspecies/cmCavg/results/',full.names=TRUE),read.csv)
cmSppXavg <- lapply(dir('../data/conModspecies/cmXavg/results/',full.names=TRUE),read.csv)
names(cmSppCavg) <- dir('../data/conModspecies/cmCavg/results/')
names(cmSppXavg) <- dir('../data/conModspecies/cmXavg/results/')
cmSppCavg <- lapply(cmSppCavg,function(x) x[,2])
cmSppXavg <- lapply(cmSppXavg,function(x) x[,2])
cmSppCavg <- do.call(cbind,cmSppCavg)
cmSppXavg <- do.call(cbind,cmSppXavg)

cmSppCavg <- (cmSppCavg[,sapply(colnames(cmSppCavg),function(x,y) x %in% y,y=colnames(cmSppXavg))])
cmSppXavg <- (cmSppXavg[,sapply(colnames(cmSppXavg),function(x,y) x %in% y,y=colnames(cmSppCavg))])

z.cmSppC <- apply(cmSppCavg,2,function(x,obs) (obs - mean(x)) / sd(x),obs=coa[1,4])
z.cmSppX <- apply(cmSppXavg,2,function(x,obs) (obs - mean(x)) / sd(x),obs=coa[3,4])

p.cmSppC <- apply(cmSppCavg,2,function(x,obs) length(x[x > obs])/length(x),obs=coa[1,4])
p.cmSppX <- apply(cmSppXavg,2,function(x,obs) length(x[x > obs])/length(x),obs=coa[3,4])

par(mfrow=c(1,2),mai=c(1.02*1.5, 0.82, 0.82, 0.42),cex.axis=0.5)
barplot(sort(z.cmSppC),las=2,ylab='Species CM (Control)')
barplot(sort(z.cmSppX),las=2,ylab='Species CM (Exclusion)')

par(mfrow=c(1,1),mai=c(1.02, 0.82, 0.82, 0.42))
axis.min <- floor(min(c(z.cmSppC,z.cmSppX)))
axis.max <- ceiling(max(c(z.cmSppC,z.cmSppX)))
plot(z.cmSppX~z.cmSppC,xlim=c(axis.min,axis.max),ylim=c(axis.min,axis.max),xlab='Species CM (Control)',ylab='Species CM (Exclusion)',pch='')
lines((axis.min:axis.max),(axis.min:axis.max),lwd=0.75)
text(y=z.cmSppX,x=z.cmSppC,labels=names(z.cmSppX),cex=0.75)

t.test(I(z.cmSppX-z.cmSppC))

MU <- c(mean(z.cmSppC),mean(z.cmSppX))
SE <- c(se(z.cmSppC),se(z.cmSppX))

barplot2(MU,plot.ci=TRUE,ci.l= MU-SE,ci.u= MU+SE,names=c('Control','Exclusion'),ylab='Contribution to Modularity (z)')

plot(density(z.cmSppC-z.cmSppX),main='',xlab='Change in Contribution to Modularity (z_C - z_E)')

### Trees
cmCavg <- lapply(dir('../data/conMod/cmCavg/results/',full.names=TRUE),read.csv)
cmXavg <- lapply(dir('../data/conMod/cmXavg/results/',full.names=TRUE),read.csv)
names(cmCavg) <- dir('../data/conMod/cmCavg/results/')
names(cmXavg) <- dir('../data/conMod/cmXavg/results/')
cmCavg <- lapply(cmCavg,function(x) x[,2])
cmXavg <- lapply(cmXavg,function(x) x[,2])
cmCavg <- do.call(cbind,cmCavg)
cmXavg <- do.call(cbind,cmXavg)
cmCavg <- cmCavg[,order(as.numeric(colnames(cmCavg)))]
cmXavg <- cmXavg[,order(as.numeric(colnames(cmXavg)))]
colnames(cmCavg) <- 1:ncol(cmCavg)
colnames(cmXavg) <- 1:ncol(cmXavg)

cmCavg <- (cmCavg[,sapply(colnames(cmCavg),function(x,y) x %in% y,y=colnames(cmXavg))])
cmXavg <- (cmXavg[,sapply(colnames(cmXavg),function(x,y) x %in% y,y=colnames(cmCavg))])

z.cmC <- apply(cmCavg,2,function(x,obs) (obs - mean(x)) / sd(x),obs=coa[1,4])
z.cmX <- apply(cmXavg,2,function(x,obs) (obs - mean(x)) / sd(x),obs=coa[3,4])

par(mfrow=c(1,2),mai=c(1.02*2.5, 0.82, 0.82, 0.42))
barplot(sort(z.cmC),las=2,ylab='Species CM (Control)')
barplot(sort(z.cmX),las=2,ylab='Species CM (Exclusion)')

par(mfrow=c(1,1),mai=c(1.02, 0.82, 0.82, 0.42))
axis.min <- floor(min(c(z.cmC,z.cmX)))
axis.max <- ceiling(max(c(z.cmC,z.cmX)))
plot(z.cmX~z.cmC,xlim=c(axis.min,axis.max),ylim=c(axis.min,axis.max),xlab='Species CM (Control)',ylab='Species CM (Exclusion)',pch=19)
lines((axis.min:axis.max),(axis.min:axis.max),lwd=0.75)

### genotypes?
cmgeno <- factor(as.character(geno.09$x[as.numeric(colnames(cmCavg))]))
z.cm <- c(z.cmX,z.cmC)
trt <- c(rep('X',length(z.cmX)),rep('C',length(z.cmC)))
geno <- factor(c(as.character(cmgeno),as.character(cmgeno)))
pb <- unlist(pb.A)

cm.data <- list(c=data.frame(geno=cmgeno,pb=pb.A$c,cm=z.cmC),x=data.frame(geno=cmgeno,pb=pb.A$x,cm=z.cmX))
anova(glm(I(cm^2) ~ geno * pb,data=cm.data$c),test='Chi')
anova(glm(I(cm^2) ~ geno * pb,data=cm.data$x),test='Chi')


par(cex.lab=1.25,cex.axis=1.25)
interaction.plot(trt,geno,z.cm,type='b')

### cross hair plots
ch.col <- rep('black',length(z.cmC))
ch.col <- rainbow(max(as.numeric(cmgeno)))[as.numeric(cmgeno)]
ch.pch <- rep(19,length(z.cmC))
chPlot(cbind(z.cmC,z.cmX),f=cmgeno,col=ch.col,pch=ch.pch,xlim=c(0,1.5),ylim=c(0,1.5))
lines(seq(0,1.5,by=0.5),seq(0,1.5,by=0.5),lty=2)
legend('bottomright',legend=levels(cmgeno),col=unique(ch.col),pch=19)

### Does the contribution to modularity of pb vary by genotype?
### Correlate PB abundance, Removal effect, richness and Z 

pb.A <- list((pbr.08$c[,1]+pbr.09$c[,1])/2,(pbr.08$x[,1]+pbr.09$x[,1])/2)
names(pb.A) <- c('c','x')
richness <- list(apply(pbr.09$c,1,function(x) sum(sign(x))),apply(pbr.09$x,1,function(x) sum(sign(x))))
names(richness) <- c('c','x')

### ANCOVA
abundance <- unlist(pb.A)
Z <- c(z.cmC,z.cmX)
genotype <- unlist(geno.09)

X <- cbind(Z,A=abundance)
X <- apply(X,2,function(x) x/max(x))
tree <- rep(1:length(z.cmC),2)

summary(lm(I((Z^2)) ~ trt*tree/abundance))

z.data <- data.frame(tree=1:length(z.cmC),geno=geno.09$c,pb.a=I(pb.A$c-pb.A$x),cm.c=z.cmC,cm.x=z.cmX)
z.data <- z.data[!(z.data$geno %in% c(1008,1020)),]
z.data <- make.rm(constant=c("tree","pb.a","geno"),repeated=c("cm.c","cm.x"),data=z.data)

summary(aov(repdat~contrasts*geno*pb.a+Error(tree),z.data))

z.data <- data.frame(tree=1:length(z.cmC),geno=geno.09$c,pb.a=I(pb.A$c-pb.A$x),cm.c=z.cmC,cm.x=z.cmX,z.d=I(abs(z.cmC-z.cmX)))

lmer(cm ~ pb * (1 | geno))

### Plots
ch.dat <- list(c=data.frame(cmgeno,pb.A$c,z.cmC),x=data.frame(cmgeno,pb.A$x,z.cmX),s=data.frame(cmgeno,pb.A$c,z.cmC)[cmgeno %in% c('996','1008','1020') == FALSE,],d=data.frame(cmgeno,(pb.A$c - pb.A$x),(z.cmX - z.cmC)))

g.pch <- toupper(as.character(levels(cmgeno))[c(2,4,5,1,8,9,7,10,6,3)])
g.pch[g.pch == 'COAL 3'] <- 'Coal-3'
g.pch <- cbind(g.pch,pch = c(), col = c())
g.pch 
g.pch 

### ch.dat <- lapply(ch.dat,function(x) x[x[,1] %in% c(1008,1020,996) == FALSE,])
ch.col <- brewer.pal(n=max(as.numeric(ch.dat$c[,1]))+2, name='Set3')[as.numeric(ch.dat$c[,1])+2]
ch.key <- data.frame(pch=c(19,19,25,24,22,22,23,23,24,25),geno=c(996,1000,1008,1017,1020,'coal 3','HE-10','Rm-2','T-15','WC-5'),col=c('black','grey','darkgrey','grey','black','darkgrey','black','grey','darkgrey','grey'))
ch.col <- as.character(cmgeno)
ch.pch <- as.character(cmgeno)
for (i in 1:nrow(ch.key)){
    ch.col[cmgeno == ch.key$geno[i]] <- as.character(ch.key$col[i])
    ch.pch[cmgeno == ch.key$geno[i]] <- ch.key$pch[i]
}
ch.pch <- as.numeric(ch.pch)


par(mfrow=c(1,1))
chPlot(ch.dat$x[,2:3],f=ch.dat$x[,1],col=rep('darkgrey',length(ch.col)),pch=rep(19,length(ch.pch)),xlim=c(0,75),ylim=c(0,1.5),se=TRUE,line.lm=TRUE,line.col='darkgrey',line.lty=1,xlab=expression(italic('Pemphigus betae')~' abundance'),ylab='Tree genotype contribution to modularity (Z)',cex=1.5)
abline(v=max(tapply(pb.A$x,cmgeno,mean)+tapply(pb.A$x,cmgeno,se)),col='lightgrey',lty=1)
chPlot(ch.dat$s[,2:3],f=ch.dat$s[,1],col=rep(1,nrow(ch.dat$s)),pch=rep(19,nrow(ch.dat$s)),xlim=c(0,75),ylim=c(0,1.5),se=TRUE,line.lm=TRUE,line.col='darkgrey',line.lty=1,xlab='Aphid Abundance',ylab='z (modularity)',cex=1.5,add=TRUE)
chPlot(ch.dat$c[,2:3],f=ch.dat$c[,1],col=rep(1,length(ch.col)),pch=rep(19,length(ch.pch)),xlim=c(0,75),ylim=c(0,1.5),se=TRUE,line.lm=TRUE,line.col='black',xlab='Aphid Abundance',ylab='z (modularity)',cex=1.5,add=TRUE,line.lty=2)
legend('topright',legend=c(expression(italic('P. betae')),'Present','Excluded'),col=c(1,1,'darkgrey'),pch=c(30,19,19),bg='white',box.col='black')

### network plots
net.c <- floor(meanMat(pbr.08$c,pbr.09$c))
net.x <- floor(meanMat(pbr.08$x,pbr.09$x))

### unimodal representation of the bipartite networks
## uni <- lapply(list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x),unipart)
uni <- lapply(list(c=net.c,x=net.x),unipart,rm.zero=TRUE,std=FALSE,thresh=0.0001)
## uni <- lapply(list(c08=pbr.08$c,x08=pbr.08$x,c09=pbr.09$c,x09=pbr.09$x),cdNet,alpha=0.001)
uni <- lapply(uni,function(x) x[order(apply(x,1,sum),decreasing=TRUE),order(apply(x,2,sum),decreasing=TRUE)])

uni.sub <- lapply(uni,function(x,n) x[1:n,1:n],n=35)
uni.col <- lapply(uni,colnames)
for (i in 1:length(uni.col)){
    uni.col[[i]][tolower(uni.col[[i]]) == 'pb'] <- 'black';uni.col[[i]][tolower(uni.col[[i]]) != 'black'] <- 'darkgrey'
}

cen <- lapply(uni,evcent,rescale=TRUE)
cen.sub <- lapply(uni,function(x,n) x[1:n],n=35)

par(mfrow=c(1,2),mai=c(0,0,0.5,0))
pc <- gplot.target(uni.sub$c,cen.sub$c,gmode='graph',circ.col='darkgrey',circ.lab=FALSE,vertex.col='white',displaylabels=FALSE,edge.col='lightgrey',edge.lwd=0.01,vertex.border='white',main='Aphid Present',vertex.cex=1)
points(pc,col=uni.col$c,pch=19,cex=1.2)
points(pc,col=uni.col$c,pch=19,cex=as.numeric(uni.col$c == 'black'))

px <- gplot.target(uni.sub$x,cen.sub$x,gmode='graph',circ.col='darkgrey',circ.lab=FALSE,vertex.col='white',displaylabels=FALSE,edge.col='lightgrey',edge.lwd=0.01,vertex.border='white',main='Aphid Present',vertex.cex=1)
points(px,col=uni.col$x,pch=19,cex=1.2)
points(px,col=uni.col$x,pch=19,cex=as.numeric(uni.col$x == 'black'))

deg <- lapply(uni,degree,rescale=TRUE)
d.deg <- deg[[2]] - deg[[1]]
d.deg <- (d.deg - mean(d.deg))/sd(d.deg)
plot(d.deg,pch='')
text(d.deg,labels=colnames(uni[[1]]))

### bipartite representation
net.thresh <- 2

net <- floor(meanMat(pbr.08$c,pbr.09$c))
net <- as.matrix(net)
net <- net[,na.omit(match(names(z.cmSppC),colnames(net)))] ###match species

net[net < net.thresh] <- 0

rownames(net) <- paste(as.character(geno.09$c),1:nrow(net),sep='_')

tree.net <- net %*% t(net)
tree.net <- tree.net/max(tree.net)
spp.net <- t(net) %*% net
spp.net <- spp.net/max(spp.net)

tree.net[tree.net < 0.01] <- 0

tree.col <- ch.col 
spp.col <- 'red' ###

cm.top3 <- names(z.cmSppC)[order(I(z.cmSppC-z.cmSppX))][1:3]
spp.lab <- colnames(net)
spp.lab[spp.lab%in%cm.top3 == FALSE] <- ''

gplot(tree.net,gmode='graph',edge.lwd=tree.net,vertex.col=tree.col)
gplot(spp.net,gmode='graph',edge.lwd=spp.net,vertex.col=spp.col,displaylabels=TRUE,label=spp.lab,edge.col='lightgrey')

plotweb(sortMat(net),text.rot=90,
        col.low=tree.col[order(apply(net,1,sum),decreasing=TRUE)],
        col.high=spp.col[order(apply(net,2,sum),decreasing=TRUE)],
        method='normal',
        )


### Final stats 8Feb2016

## 1. Network modularity was 
coa[c(2,4),]

## 2. PB increased modularity at the scale of all trees
## randomizing PB leads to community networks that are
## less modular
sort(z.cmSppC[abs(z.cmSppC) > 1.5])
sort(z.cmSppX[abs(z.cmSppX) > 1.5])

sort(p.cmSppC[abs(z.cmSppC) > 1.5])
sort(p.cmSppX[abs(z.cmSppX) > 1.5])


## 3. PB impacts tree genotype modularity

anova(glm(I(cm^2) ~ geno * pb,data=cm.data$c),test='F')
anova(glm(I(cm^2) ~ geno * pb,data=cm.data$x),test='F')

