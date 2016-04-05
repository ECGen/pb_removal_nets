library(grid)
library(bipartite)
library(RCurl)
library(XML)
library(vegan)
library(ggnet)

unipart <- function(x,rm.zero=TRUE,std=TRUE,thresh=0.05){
    out <- t(as.matrix(sign(x))) %*% 
        as.matrix(sign(x))
    if (std){out <- out / max(out)}
    if (thresh == FALSE){thresh <- 0}else{
        if (thresh > 1 | thresh < 0){
            warning('Threshold outside bounds (0-1).')
            thresh[thresh > 1] <- 1
            thresh[thresh < 0] <- 0
        }
        h <- hist(out,plot=FALSE)
        h <- head(h$breaks[h$density > thresh],2)
        if (length(h) < 2 | any(is.na(h))){h <- 0}else{
            h <- tail(h,1)
        }
    }
    out[out < h] <- 0
    if (rm.zero){out <- out[apply(out,1,sum) != 0, apply(out,2,sum) != 0]}
    return(out)
}

se <- function(x){sd(x)/sqrt(length(x))}

gitPush <- function(x='.',message='update'){
    system(paste('git add -A ',x,sep=''))
    system(paste('git commit -am ',message,sep=''))
    system('git push')
}

alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3], alpha=alpha))
}

getElev <- function(latitude=52.4822,longitude=-1.8946){
    url <- paste(
        "http://www.earthtools.org/height",
        latitude,
        longitude,
        sep = "/"
        )
    page <- getURL(url)
    ans <- xmlTreeParse(page, useInternalNodes = TRUE)
    heightNode <- xpathApply(ans, "//meters")[[1]]
    as.numeric(xmlValue(heightNode))
}

chPlot <- function(x,f,col,pch,se=FALSE,xlim=c(-1,1),ylim=c(-1,1),xlab='X',ylab='Y',add.plot=FALSE,line.lm=FALSE,line.lty=1,line.col=1,cex=1){
    col <- tapply(col,f,function(x) x[1])
    pch <- tapply(pch,f,function(x) x[1])
    mu <- apply(x,2,function(x,f) tapply(x,f,mean),f=f)
    if (se){
        bars <- apply(x,2,function(x,f) 
            tapply(x,f,function(x) sd(x)/sqrt(length(x))),f=f)
    }else{
        bars <- apply(x,2,function(x,f) 
            tapply(x,f,sd),f=f)
    }
    bar.up <- mu + bars
    bar.lo <- mu - bars
### make the plot
    if (add.plot){
        points(mu,col=col,pch=pch,cex=cex)
    }else{
        plot(mu,col=col,pch=pch,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,cex=cex)
    }
    for (i in 1:nrow(mu)){
        lines(c(bar.up[i,1],bar.lo[i,1]),rep(mu[i,2],2),col=col[i])
        lines(rep(mu[i,1],2),c(bar.up[i,2],bar.lo[i,2]),col=col[i])
    }
    if (line.lm){
        abline(lm(mu[,2]~mu[,1]),lty=line.lty,col=line.col)
    }
    return(list(col=col,pch=pch))
}

ch.point <- function(x,f){
    list(mu=tapply(x,f,mean),se=tapply(x,f,function(x) sd(x)/sqrt(length(x))))
}

ch.vector <- function(x,t,pch,col){
    mu <- apply(x,2,function(x,t) tapply(x,t,mean),t=t)
    se <- apply(x,2,function(x,t) tapply(x,t,function(x) sd(x)/sqrt(length(x))),t=t)
    cbind(mu,se)
}


coStats <- function(obs,null){
    z <- (obs - mean(null)) / sd(null)
    if (z <= 0){
        pvalue <- length(null[null <= obs])/length(null)
    }else if (z > 0){
        pvalue <- length(null[null >= obs])/length(null)
    }
    out <- c(z=z,pvalue=pvalue,obs=obs,mean=mean(null),sd=sd(null))
    return(out)
}

meanMat <- function(x,y){
    (equalMat(x,y)[[1]] + equalMat(x,y)[[2]]) / 2
}

equalMat <- function(x,y){
    x <- list(x,x)
    y <- list(y,y)
    x[[2]] <- matrix(0,nrow=nrow(x[[1]]),
                           ncol=length(
                               colnames(y[[1]])[colnames(y[[1]]) 
                                                %in% 
                                                colnames(x[[1]]) == FALSE])
                     )
    colnames(x[[2]]) <- colnames(y[[1]])[colnames(y[[1]]) 
                                         %in% 
                                         colnames(x[[1]]) == FALSE]
    x[[2]] <- cbind(x[[1]],x[[2]])

    y[[2]] <- matrix(0,nrow=nrow(y[[1]]),
                           ncol=length(
                               colnames(x[[1]])[colnames(x[[1]]) 
                                                %in% 
                                                colnames(y[[1]]) == FALSE])
                     )
    colnames(y[[2]]) <- colnames(x[[1]])[colnames(x[[1]]) 
                                         %in% 
                                         colnames(y[[1]]) == FALSE]
    y[[2]] <- cbind(y[[1]],y[[2]])
    return(list(x=x[[2]],y=y[[2]]))
}

mkNull <- function(x,n,null.dir){
    for (i in 1:n){
        write.table(nullmodel(x,method='r2d',N=1),
                    file=paste(null.dir,i,sep='/'),
                    row.names=FALSE,col.names=FALSE,
                    sep=',')
    }
    print("Done!")
}

sortMat <- function(x){
    x[order(apply(x,1,sum),decreasing=TRUE),order(apply(x,2,sum),decreasing=TRUE)]
}

## calculate modularity
getMods <- function(null){
    out <- lapply(null,computeModules)
    out <- unlist(lapply(out,function(x) 
        slot(x,'likelihood')))
    return(out)
}
    
## calculate contribution to modularity
getCM <- function(x,null){
    out <- list()
    for (i in 1:nrow(x)){
        m <- x
        mod <- numeric()
        for (j in 1:length(null)){
            m[,i] <- null[[j]][,i]
            mod[j] <- slot(computeModules(m),
                        "likelihood")
        }
        out[[i]] <- mod
    }
    return(out)
}
    
## calculate statistics
getStats <- function(obs,sim){
    mu.sim <- mean(sim)
    sd.sim <- sd(sim)
    z <- (obs - mu.sim) / sd.sim
    p <- length(sim[sim >= obs]) / length(sim)
    out <- c(obs,mu.sim,sd.sim,z,p)
    return(out)
}

## Multiple plot function
##
## ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {


                                        # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

                                        # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
                                        # Make the panel
                                        # ncol: Number of columns of plots
                                        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {
                                        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                                        # Make each plot, in the correct location
        for (i in 1:numPlots) {
                                        # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
        }
    }
}

                                        # creates a new data frame from an existing one containing repeated
                                        # measures as columns, e.g. a data frame named exp1.df:
                                        # subjectagerepeat1repeat2
                                        # 001341.451.67
                                        # 002381.201.54
                                        # ...
                                        # called like:
                                        # make.rm(1:2,3:4,exp1.df)
                                        # would multiply the "constant" variables by the number of
                                        # repeats and reformat the repeats to a single column.
                                        # subjectagerepdatcontrasts
                                        # 001341.45T1
                                        # 002381.20T1
                                        # ...
                                        # 001341.67T2
                                        # 002381.54T2
                                        # ...
                                        # this allows a "univariate" repeated measures analysis of the data.

# creates a new data frame from an existing one containing repeated
# measures as columns, e.g. a data frame named exp1.df:
# subject  age	       repeat1	   repeat2
# 001	   	       34	   1.45	1.67
# 002		       		   38	1.20	1.54
# ...
# called like:
# make.rm(1:2,3:4,exp1.df)
# would multiply the "constant" variables by the number of
# repeats and reformat the repeats to a single column.
# subject age repdat   contrasts
# 001	      34       1.45	T1
# 002	      	       38	1.20	T1
# ...
# 001		34	1.67	T2
# 002			38	1.54	T2
# ...
# this allows a "univariate" repeated measures analysis of the data.

make.rm<-function(constant,repeated,data,contrasts) {
    if(!missing(constant) && is.vector(constant)) {
        if(!missing(repeated) && is.vector(repeated)) {
            if(!missing(data)) {
                dd<-dim(data)
                replen<-length(repeated)
                if(missing(contrasts))
                    contrasts<-
                        ordered(sapply(paste("T",1:length(repeated),sep=""),rep,dd[1]))
                else
                    contrasts<-matrix(sapply(contrasts,rep,dd[1]),ncol=dim(contrasts)[2])
                if(length(constant) == 1)
                    cons.col<-rep(data[,constant],replen)
                else cons.col<-lapply(data[,constant],rep,replen)
                new.df<-data.frame(cons.col,
                                   repdat=as.vector(data.matrix(data[,repeated])),
                                   contrasts)
                return(new.df)
            }
        }
    }
    cat("Usage: make.rm(constant, repeated,
 data [, contrasts])\n")
    cat("\tWhere 'constant' is a vector of indices of non-repeated data
 and\n")
    cat("\t'repeated' is a vector of indices of the repeated measures
 data.\n")
}
