library(bipartite)
library(RCurl)
library(XML)
library(vegan)
library(ggnet)

unipart <- function(x,rm.zero=TRUE,plot=TRUE){
    out <- t(as.matrix(sign(x))) %*% as.matrix(sign(x))
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

chPlot <- function(x,f,col,pch,se=FALSE,xlim=c(-1,1),ylim=c(-1,1),xlab='X',ylab='Y',add.plot=FALSE,line.lm=FALSE,line.lty=1,line.col=1){
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
        points(mu,col=col,pch=pch)
    }else{
        plot(mu,col=col,pch=pch,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    }
    for (i in 1:nrow(mu)){
        lines(c(bar.up[i,1],bar.lo[i,1]),rep(mu[i,2],2),col=col[i])
        lines(rep(mu[i,1],2),c(bar.up[i,2],bar.lo[i,2]),col=col[i])
    }
    if (line.lm){
        abline(lm(mu[,2]~mu[,1]),lty=line.lty,col=line.col)
    }
    return(col)
}

ch.point <- function(x,f){
    list(mu=tapply(x,f,mean),se=tapply(x,f,function(x) sd(x)/sqrt(length(x))))
}

ch.vector <- function(x,t,pch,col){
    mu <- apply(x,2,function(x,t) tapply(x,t,mean),t=t)
    se <- apply(x,2,function(x,t) tapply(x,t,function(x) sd(x)/sqrt(length(x))),t=t)
    cbind(mu,se)
}

ch.plot <-
function(x='ordination matrix',g='groupings',cex=1,plot.legend=FALSE,loc='topleft',mu.pch=19){
  mu <- apply(x,2,function(x,g) tapply(x,g,mean),g=g)
  se <- apply(x,2,function(x,g) tapply(x,g,function(x) sd(x)/sqrt(length(x))),g=g)
  mu <- na.omit(mu)
  se <- na.omit(se)
                                        #error bars
  cl.xu <- mu[,1] + se[,1]
  cl.xl <- mu[,1] - se[,1]
  cl.yu <- mu[,2] + se[,2]
  cl.yl <- mu[,2] - se[,2]
    if (plot.legend){
                                        #coloring
      mu.col <- rainbow(length(unique(g)))[as.numeric(unique(g))]
      plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
      for (i in 1:nrow(mu)){
        lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
        lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
      }    
      legend(loc,legend=rownames(se),cex=cex*0.5,pch=mu.pch,col=mu.col,border='grey')
  }else{
                                        #coloring
    mu.col <- 'black'
    plot(mu,pch=mu.pch,cex=cex,xlim=c(min(cl.xl),max(cl.xu)),ylim=c(min(cl.yl),max(cl.yu)),col=mu.col)
    for (i in 1:nrow(mu)){
      lines(x=c(cl.xl[i],cl.xu[i]),y=c(mu[i,2],mu[i,2]))
      lines(x=c(mu[i,1],mu[i,1]),y=c(cl.yl[i],cl.yu[i]))
    }
  }
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
