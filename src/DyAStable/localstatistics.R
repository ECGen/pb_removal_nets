library(igraph)
library(proxy)

"+" <- function(...) UseMethod("+") 
"+.default" <- .Primitive("+") 
"+.character" <- function(...) paste(...,sep="") 
assign.values.to.nodes <- function(graphs,years=NULL,func="betweenness")
{
tryCatch({
	for(i in 1:length(graphs))
	{
		if(is.null(years) || length(years[years==graphs[[i]]$year]) != 0)
		{
			graph <- list()
			graph[[1]] <- graphs[[i]]
			for(f in func)
			{
				res <- do.call(f,graph)
				graphs[[i]] <- set.vertex.attribute(graphs[[i]],f,value=res)
			}
		}
	}
	graphs
	},error=function(ex){print("error")})
}

eigen <- function(graph)
{
tryCatch({
	res <- evcent(graph)
	res$vector
},error=function(ex){print("error")})	
}

pr.cent <- function(graph)
{
tryCatch({
	res <- page.rank(graph)
	res$vector
},error=function(ex){print("error")})
}

kleinberg <- function(graph)
{
tryCatch({
	res <- hub.score(graph)
	res$vector
},error=function(ex){print("error")})
}
# !!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!!! #

set.degree.attributes <- function(graphs,years=NULL,mode=c("all","out","in","total"),loops=TRUE)
{
tryCatch({
	for(i in 1:length(graphs))
	{
		if(is.null(years) || length(years[years==graphs[[i]]$year]) != 0)
		{
			graphs[[i]] <- set.vertex.attribute(graphs[[i]],"degree",value=degree(graphs[[i]],mode=mode,loops=loops))
		}
	}
	graphs
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

cut.irrelevant.nodes <- function(graphs,years=NULL,attr,th,cutmode="fix")
{
	tryCatch({
	if(length(attr) > length(th))
	{
		stop("not enough threshold value")
	}
	
	newgraphs <- NULL
	left <- NULL
	for(i in 1:length(graphs))
	{
		
		if(is.null(years) || length(years[years==graphs[[i]]$year]) != 0)
		{
			k <- 0
			newgraphs[[i]] <- graphs[[i]]
			for(a in attr)
			{
				attributes <- get.vertex.attribute(newgraphs[[i]],a)
				if(is.null(attributes))stop("this attribute is not set yet")
				if(cutmode == "fix")
				{
					k <- k + 1
					relevantNodes <- V(newgraphs[[i]])[attributes>th[k]]
				
				}else if(cutmode == "rate")
				{
					k <- k + 1
					relevantNodes <- V(newgraphs[[i]])[attributes>max(attributes) / th[k]]
				}
				else
				{
					stop("invalid cutmode parameter")
				}
				newgraphs[[i]] <- subgraph(newgraphs[[i]],relevantNodes)
			}
			nodes <-  length(V(graphs[[i]])) - length(V(newgraphs[[i]]))
			edges <- length(E(graphs[[i]])) - length(E(newgraphs[[i]]))
			stats <- NULL
			stats$nodes <- nodes
			stats$edges <- edges
			left[[i]] <- stats
			graph <- list()
			graph[[1]] <- newgraphs[[i]]
			for(a in attr)
			{
				newgraphs[[i]] <- set.vertex.attribute(newgraphs[[i]],a,value=do.call(a,graph))
			}
		}else 
		{
			newgraphs[[i]] <- graphs[[i]]
		}
	}
	result <- NULL
	result$graph <- newgraphs
	result$left <- left
	result
	},error=function(ex){print("error")})
}


# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

cut.irrelevant.attributes <- function(graphs,years=NULL,attr,th,cutmode="fix")
{
tryCatch({
	if(length(attr) > length(th))stop("not enough threshold value")
	for(i in 1:length(graphs))
	{
		if(is.null(years) || length(years[years==graphs[[i]]$year]) != 0)
		{
			k <- 0
			for(a in attr)
			{
				
				attributes <- get.vertex.attribute(graphs[[i]],a)
				if(is.null(attributes))stop("this attribute is not set yet")
				if(cutmode == "fix")
				{
					k <- k + 1
					graphs[[i]] <- set.vertex.attribute(graphs[[i]],a,index=V(graphs[[i]])[attributes<th[k]],value=0)
				}
				else if(cutmode == "rate")
				{
					k <- k + 1
					graphs[[i]] <- set.vertex.attribute(graphs[[i]],a,index=V(graphs[[i]])[attributes<max(attributes)/th[k]],value=0)
				}
				else
				{
					stop("invalid cutmode parameter")
				}
			}
		}
	}
	graphs
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

dynamic.node.statistics <- function(graphs,years=NULL,attr,th=NULL,cutmode="rate")
{
tryCatch({
	cut <- TRUE
	if(!is.null(th) && cutmode=="rate" && th<=1)stop("threshold value must be bigger than 1")
	attrtable <- NULL
	ys <- NULL
	d <- FALSE
	for (graph in graphs) 
	{
		if(!is.igraph(graph))stop("not a graph object")
		if(is.null(years) || length(intersect(years,graph$interval)) != 0)
		{
			ys <- c(ys,graph$year)
			attribute <- get.vertex.attribute(graph,attr)
			if(is.null(attribute)) stop("this attribute is not set yet")
			if(is.null(th))
			{
				if(cut)
				{
					relevantNodes <- V(graph)[attribute>summary(attribute)[[4]]]
					newgraph <- subgraph(graph,relevantNodes)
					bp_filtered=attribute[attribute>summary(attribute)[[4]]]
					names(bp_filtered) <- V(newgraph)$id
				}
				else
				{
					names(attribute) <- V(graph)$id
					bp_filtered=attribute[attribute>summary(attribute)[[4]]]
				}
			}
			else
			{
				if(cut)
				{
					if(cutmode == "rate")
					{
						names(attribute) <- V(graph)$id
						relevantNodes <- V(graph)[attribute>(max(attribute)/th)]
						newgraph <- subgraph(graph,relevantNodes)
						bp_filtered <- attribute[attribute>(max(attribute)/th)]
						
					}else if(cutmode == "fix")
					{
						relevantNodes <- V(graph)[attribute>th]
						newgraph <- subgraph(graph,relevantNodes)
						bp_filtered=attribute[attribute>th]
						names(bp_filtered) <- V(newgraph)$id
					}
					else
					{
						stop("invalid cutmode")
					}
				}
				else
				{
					if(cutmode == "rate")
					{
						names(attribute) <- V(graph)$id
						bp_filtered <- attribute[attribute>(max(attribute)/th)]
					}else if(cutmode == "fix")
					{
						names(attribute) <- V(graph)$id
						bp_filtered=attribute[attribute>th]
					}
				}
			}
			
			if(length(bp_filtered) != 0)
			{
				d <- TRUE
				bp_filtered_sorted=sort(bp_filtered,decreasing=TRUE);
				uj=data.frame(names(bp_filtered_sorted),bp_filtered_sorted,graph$year)
				attrtable=rbind(attrtable,uj)
			}
				
		}
	}
	centtab <- NULL
	if(d)centtab=xtabs(bp_filtered_sorted~names.bp_filtered_sorted.+graph.year,data=attrtable)
	centtab
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

plot.dynamic.node.statistics <- function(graphs,years=NULL,attr,th=NULL,cutmode="rate",ranked=TRUE,plotmode="hist",log=FALSE)
{
tryCatch({	
	attrtable <- dynamic.node.statistics(graphs,years,attr,th,cutmode)	
	years <- as.numeric(colnames(attrtable))
	allys <- NULL
	c <- 1
	periods <- NULL
	lastyears <- NULL
	a <- 0
	for(g in graphs)
	{
		allys <- c(allys,g$year)
	}	
	l <- NULL
	newColor <- "green"
	colors <- data.frame(year=years,color = sample(colors(),size=length(years),replace=FALSE)) 
	oldVertices <- NULL
	actualColors <- NULL
	newVertex <- NULL
	weights <- NULL
	filenames <- NULL
	n <- 0
	f <- 0
	periods <- NULL
	for(y in years)
	{
		period <- NULL
		for(g in graphs)
		{
			if(any(g$interval == y))
			{
				period$l <- g$label
				period$i <- g$interval
				break
			}
		}
		periods[[f <- f+1]] <- period
		verticesInYear <- NULL
		centFrame <- data.frame(attrtable)
		centFrame <- data.frame(content = centFrame$names.bp_filtered_sorted.,year=centFrame$graph.year,cent=centFrame$Freq)
		centInYear <- subset(centFrame,year==y,select=c(content,cent))
		centNonZero <- subset(centInYear,cent>0,select=c(content,cent))
		actualColors <- NULL
		centOrderedNodes <- centNonZero[sort.list(centNonZero$cent),]
		if(ranked)
		{
			terms <- centOrderedNodes
		}
		else
		{
			terms <- centNonZero
		}
		for(j in 1:length(terms$content)){
			if(!is.null(oldVertices))
				newVertex <- subset(oldVertices,content==terms$content[j],select=content)
			if(length(newVertex$content) == 0){
				oldVertices <- rbind(oldVertices,data.frame(content=terms$content[j],year=y,age=1,weight=terms$cent[j]))
				actualColors <- c(actualColors,newColor)
			}else{
				Id <- which(oldVertices$content==terms$content[j])
				oldVertices$age[Id] = oldVertices$age[Id] + 1
				oldVertices$weight[Id] = oldVertices$weight[Id] + oldVertices$age[Id]*terms$cent[j]
				oldyear <- as.numeric(subset(oldVertices,content == terms$content[j],select=year))
				col <- subset(colors,year==oldyear,select=color)
				actualColors <- c(actualColors,toString(col$color))
				verticesInYear <- c(verticesInYear,oldyear)  
			}
			ordinate <- terms$content
			absciss <- terms$cent  
		}
		if(plotmode=="hist")
		{
			if(ranked)
			{
				filename <- attr + "Dist_" + period$l + "_" +unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height=absciss,names.arg=ordinate,horiz=TRUE,las=2,cex.names=1,cex.axis=0.5,col=actualColors,main=paste(paste(paste("distribution for the node attribute: ",attr)," in the year "), period$l))
			}
			else
			{
				
				filename <- attr + "Dist_" + period$l + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height=absciss,names.arg=ordinate,legend.text=verticesInYear,horiz=TRUE,las=2,cex.names=1,cex.axis=0.5,col=actualColors,main=paste(paste(paste("distribution for the node attribute: ",attr)," in the year "), period$l))
			}
		}else if(plotmode == "point")
		{
			if(log)
			{
				filename <- attr + "_log_Dist_" + period$l + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				plot(x=1:length(absciss),y=rev(absciss), type='p',log="xy",main=paste(paste(paste("log distribution for the node attribute: ",attr)," in the year "), period$l),ylab=paste(attr,"scores"))
			}
			else
			{
				filename <- attr + "Dist_" + period$l + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(family="serif",mfrow=c(1,1))
				plot(x=1:length(absciss),y=rev(absciss),type='p',main=paste(paste(paste("distribution for the node attribute: ",attr)," in the year "), period$l),xlab="nodes",ylab=paste(attr,"scores"))	
			}
		}
		else
		{
			stop("invalid plotmode")
		}
		dev.off()
		#filenames <- c(filenames,filename)
		t <- NULL
		t$year <- y
		t$filename <- filename
		l[[n <- n+1]] <- t
	}
	if(is.null(l))return(rep("NaP",length(allys)))
	i <- 1
	for(y in allys)
	{
		if(y == l[[i]]$year)
		{
			filenames <- c(filenames,l[[i]]$filename)
			i <- i+1
		}
		else if(y > l[[i]]$year)
		{
			
		}
		else if(y < l[[i]]$year)
		{
			filenames <- c(filenames, "NaP")
		}
	}
	paste(paste(getwd(), "/",sep=""), filenames,sep="")
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

vertex.betweenness.statistics <- function(graphs,years=NULL,th=NULL,cutmode="rate")
{
tryCatch({
	vertex.centrality.statistics(graphs,years,method="betweenness",th,cutmode)
},error=function(ex){print("error")})
}
	
# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.betweenness.statistics <- function(graphs,years=NULL,th=8,ranked=TRUE,cutmode="rate",plotmode="hist",log=FALSE)
{
tryCatch({
	plot.centrality.statistics(graphs,years,method="betweenness",th,ranked,cutmode,plotmode,log)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	cdn.statistics <- function(tab,th=5)
	{
	tryCatch({
		if(is.null(tab))return(NULL)
		years <- as.numeric(colnames(tab))
		oldVertices <- NULL
		newVertex <- NULL
		weights <- NULL
		centralities <- NULL
		bwFrame <- data.frame(tab)
		bwFrame <- data.frame(content = bwFrame$names.bp_filtered_sorted.,year=bwFrame$graph.year,bw=bwFrame$Freq)
		for(y in years)
		{
			verticesInYear <- NULL
			bwInYear <- subset(bwFrame,year==y,select=c(content,bw))
			bwNonZero <- subset(bwInYear,bw>0,select=c(content,bw))
			for(j in 1:length(bwNonZero$content)){
				if(!is.null(oldVertices))
					newVertex <- subset(oldVertices,content==bwNonZero$content[j],select=content)
				if(length(newVertex$content) == 0){
					oldVertices <- rbind(oldVertices,data.frame(content=bwNonZero$content[j],year=y,age=1,weight=bwNonZero$bw[j]))
				}else{
					Id <- which(oldVertices$content==bwNonZero$content[j])
					oldVertices$age[Id] = oldVertices$age[Id] + 1
					oldVertices$weight[Id] = oldVertices$weight[Id] + oldVertices$age[Id]*bwNonZero$bw[j]
					oldyear <- as.numeric(subset(oldVertices,content == bwNonZero$content[j],select=year))
					verticesInYear <- c(verticesInYear,oldyear)  
				}
			}
		}
		winner <- max(oldVertices$weight)
		centralities[[1]] <- subset(oldVertices,weight>winner/th)
		centralities[[2]] <- years
		centralities
	},error=function(ex){print("error")})
	}
	
# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	betweenness.cdn.statistics <- function(graphs,years=NULL,staticth=8,dynamicth=5,cutmode="rate")
	{
	tryCatch({
		if(staticth <= 1) stop("static threshold must be bigger than 1")
		bwtab <- vertex.betweenness.statistics(graphs,years,th=staticth,cutmode=cutmode)
		cdn.statistics(bwtab,dynamicth)
	},error=function(ex){print("error")})
	}

	# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	cdn.betweenness.stability <- function(graphs,years=NULL,staticth=8,dynamicth=5,cutmode="rate")
	{
	tryCatch({
		if(staticth <= 1) stop("static threshold must be bigger than 1")
		centralities <- betweenness.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
		get.cdn.stability(graphs,centralities)
	},error=function(ex){print("error")})
	}	

	# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	get.cdn.stability <- function(graphs,centralities)
	{
	tryCatch({
		stability.scores <- NULL
		years <-  centralities[[2]]
		centralities <- centralities[[1]]
		allNeighbors <- NULL
		g <- 0
		for(graph in graphs)
		{
			if(is.null(years) || length(years[years==graph$year]) != 0)
			{
				for(c in 1:(length(centralities$content)+1))
				{ 
					id <- which(V(graph)$id==centralities$content[c])
					if(length(id) != 0){
						neighbors <- neighbors(graph,id-1)
						neighborsText <- NULL
						weights <-  NULL
						for(i in 1:length(neighbors))
						{
							neighborsText <- c(neighborsText,V(graph)$id[neighbors[i]+1])
							weights <- c(weights,E(graph,P=c(id-1,neighbors[i]))$weight)
						}	
						cent <- NULL
						cent[1:length(neighbors)] <- toString(centralities$content[c])
						allNeighbors <- rbind(allNeighbors,data.frame(cent,neighborsText,graph$year,weights))
					}     
				}
			}	
		}
		stability <- NULL
		jaccards <- NULL
		jaccards2 <- NULL
		i <- 0
		for(c in 1:length(centralities$content))
		{	
			contentCounter <- NULL
			filteredByContent <- subset(allNeighbors,cent == toString(centralities$content[c]),select=c(neighborsText,graph.year,weights))	
			#contentTab <- table(filteredByContent$cent,filteredByContent$neighborsText)
			diversities <- NULL
			for(i in 1:length(filteredByContent$neighborsText))
			{
				item <- subset(contentCounter,contentCounter$content == toString(filteredByContent$neighborsText[i]))
				if(length(item$content)==0){
					contentCounter <- rbind(contentCounter,data.frame(content = toString(filteredByContent$neighborsText[i]),count = 1,weight=filteredByContent$weights[i],lastWeight=filteredByContent$weights[i],diversity=0,monotonity=0))
				}else{
					id <- which(contentCounter$content == toString(filteredByContent$neighborsText[i]))
					contentCounter$count[id] <- contentCounter$count[id] + 1
					contentCounter$weight[id] <- contentCounter$weight[id] + filteredByContent$weights[i]
					contentCounter$diversity[id] <- contentCounter$diversity[id] + abs(filteredByContent$weights[i] - contentCounter$lastWeight[id])
					if(contentCounter$monotonity[id] == 0){
						if(filteredByContent$weight[i] < contentCounter$lastWeight[id])
							contentCounter$monotonity[id] <- 1
						else
							contentCounter$monotonity[id] <- 2
					}
					if(contentCounter$monotonity[id] == 1)
					{
						if(filteredByContent$weight[i] > contentCounter$lastWeight[id])
							contentCounter$monotonity[id] <- 3
					}
					if(contentCounter$monotonity[id] == 2)
					{
						if(filteredByContent$weight[i] < contentCounter$lastWeight[id])
							contentCounter$monotonity[id] <- 3
					}
					contentCounter$lastWeight[id] <- filteredByContent$weights[i]           
				}	
			}
			contentCounter$weight <- contentCounter$weight / contentCounter$count
			contentCounter$diversity <- contentCounter$diversity / (contentCounter$count-1)
			longLifeNeighbors <- subset(contentCounter,count>1)
			longLifeNeighbors_sorted <- longLifeNeighbors[sort.list(longLifeNeighbors$weight),]
			
			if(length(longLifeNeighbors_sorted$content) != 0)
			{
				if(sum(longLifeNeighbors$diversity)==0)
				{
					stability = c(stability,100)
				}else
				{
					stability <- c(stability,sum(longLifeNeighbors$count)/(sum(longLifeNeighbors$diversity)/length(longLifeNeighbors$diversity)))
				}
				t <- table(filteredByContent$neighborsText,filteredByContent$graph.year)
				fr <- data.frame(t)
				jac <- 0
				jac2 <- 0
				for(ys in (length(years)-1))
				{
					fr1 <- subset(fr,fr[2] == years[ys])
					fr2 <- subset(fr,fr[2] == years[ys+1])
					if(!(length(fr1$Freq) == 0 || length(fr2$Freq) == 0)){
						freq1 <- fr1$Freq
						freq2 <- fr2$Freq
						inter <- freq1 & freq2
						union <- freq1 | freq2
						jaccard <- sum(inter) / sum(union)
						jac <- jac + jaccard
						cont <- fr2$Var1[inter]
						if(length(cont)!=0){
							divsum <- 0
							for(i in 1:length(cont)){
								filt <- subset(allNeighbors,cent==toString(centralities$content[c]) & neighborsText==cont[i] & (graph.year==years[ys] | graph.year==years[ys+1]))
								divsum <- divsum + abs(filt$weights[2] - filt$weights[1])
							}	
							jaccard2 <- (1-divsum)/sum(union)
							jac2 <- jac2 + jaccard2
						}	
					}       
				}
			
				jaccards <- c(jaccards,jac)
				jaccards2 <- c(jaccards2,jac2)
			}
			else
			{
				jaccards <- c(jaccards,0)
				jaccards2 <- c(jaccards2,0)
				stability <- c(stability,0)
			}
		}	
		names <- as.vector(centralities$content)
		stability.scores$names <- names
		stability.scores$measure1 <- stability
		stability.scores$measure2 <- jaccards
		stability.scores$measure3 <- jaccards2
		stability.scores
	},error=function(ex){print("error")})
	}

	# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
vertex.closeness.statistics <- function(graphs,years=NULL,th=NULL,cutmode="rate")
{
	tryCatch({
	cltab <- vertex.centrality.statistics(graphs,years,method="closeness",th,cutmode)
	cltab
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

plot.closeness.statistics <-function(graphs,years=NULL,th=1.001,ranked=TRUE,cutmode="rate",plotmode="hist",log=FALSE)
{
	tryCatch({
	plot.centrality.statistics(graphs,years,method="closeness",th,ranked,cutmode,plotmode,log)
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

eigen.vector.statistics <- function(graphs,years=NULL,th=NULL,cutmode="rate")
{
	tryCatch({
	evcenttab <- vertex.centrality.statistics(graphs,years,method="eigen",th,cutmode)
	evcenttab
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.evcent.statistics <- function(graphs,years=NULL,th=1.001,ranked=TRUE,cutmode="rate",plotmode="hist",log=FALSE)
{
	tryCatch({
	plot.centrality.statistics(graphs,years,method="eigenvector",th,ranked,cutmode,plotmode,log)
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
kleinberg.statistics <- function(graphs,years=NULL,th=NULL,cutmode="rate")
{
	tryCatch({
	vertex.centrality.statistics(graphs,years,method="kleinberg",th,cutmode)
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.kb.statistics <- function(graphs,years=NULL, th=4,ranked=TRUE,cutmode="rate",plotmode="hist",log=FALSE)
{
	tryCatch({
	plot.centrality.statistics(graphs,years,method="kleinberg",th,ranked,cutmode,plotmode,log)
	},error=function(ex){print("error")})
}
	
# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
page.rank.statistics <- function(graphs,years=NULL,th=NULL,cutmode="rate")
{
tryCatch({
	vertex.centrality.statistics(graphs,years,method="pagerank",th,cutmode)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.pr.statistics <- function(graphs,years=NULL,th=1.05,ranked=TRUE,cutmode="rate",plotmode="hist",log=FALSE)
{	
tryCatch({
	plot.centrality.statistics(graphs,years,method="pagerank",th,ranked,cutmode,plotmode,log)
},error=function(ex){print("error")})
}	


# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

plot.betweenness.stability <- function(graphs,years=NULL,staticth=4,dynamicth=6,cutmode="rate",plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- betweenness.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
	plot.stability(graphs,centralities,plotneighbors,cdnplotmode,stabilitymeasure)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

closeness.cdn.statistics <- function(graphs,years=NULL,staticth=1.001,dynamicth=4,cutmode="rate")
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	cltab <- vertex.centrality.statistics(graphs,years,method="closeness",th=staticth,cutmode=cutmode)
	cdn.statistics(cltab,dynamicth)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #

cdn.closeness.stability <- function(graphs,years=NULL,staticth=1.001,dynamicth=4,cutmode="rate")
{	
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- closeness.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
	get.cdn.stability(graphs,centralities)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.closeness.stability <- function(graphs,years=NULL,staticth=1.001,dynamicth=4,cutmode="rate",plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- closeness.cdn.statistics(graphs,years,staticth,dynamicth,cutmode) 
	plot.stability(graphs,centralities,plotneighbors,cdnplotmode,stabilitymeasure)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
evcent.cdn.statistics <- function(graphs,years=NULL,staticth=1.0001,dynamicth=2.5,cutmode="rate")
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	evcenttab <- vertex.centrality.statistics(graphs,years,method="eigenvector",th=staticth,cutmode=cutmode)
	cdn.statistics(evcenttab,dynamicth)
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
cdn.evcent.stability <- function(graphs,years=NULL, staticth=1.0001,dynamicth=2.5,cutmode="rate")
{	
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- evcent.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
	get.cdn.stability(graphs,centralities)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.evcent.stability <- function(graphs,years=NULL,staticth=1.001,dynamicth=2.5,cutmode="rate",plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- evcent.cdn.statistics(graphs,years,staticth,dynamicth,cutmode) 
	plot.stability(graphs,centralities)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
kleinberg.cdn.statistics <- function(graphs,years=NULL,staticth=2,dynamicth=3,cutmode="rate")
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	kbtab <- vertex.centrality.statistics(graphs,years,method="kleinberg",th=staticth,cutmode=cutmode)
	cdn.statistics(kbtab,dynamicth)		
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
cdn.kleinberg.stability <- function(graphs,years=NULL,staticth=3,dynamicth=3,cutmode="rate")
{	
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- kleinberg.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
	get.cdn.stability(graphs,centralities)
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.kleinberg.stability <- function(graphs,years=NULL,staticth=3,dynamicth=3,cutmode="rate",plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- kleinberg.cdn.statistics(graphs,years,staticth,dynamicth,cutmode) 
	plot.stability(graphs,centralities,plotneighbors,cdnplotmode,stabilitymeasure)		
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
page.rank.cdn.statistics <- function(graphs,years=NULL,staticth=1.07,dynamicth=1.5,cutmode="rate")
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	prtab <- vertex.centrality.statistics(graphs,years,method="pagerank",th=staticth,cutmode=cutmode)
	cdn.statistics(prtab,dynamicth)	
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
cdn.page.rank.stability <- function(graphs,years=NULL,staticth=1.07,dynamicth=1.5,cutmode="rate")
{	
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- page.rank.cdn.statistics(graphs,years,staticth,dynamicth,cutmode)
	get.cdn.stability(graphs,centralities)	
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.page.rank.stability <- function(graphs,years=NULL,staticth=1.07,dynamicth=1.5,cutmode="rate",plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(staticth <= 1 && cutmode=="rate") stop("static threshold must be bigger than 1")
	if(dynamicth <= 1 && cutmode == "rate") stop("dynamic threshold must be bigger than 1")
	centralities <- page.rank.cdn.statistics(graphs,years,staticth,dynamicth,cutmode) 
	plot.stability(graphs,centralities,plotneighbors,cdnplotmode,stabilitymeasure)	
},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
		
vertex.centrality.statistics <- function(graphs,years=NULL,method="betweenness",th=NULL,cutmode="rate")
{
tryCatch({
	t.cent <- NULL
	ys <- NULL
	cut <- TRUE
	d <- FALSE
	for (graph in graphs) 
	{
		if(!is.igraph(graph))stop("not a graph object")
		if(is.null(years) || length(years[years==graph$year]) != 0)
		{
			if(method=="betweenness")
			{
				cent <- betweenness(graph,v=V(graph),directed = FALSE)
				n1 <- set.vertex.attribute(graph,"bw",index=V(graph),cent)
			}
			else if(method == "closeness")
			{
				cent <- closeness(graph)
				n1 <- set.vertex.attribute(graph,"cl",index=V(graph),cent)
			}
			else if(method == "eigenvector")
			{
				cent <- evcent(graph)$vector
				n1 <- set.vertex.attribute(graph,"eigen",index=V(graph),cent)
			}
			else if(method == "kleinberg")
			{
				cent <- hub.score(graph)$vector
				n1=set.vertex.attribute(graph,"kb",index=V(graph),cent)
			}
			else if(method == "pagerank")
			{
				cent <- page.rank(graph)$vector
				n1=set.vertex.attribute(graph,"pr",index=V(graph),cent)
			}
			else
			{
				stop("invalid method")
			}
			ys <- c(ys,graph$year)
			
			if(is.null(th))
			{
				if(cut)
				{
					relevantNodes <- V(graph)[cent>summary(cent)[[4]]]
					newgraph <- subgraph(graph,relevantNodes)
					bp_filtered=cent[cent>summary(cent)[[4]]]
					names(bp_filtered) <- V(newgraph)$id
				}
				else
				{
					names(cent) <- V(graph)$id
					bp_filtered=cent[cent>summary(cent)[[4]]]
				}
			}
			else
			{
				if(cut)
				{
					if(cutmode == "rate")
					{
						names(cent) <- V(graph)$id
						relevantNodes <- V(graph)[cent>(max(cent)/th)]
						newgraph <- subgraph(graph,relevantNodes)
						bp_filtered <- cent[cent>(max(cent)/th)]
						
					}else if(cutmode == "fix")
					{
						relevantNodes <- V(graph)[cent>th]
						newgraph <- subgraph(graph,relevantNodes)
						bp_filtered=cent[cent>th]
						names(bp_filtered) <- V(newgraph)$id
					}
					else
					{
						stop("invalid cut mode")
					}
				}
				else
				{
					if(cutmode == "rate")
					{
						names(cent) <- V(graph)$id
						bp_filtered <- cent[cent>(max(cent)/th)]
					}else if(cutmode == "fix")
					{
						names(cent) <- V(graph)$id
						bp_filtered=cent[cent>th]
					}
					else
					{
						stop("invalid cut mode")
					}
				}
			}
			if(length(bp_filtered) != 0)
			{
				d <- TRUE
				bp_filtered_sorted=sort(bp_filtered,decreasing=TRUE);
				uj=data.frame(names(bp_filtered_sorted),bp_filtered_sorted,graph$year)
				t.cent=rbind(t.cent,uj)
			}
				
		}
	}
	centtab <- NULL
	if(d)centtab=xtabs(bp_filtered_sorted~names.bp_filtered_sorted.+graph.year,data=t.cent)
	centtab
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
plot.centrality.statistics <- function(graphs,years=NULL,method="betweenness",th=NULL,ranked=FALSE,cutmode="rate",plotmode="hist",log=FALSE)
{
tryCatch({
	centtab <- vertex.centrality.statistics(graphs,years,method,th,cutmode)
	par(mar=c(3,8,1,1),family="serif",mfrow=c(2,2))
	years <- as.numeric(colnames(centtab))
	#years <- colnames(centtab)
	allys <- NULL
	for(g in graphs)
	{
		allys <- c(allys,g$year)
	}
	l <- list()
	n <- 0
	newColor <- "green"
	colors <- data.frame(year=years,color = sample(colors(),size=length(years),replace=FALSE)) 
	oldVertices <- NULL
	actualColors <- NULL
	newVertex <- NULL
	weights <- NULL
	filenames <- NULL
	centFrame <- data.frame(centtab)
	centFrame <- data.frame(content = centFrame$names.bp_filtered_sorted.,year=centFrame$graph.year,cent=centFrame$Freq)
	for(y in years)
	{
		verticesInYear <- NULL
		centInYear <- subset(centFrame,year==y,select=c(content,cent))
		centNonZero <- subset(centInYear,cent>0,select=c(content,cent))
		actualColors <- NULL
		centOrderedNodes <- centNonZero[sort.list(centNonZero$cent),]
		if(ranked)
		{
			terms <- centOrderedNodes
		}
		else
		{
			terms <- centNonZero
		}
		for(j in 1:length(terms$content)){
			if(!is.null(oldVertices))
				newVertex <- subset(oldVertices,content==terms$content[j],select=content)
			if(length(newVertex$content) == 0){
				oldVertices <- rbind(oldVertices,data.frame(content=terms$content[j],year=y,age=1,weight=terms$cent[j]))
				actualColors <- c(actualColors,newColor)
			}else{
				Id <- which(oldVertices$content==terms$content[j])
				oldVertices$age[Id] = oldVertices$age[Id] + 1
				oldVertices$weight[Id] = oldVertices$weight[Id] + oldVertices$age[Id]*terms$cent[j]
				oldyear <- as.numeric(subset(oldVertices,content == terms$content[j],select=year))
				#oldyear <- subset(oldVertices,content == terms$content[j],select=year)
				col <- subset(colors,year==toString(oldyear),select=color)
				actualColors <- c(actualColors,toString(col$color))
				verticesInYear <- c(verticesInYear,oldyear)  
			}
			ordinate <- terms$content
			absciss <- terms$cent  
		}
		if(plotmode=="hist")
		{
			if(ranked)
			{
				filename <- method + "Dist_" + y + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height=absciss,names.arg=ordinate,horiz=TRUE,las=2,cex.names=1,cex.axis=0.5,col=actualColors,main=paste(method,y,sep=" "))
			}
			else
			{
				filename <- method + "Dist_" + y + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height=absciss,names.arg=ordinate,legend.text=verticesInYear,horiz=TRUE,las=2,cex.names=1,cex.axis=0.5,col=actualColors,main=paste(method,y,sep=" "))
			}
		}else if(plotmode == "point")
		{
			if(log)
			{
				filename <- method + "_log_Dist_" + y + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				plot(x=1:length(absciss),y=rev(absciss), type='p',log="xy",main=paste(method,paste(" log distribution ",y)),ylab=paste(method," scores"))
			}
			else
			{
				filename <- method + "Dist_" + y + "_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(family="serif",mfrow=c(1,1))
				plot(x=1:length(absciss),y=rev(absciss),type='p',main=paste(method,paste(" distribution ",y)),xlab="nodes",ylab=paste(method," scores"))	
			}
		}
		else
		{
			stop("invalid plotmode")
		}
		dev.off()
		t <- NULL
		t$year <- y
		t$filename <- filename
		l[[n <- n+1]] <- t
	}
	
	if(is.null(l))return(rep("NaP",length(allys)))
	i <- 1
	for(y in allys)
	{
		if(y == l[[i]]$year)
		{
			filenames <- c(filenames,l[[i]]$filename)
			i <- i+1
		}
		else if(y > l[[i]]$year)
		{
			
		}
		else if(y < l[[i]]$year)
		{
			filenames <- c(filenames, "NaP")
		}
	}
	paste(paste(getwd(), "/",sep=""), filenames,sep="")
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
		
plot.stability <- function(graphs,centralities,plotneighbors=TRUE,cdnplotmode=1,stabilitymeasure=2)
{
tryCatch({
	if(is.null(centralities))return("NaP")	
	filenames <- NULL
	filename <- NULL
	years <-  centralities[[2]]
	centralities <- centralities[[1]]
	allNeighbors <- NULL
	skeleton <- NULL
	for(graph in graphs)
	{
		if(is.null(years) || length(intersect(years,graph$interval)) != 0)
		{
			for(c in 1:(length(centralities$content)+1))
			{ 
				id <- which(V(graph)$id==centralities$content[c])
				if(length(id) != 0){
					neighbors <- neighbors(graph,id-1)
					if(length(neighbors) != 0)
					{
						neighborsText <- NULL
						weights <-  NULL
						for(i in 1:length(neighbors))
						{
							neighborsText <- c(neighborsText,V(graph)$id[neighbors[i]+1])
							weights <- c(weights,E(graph,P=c(id-1,neighbors[i]))$weight)
						}	
						cent <- NULL
						cent[1:length(neighbors)] <- toString(centralities$content[c])
						allNeighbors <- rbind(allNeighbors,data.frame(cent,neighborsText,graph$year,weights))
					}
				}
			}     
		}
	}
	
	stability <- NULL
	jaccards <- NULL
	jaccards2 <- NULL
	i <- 0

	for(c in 1:length(centralities$content))
	{
		contentCounter <- NULL
		print(allNeighbors)
		filteredByContent <- subset(allNeighbors,cent == toString(centralities$content[c]),select=c(neighborsText,graph.year,weights))	
		#contentTab <- table(filteredByContent$cent,filteredByContent$neighborsText)
		diversities <- NULL
		if(length(centralities$content) == 0)
		{
			stability <- c(stability,0)
			jaccards <- c(jaccards, 0)
			jaccards2 <- c(jaccards2, 0)
		}
		for(i in 1:length(filteredByContent$neighborsText))
		{
			item <- subset(contentCounter,contentCounter$content == toString(filteredByContent$neighborsText[i]))
			if(length(item$content)==0){
				contentCounter <- rbind(contentCounter,data.frame(content = toString(filteredByContent$neighborsText[i]),count = 1,weight=filteredByContent$weights[i],lastWeight=filteredByContent$weights[i],diversity=0,monotonity=0))
			}else{
				id <- which(contentCounter$content == toString(filteredByContent$neighborsText[i]))
				contentCounter$count[id] <- contentCounter$count[id] + 1
				contentCounter$weight[id] <- contentCounter$weight[id] + filteredByContent$weights[i]
				contentCounter$diversity[id] <- contentCounter$diversity[id] + abs(filteredByContent$weights[i] - contentCounter$lastWeight[id])
				if(contentCounter$monotonity[id] == 0){
					if(filteredByContent$weight[i] < contentCounter$lastWeight[id])
						contentCounter$monotonity[id] <- 1
					else
						contentCounter$monotonity[id] <- 2
				}
				if(contentCounter$monotonity[id] == 1)
				{
					if(filteredByContent$weight[i] > contentCounter$lastWeight[id])
						contentCounter$monotonity[id] <- 3
				}
				if(contentCounter$monotonity[id] == 2)
				{
					if(filteredByContent$weight[i] < contentCounter$lastWeight[id])
						contentCounter$monotonity[id] <- 3
				}
				contentCounter$lastWeight[id] <- filteredByContent$weights[i]           
			}	
		}
		contentCounter$weight <- contentCounter$weight / contentCounter$count
		contentCounter$diversity <- contentCounter$diversity / (contentCounter$count-1)
		longLifeNeighbors <- subset(contentCounter,count>1)
		longLifeNeighbors_sorted <- longLifeNeighbors[sort.list(longLifeNeighbors$weight),]
		if(length(longLifeNeighbors$count) !=0)
		{
		if(plotneighbors)
		{
			ncount <- length(longLifeNeighbors$count)
			level <- 255/ncount*2
			if(level > 30) level <- 30
			if(any(cdnplotmode == 1))
			{
				w2 <- longLifeNeighbors_sorted$weight
				colors <- NULL
				w3 <- NULL
				col <- 255
				col2 <- 255
				tomb <- NULL
				contents <- NULL
				for(i in 1:length(w2)){
					e <- w3[w3 == w2[i]]
					if(length(e) == 0){
						if(col2 > level){
							col2 <- col2 - level
						}else{
							if(col > 20){
								col <- col - level
							}
						}
						colors <- c(colors,rgb(col,col2,0,maxColorValue=255))
						w3 <- c(w3,w2[i])
					}else{
						colors <- c(colors,rgb(col,col2,0,maxColorValue=255))
					}
				}
		
				w2 <- round(w2,digits=2)
				filename <- "cdn_stat_" + centralities$content[c] + "_mode1_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height = longLifeNeighbors_sorted$count,names.arg = longLifeNeighbors_sorted$content,horiz=TRUE,las=2,cex.names=1,cex.axis=0.5,main=paste("statistics of the neighbors of the CDN node: ",centralities$content[c]),col=colors,legend.text=w2)
				dev.off()
		}

			inter <- intersect(longLifeNeighbors_sorted$content,centralities$content)
			if(length(inter) != 0)
				skeleton <- rbind(skeleton,data.frame(vertex1=centralities$content[c],vertex2=inter))
			if(any(cdnplotmode==2) || any(cdnplotmode==3)){
				longLifeNeighbors_sorted <- longLifeNeighbors[sort.list(longLifeNeighbors$diversity),]
				w2 <- longLifeNeighbors_sorted$diversity
				colors <- NULL
				colors2 <- NULL
				w3 <- NULL
				col <- 255
				col2 <- 255
				tomb <- NULL
				w4 <- longLifeNeighbors_sorted$monotonity
				contents <- NULL
				monotonities <- NULL
				for(i in 1:length(w2)){
					if(w4[i] == 1){
						colors2 <- c(colors2,"lightblue")
						monotonities <- c(monotonities,"decreasing")
					}	
					if(w4[i] == 2){
						colors2 <- c(colors2,"tomato")
						monotonities <- c(monotonities, "increasing")
					}
					if(w4[i] == 3){
						colors2 <- c(colors2, "yellow4")
						monotonities <- c(monotonities,"not monotone")
					}
					e <- w3[w3 == w2[i]]
					if(length(e) == 0){
						if(col2 > level){
							col2 <- col2 - level
						}else{
							if(col > 20){
								col <- col - level
							}
						}
						if(col2 < 0) col2 <- 1
						if(col < 0) col <- 1
						colors <- c(colors,rgb(col,col2,0,maxColorValue=255))
						w3 <- c(w3,w2[i])
					}else{
						colors <- c(colors,rgb(col,col2,0,maxColorValue=255))
					}
				}
			}
			w2 <- round(w2,digits=2)
			if(any(cdnplotmode == 2))
			{
				filename <- "cdn_stat_" + centralities$content[c] + "_mode2_" + unclass(Sys.time()) + ".png"
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height = longLifeNeighbors_sorted$count,names.arg = longLifeNeighbors_sorted$content,horiz=TRUE,las=2,cex.names=1,cex.axis=1,main=paste("statistics of the neigbors of CDN node: ",centralities$content[c]),col=colors,legend.text=w2)
				dev.off()
			}
			if(any(cdnplotmode == 3))	
			{
				filename <- "cdn_stat" + centralities$content[c] + "_mode3_" + unclass(Sys.time()) + ".png" 
				png(width=800,height=600,filename=filename)
				par(mar=c(3,8,1,1),family="serif",mfrow=c(1,1))
				barplot(height = w2, names.arg = longLifeNeighbors_sorted$content,horiz = TRUE,main=paste("statistics of the neighbors of CDN node: ",centralities$content[c]),col = colors2,cex.names=1,las=2,cex.axis=1, legend.text = monotonities)
				dev.off()
			}
		}
			if(sum(longLifeNeighbors$diversity == 0))
				stability <- c(stability,100)
			else
				stability <- c(stability,sum(longLifeNeighbors$count)/(sum(longLifeNeighbors$diversity)/length(longLifeNeighbors$diversity)))
			t <- table(filteredByContent$neighborsText,filteredByContent$graph.year)
			fr <- data.frame(t)
			jac <- 0
			jac2 <- 0
			for(ys in (length(years)-1))
			{
				fr1 <- subset(fr,fr[2] == years[ys])
				fr2 <- subset(fr,fr[2] == years[ys+1])
				if(!(length(fr1$Freq) == 0 || length(fr2$Freq) == 0)){
					freq1 <- fr1$Freq
					freq2 <- fr2$Freq
					inter <- freq1 & freq2
					union <- freq1 | freq2
					jaccard <- sum(inter) / sum(union)
					jac <- jac + jaccard
					cont <- fr2$Var1[inter]
					if(length(cont)!=0){
						divsum <- 0
						for(i in 1:length(cont)){
							filt <- subset(allNeighbors,cent==toString(centralities$content[c]) & neighborsText==cont[i] & (graph.year==years[ys] | graph.year==years[ys+1]))
							divsum <- divsum + abs(filt$weights[2] - filt$weights[1])
						}	
						jaccard2 <- (1-divsum)/sum(union)
						jac2 <- jac2 + jaccard2
					}	
				}       
			}
			jaccards <- c(jaccards,jac)
			jaccards2 <- c(jaccards2,jac2)
			filenames <- c(filenames,filename)
		}else{
			jaccards <- c(jaccards,0)
			jaccards2 <- c(jaccards2,0)
			stability <- c(stability,0)
		}
		
	}

	names <- as.vector(centralities$content)
	if(!is.null(skeleton ))
		graphSkeleton <- graph.data.frame(skeleton,directed=FALSE)
	label <- c(as.vector(skeleton$vertex1),as.vector(skeleton$vertex2))
	labelu <- unique(label)
	if( !is.null(skeleton))
	{
		filename <- "dyn_skeleton_" + unclass(Sys.time()) + ".png"
		filenames <- c(filenames,filename)
		png(width=800,height=600,filename=filename)
		plot.igraph(graphSkeleton,vertex.size=60,vertex.label=labelu,vertex.color="green")
		dev.off()
	}
	#if(length(stability) > 50)
		#stop("too much CDN to plot the stability diagrams. Set a correct threshold")
	if(any(stabilitymeasure == 1))
	{
		filename <- "cdn_stab_m_1_" + unclass(Sys.time()) + ".png"
		filenames <- c(filenames,filename)
		png(width=800,height=600,filename=filename)
		barplot(height = stability, names.arg = names,horiz = TRUE,main="stability of the CDNs, measure 1",cex.names=1,las=2,cex.axis=1)
		dev.off()
	}
	
	if(any(stabilitymeasure == 2))
	{
		filename <- "cdn_stab_m_2_" + unclass(Sys.time()) + ".png"
		filenames <- c(filenames,filename)
		png(width=800,height=600,filename=filename)
		barplot(height = jaccards, names.arg = names, horiz = TRUE, main = "stability of the CDNs, measure 2",cex.names=1,las=2,cex.axis=1)
		dev.off()
	}
	if(any(stabilitymeasure == 3))
	{
		filename <- "cdn_stab_m_3_" + unclass(Sys.time()) + ".png"
		filenames <- c(filenames,filename)
		png(width=800,height=600,filename=filename)
		barplot(height = jaccards2, names.arg = names, horiz = TRUE, main = "stability of the CDNs, measure 3",cex.names=1,las=2,cex.axis=1)
		dev.off()
	}
	if(is.null(filenames))"NaP"
	else {
		 paste(paste(getwd(), "/",sep=""), filenames,sep="")
	}
	},error=function(ex){print("error")})
}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	timeplot.node.statistics <-function(graphs,years=NULL,attr = "betweenness", th=NULL,cutmode="rate")
	{
	tryCatch({
		cut <- TRUE
		t.bw=NULL
		g <- 0
		ys <- NULL
		d <- FALSE
		filename <- attr + "_timeplot_" + unclass(Sys.time()) + ".png"
		png(width=800,height=600,filename=filename)
		for (graph in graphs) 
		{
			if(!is.igraph(graph))stop("not a graph object")
			if(is.null(years) || length(years[years==graph$year]) != 0)
			{
				if(attr=="betweenness")
				{
					cent <- betweenness(graph,v=V(graph),directed = FALSE)
				}
				else if(attr == "closeness")
				{
					cent <- closeness(graph)
				}
				else if(attr == "eigen")
				{
					cent <- evcent(graph)$vector
				}
				else if(attr == "kleinberg")
				{
					cent <- hub.score(graph)$vector
				}
				else if(attr == "page.rank")
				{
					cent <- page.rank(graph)$vector
				}else
				{
					g <- list()
					g[[1]] <- graph
					cent <- do.call(attr,g)
				}
				
				n1 <- set.vertex.attribute(graph,"cent",index=V(graph),cent)
				if(is.null(th))
				{
					if(cut)
					{
						relevantNodes <- V(graph)[cent>summary(cent)[[4]]]
						newgraph <- subgraph(graph,relevantNodes)
						bp_filtered=cent[cent>summary(cent)[[4]]]
						names(bp_filtered) <- V(newgraph)$id
					}
					else
					{
						names(cent) <- V(graph)$id
						bp_filtered=cent[cent>summary(cent)[[4]]]
					}
				}
				else
				{
					if(cut)
					{
						if(cutmode == "rate")
						{
							names(cent) <- V(graph)$id
							relevantNodes <- V(graph)[cent>(max(cent)/th)]
							newgraph <- subgraph(graph,relevantNodes)
							bp_filtered <- cent[cent>(max(cent)/th)]
						}else if(cutmode == "fix")
						{
							relevantNodes <- V(graph)[cent>th]
							newgraph <- subgraph(graph,relevantNodes)
							bp_filtered=cent[cent>th]
							names(bp_filtered) <- V(newgraph)$id
						}
						else 
						{
							stop("invalid cutmode")
						}
					}
					else
					{
						if(cutmode == "rate")
						{
							names(cent) <- V(graph)$id
							bp_filtered <- cent[cent>(max(cent)/th)]
						}else if(cutmode == "fix")
						{
							names(cent) <- V(graph)$id
							bp_filtered=cent[cent>th]
						}
						else
						{
							stop("invalid cutmode")
						}
					}
				}
				if(length(bp_filtered) != 0)
				{
					d <- TRUE
					ys <- c(ys,graph$year)
					bp_filtered_sorted=sort(bp_filtered,decreasing=TRUE);
					uj=data.frame(names(bp_filtered_sorted),bp_filtered_sorted,graph$year)
					t.bw=rbind(t.bw,uj)
				}
			}
		}
		bwtab <- NULL
		if(d)bwtab=xtabs(bp_filtered_sorted~names.bp_filtered_sorted.+graph.year,data=t.bw)
		if(is.null(bwtab))
		{
			print("threshold too high. Nothing plotted")
			return("NaP")
		}else
		{
			years <- ys
			y <- years
			bwFrame <- data.frame(bwtab)
			bwFrame <- data.frame(content = bwFrame$names.bp_filtered_sorted.,year=bwFrame$graph.year,bw=bwFrame$Freq)
			maxBw <-subset(bwFrame,bw==max(bwFrame$bw))
			maxNodeDynamics <- subset(bwFrame,content==maxBw$content[1])
			plot(y,maxNodeDynamics$bw,type='o',col=sample(colors()),main=paste("dynamic statistics of the node attribute: ",attr),ylab=paste(attr, " scores"))
			text(as.numeric(as.vector(maxBw$year)),maxBw$bw+maxBw$bw/50,maxBw$content[1],cex=0.7)
			bwFrame <- subset(bwFrame,content!=maxBw$content[1])
			for(term in bwFrame$content)
			{
				nodeDynamics<-subset(bwFrame,content==term)
				bwFrame <- subset(bwFrame,content!=term)
				localMax<- subset(nodeDynamics,bw==max(nodeDynamics$bw))
				if(length(nodeDynamics$bw[nodeDynamics$bw>0])>0)
				{
					lines(y,nodeDynamics$bw,type='o',col=sample(colors()))
					if(localMax$bw>maxBw$bw/6)
					{
						filt<-subset(nodeDynamics,bw>localMax$bw/4)
						xpos<- sample(as.numeric(as.vector(filt$year)),size=1)
						ypos<-subset(nodeDynamics,year==xpos)$bw
						text(xpos,ypos+localMax$bw/50,term,cex=0.7)
					}
				}
			}
			dev.off()
		}
		if(is.null(filename))"NaP"
		else  {
			paste(paste(getwd(), "/",sep=""), filename,sep="")
		}
		},error=function(ex){print("error")})
	}	

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	dynamic.page.rank.statistics <- function(graphs,years=NULL,th=1.001,cutmode="rate")
	{
		tryCatch({
			timeplot.node.statistics(graphs,years,attr="pagerank",th,cutmode)
		},error=function(ex){print("error")})
	}
	
# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	dynamic.kleinberg.statistics <- function(graphs,years=NULL,th=3,cutmode="rate")
	{
		tryCatch({
			timeplot.node.statistics(graphs,years,attr="kleinberg",th,cutmode)
		},error=function(ex){print("error")})
	}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	dynamic.evcent.statistics <- function(graphs,years=NULL,th=0.1,cutmode="fix")
	{
		tryCatch({
			timeplot.node.statistics(graphs,years,attr="eigen",th,cutmode)
		},error=function(ex){print("error")})
	}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	dynamic.closeness.statistics <- function(graphs,years=NULL,th=1.001,cutmode="rate")
	{
		tryCatch({
			timeplot.node.statistics(graphs,years,attr="closeness",th,cutmode)
		},error=function(ex){print("error")})
	}

# !!!!!!!!!!!!!!!!! TEST OK !!!!!!!!!!!!!!!!! #
	
	dynamic.betweenness.statistics <- function(graphs,years=NULL,th=4,cutmode="rate")
	{	
		tryCatch(
			timeplot.node.statistics(graphs,years,attr="betweenness",th,cutmode)
		},error=function(ex){print("error")})
	}
	
	
