build.graph3 <- function(nodeIndex="int02_id_CLEANED.csv",trendIndex="int02_py.csv",filter=NULL,sep=";",node.cut=c("intentionality","INTENTIONALITY"),time.cut=1,sim="jaccard")
{
tryCatch({
	graphs <- NULL
	t1 <- read.csv2(nodeIndex,header=TRUE,sep=";")
	t2 <- read.csv2(trendIndex,header=TRUE,sep=";")
	h1 <- colnames(t1)
	h2 <- colnames(t2)
	index <- intersect(h1,h2)
	if(length(index)==0)stop(paste(paste(paste("there are no matching columns in ",nodeIndex)," and "),trendIndex))
	if(length(index) > 1)stop(paste("multiple columns matching: ",toString(index)))
	trend <- h2[h2 != index]
	if(length(trend) ==0)stop(paste("there is no trend column in file",trendIndex))
	if(length(trend) >1)stop(paste("more than two columns in file",trendIndex))
	if(!is.null(filter))evek <- filter
	else evek <- sort(unique(eval(parse(text=paste("t2$",trend)))))
	node <- h1[h1 != index]
	k <- 0
	if(length(node) == 0)stop(paste("there is no node column in file",nodeIndex))
	if(length(node) > 1)stop(paste("more than two columns in file",nodeIndex))
	if(time.cut > 1)
	{
		
		felbont=cut(evek,time.cut)
		period=as.data.frame(cbind(evek,as.vector(felbont)))
		colnames(period)=c(trend,"period")
		t1m=merge(t1,t2,by.x=index,by.y=index)
		colnames(t1m)=c(index,node,trend)
		t1p=merge(t1m,period,by.x=trend,by.y=trend)
		idok=split(t1p,t1p$period)
		for (i in 1:length(idok)) 
		{
			terms <- eval(parse(text=paste("idok[[i]]$",node)))
			actIndex <- eval(parse(text=paste("idok[[i]]$",index)))
			actYear <- eval(parse(text=paste("idok[[i]]$",trend)))
			if(length(terms) > 0)
			{
				terms=as.character(terms)
				td=table(terms,actIndex)
				s1 <- as.matrix(simil(td,method=sim))
				write(s1, paste("graph",idok[[i]]$period[1]))
				g=graph.adjacency(s1, mode="undirected", weighted=T, diag=F, add.colnames=NULL)
				g$year <- actYear[1]
				#g$year <- as.vector(idok[[i]]$period[1])
				g$interval <- g$year:actYear[length(actYear)]
				g$label <- idok[[i]]$period[1]
				g <- set.vertex.attribute(g,"id",index=V(g),colnames(s1))
				for(c in node.cut)
				{
					sub=subset(V(g),!V(g)$id==c)
					g <- subgraph(g,sub)
				}
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- g
			}
		}
		return(graphs)
	}
	else
	{
		t1m <- merge(t1,t2,by.x=index,by.y=index)
		colnames(t1m) <- c(index,node,trend)
		for(i in evek)
		{
			t.filtered <- subset(t1m,t1m[colnames(t1m)==trend]==i,select=c(index,node))
			no <- eval(parse(text=paste("t.filtered$",node)))
			ind <- eval(parse(text=paste("t.filtered$",index)))
			if(length(no)>0 || length(t.filtered$V1) > 0)
			{
				write.table(t.filtered,file=paste(i,".csv",sep=""),sep=";",row.names=FALSE)
				t2 <- read.table(paste(i,".csv",sep=""),sep=";",col.names=c("record","term"),header=TRUE)
				unlink(paste(i,".csv",sep=""))
				m1 <- table(t2$term,t2$record)
				s1 <- as.matrix(simil(m1,method="jaccard"))
				pro.n1=graph.adjacency(s1, mode="undirected", weighted=TRUE, diag=FALSE)
				n0a=set.vertex.attribute(pro.n1,"id",index=V(pro.n1),colnames(s1))
				n0a <- set.graph.attribute(n0a,"year",i)
				n0a$interval <- i
				n0a$label <- i
				for(c in node.cut)
				{
					sub=subset(V(n0a),!V(n0a)$id==c)
					n0a <- subgraph(n0a,sub)
				}
				g <- n0a
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- g
			}
		}
	}
	graphs
	},error=function(ex){print("error")})
}

build.graph2 <- function(fileName1="int02_id_CLEANED.csv",fileName2="int02_py.csv",header1 = TRUE,header2= TRUE,index="record",trend="year",node="term",filter=NULL,sep=";",node.cut=c("intentionality","INTENTIONALITY"),time.cut=1,sim="jaccard")
{
	tryCatch({
	if(time.cut<1)stop("time.cut must be positive integer")
	graphs <- NULL
	k <- 0
	t1=read.csv2(fileName1,header=header1,sep=sep)
	t2=read.csv2(fileName2,header=header2,sep=sep)
	h1 <- NULL
	h2 <- NULL
	temp <- NULL
	if(is.null(index))
	{
		index[1] <- colnames(t1[1])
		index[2] <- colnames(t2[1])
		node <- colnames(t1[2])
	}
	if(length(index) >= 2)
	{
		temp1 <- index[1]
		temp2 <- index[2]
		if(is.numeric(index[1]))
		{
			temp1 <- colnames(t1)[index[1]]
		}
		if(is.numeric(index[2]))
		{
			temp2 <- colnames(t2)[index[2]]
		}
		index[1] <- temp1
		index[2] <- temp2
	}
	else
	{
		if(is.numeric(index))
		{
			temp <- NULL
			temp[1] <- colnames(t1)[index]
			temp[2] <- colnames(t2)[index]
			index <- temp
		}
		else
		{
			temp <- c(index,index)
			index <- temp
		}
	}
	if(is.null(node))
	{
		node <- colnames(t1[2])
	}
	if(is.numeric(node))
	{
		node <- colnames(t1[node[1]])
	}
	else
	{
		node <- node[1]
	}
	h1 <- colnames(t1)
	ind <- h1[h1 == node]
	if(length(ind) == 0) stop(paste(paste(paste("there is no column '",node,sep=""),"' in file", sep=""),fileName1))
	if(length(ind) > 1) stop(paste(paste(paste("multpile column:",node),"in file"),fileName1))
	ind <- h1[h1 == index[1]]
	if(length(ind) == 0)stop(paste(paste(paste("there is no column '",index[1],sep=""),"in file",sep=""),fileName1))
	if(length(ind) > 1) stop(paste(paste(paste("multiple column:",index[1]),"in file"),fileName1))
	h2 <- colnames(t2)
	ind <- h2[h2 == index[2]]
	if(length(ind) == 0)stop(paste(paste(paste("there is no column '",index[2],sep=""),"' in file",sep=""),fileName2))
	if(!is.null(filter))
		evek <- years
	else
	{
		if(is.null(trend))
		{
			trend <- colnames(t2[2])
		}
		if(is.numeric(trend))
		{
			trend <- colnames(t2[trend[1]])
		}
		else
		{
			trend <- trend[1]
		}
		tr <- h2[h2 == trend]
		if(length(tr) > 1)stop(paste(paste(paste("multiple column:",trend),"in file"),fileName2))
		if(length(tr) == 0)stop(paste(paste(paste("no column '",trend),"' exists in file"),fileName2))
		evek <- sort(unique(eval(parse(text=paste("t2$",trend)))))
		if(!is.numeric(evek))stop("trend parameter must be numeric")
	}
	if(time.cut > 1)
	{
		felbont=cut(evek,time.cut)
		period=as.data.frame(cbind(evek,as.vector(felbont)))
		colnames(period)=c(trend,"period")
		t1m=merge(t1,t2,by.x=index[1],by.y=index[2])
		colnames(t1m)=c(index[1],node,trend)
		t1p=merge(t1m,period,by.x=trend,by.y=trend)
		idok=split(t1p,t1p$period)
		for (i in 1:length(idok)) 
		{
			actNode <- eval(parse(text=paste("idok[[i]]$",node)))
			actIndex <- eval(parse(text=paste("idok[[i]]$",index)))
			actYear <- eval(parse(text=paste("idok[[i]]$",trend)))
			if(length(actNode) > 0)
			{
				actNode=as.character(actNode)
				td=table(actNode,actIndex)
				s1 <- as.matrix(simil(td,method=sim))
				g=graph.adjacency(s1, mode="undirected", weighted=T, diag=F, add.colnames=NULL)
				g$year <- actYear[1]
				#g$year <- as.vector(idok[[i]]$period[1])
				g$interval <- g$year:actYear[length(actYear)]
				g$label <- idok[[i]]$period[1]
				g <- set.vertex.attribute(g,"id",index=V(g),colnames(s1))
				for(c in node.cut)
				{
					sub=subset(V(g),!V(g)$id==c)
					g <- subgraph(g,sub)
				}
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- g
			}
		}	
	}
	else
	{
		t1m <- merge(t1,t2,by.x=index[1],by.y=index[2])
		colnames(t1m) <- c(index[1],node,trend)
		for(i in evek)
		{
			t.filtered <- subset(t1m,t1m[colnames(t1m)==trend]==i,select=c(index[1],node))
			no <- eval(parse(text=paste("t.filtered$",node)))
			ind <- eval(parse(text=paste("t.filtered$",index[1])))
			if(length(no)>0)
			{
				write.table(t.filtered,file=paste(i,".csv",sep=""),sep=";",col.names=c("record","term"),row.names=FALSE)
				t2 <- read.table(paste(i,".csv",sep=""),sep=";",header=TRUE)
				unlink(paste(i,".csv",sep=""))
				m1 <- table(t2$term,t2$record)
				s1 <- as.matrix(simil(m1,method="jaccard"))
				pro.n1=graph.adjacency(s1, mode="undirected", weighted=TRUE, diag=FALSE)
				n0a=set.vertex.attribute(pro.n1,"id",index=V(pro.n1),colnames(s1))
				n0a <- set.graph.attribute(n0a,"year",i)
				n0a$interval <- i
				n0a$label <- i
				for(c in node.cut)
				{
					sub=subset(V(n0a),!V(n0a)$id==c)
					n0a <- subgraph(n0a,sub)
				}
				g <- n0a
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- g
			}
		}
	}
	graphs
	},error=function(ex){print("error")})
}

build.metric.graphs <- function(fileName="WOS_kw_index.csv",years=1993:2008,sep=";",header=TRUE,index="record",trend="year",node="term",node.cut=c("intentionality","INTENTIONALITY"),time.cut=1,sim="jaccard")
{
tryCatch({
	graphs <- NULL
	k <- 0
	t=read.csv2(fileName,header=header,sep=sep)
	if(TRUE)
	{
		index <- "record"
		node <- "term"
		trend <- "year"
		colnames(t) <- c(index,node,trend)
		evek <- sort(unique(t$year))
	}
	else
	{
		evek <- t$year
	}
	if(!is.null(years))evek <- years

	if(time.cut > 1)
	{
		felbont=cut(evek,time.cut)
		period=as.data.frame(cbind(evek,as.vector(felbont)))
		colnames(period)=c(trend,"period")
		t1p=merge(t,period,by.x=trend,by.y=trend)
		idok=split(t1p,t1p$period)
		for (i in 1:length(idok)) 
		{
			
			actNode <- eval(parse(text=paste("idok[[i]]$",node)))
			actIndex <- eval(parse(text=paste("idok[[i]]$",index)))
			actYear <- eval(parse(text=paste("idok[[i]]$",trend)))
			
			if(length(idok[[i]]$term) > 0)
			{
				idok[[i]]$term=as.character(idok[[i]]$term)
				td=table(idok[[i]]$term,idok[[i]]$record)
				d=dist(td, method = "binary")
				d=as.matrix(d)
				w=1-d
				g=graph.adjacency(w, mode="undirected", weighted=T, diag=F, add.colnames=NULL)
				g$year <- idok[[i]]$year[1]
				#g$year <- as.vector(idok[[i]]$period[1])
				#g$last.year <- idok[[i]]$evek[length(idok[[i]]$evek)]
				g$label <- idok[[i]]$period[1]
				g$interval <- g$year:actYear[length(actYear)]
				g <- set.vertex.attribute(g,"id",index=V(g),colnames(w))
				for(c in node.cut)
				{
					sub=subset(V(g),!V(g)$id==c)
					g <- subgraph(g,sub)
				}
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- g
			}
		}	
	}
	else
	{
		counter <- 0
		for(i in evek)
		{
			t.filtered <- subset(t,year==i,select=c(index,node))
			if(length(t.filtered$term)>0 || length(t.filtered$V1) > 0)
			{
				write.table(t.filtered,file=paste(i,".csv",sep=""),sep=";",row.names=FALSE)
				t2 <- read.table(paste(i,".csv",sep=""),sep=";",header=TRUE)
				unlink(paste(i,".csv",sep=""))
				m1 <- table(t2$term,t2$record)
				s1 <- as.matrix(simil(m1,method="jaccard"))
				pro.n1 <- graph.adjacency(s1, mode="undirected", weighted=TRUE, diag=FALSE)
				n0a <- set.vertex.attribute(pro.n1,"id",index=V(pro.n1),colnames(s1))
				n0a <- set.graph.attribute(n0a,"year",i)
				n0a$label <- i
				n0a$interval <- i
				for(c in node.cut)
				{
					sub=subset(V(n0a),!V(n0a)$id==c)
					n0a <- subgraph(n0a,sub)
				}
				g$edges <- length(E(g))
				g$nodes <- length(V(g))
				comp <- clusters(g,mode="strong")
				g$components <- comp$no
				graphs[[k <- k + 1]] <- n0a
			}
		}
	}
	graphs
	},error=function(ex){print("error")})
}

read.dynamic.graphs <- function(fileName,years,sep=";",directed=FALSE)
{
tryCatch({
	graphs <- list()
	t1 <- read.csv(fileName,sep=sep,header=TRUE)
	counter <- 0
	for(i in years)
	{
		t.filtered <- subset(t1,year==i)
		graph.frame <- data.frame(t.filtered)
		if(length(graph.frame$year) == 0)stop("invalid year attribute")
		vertices <- c(as.vector(graph.frame$vertex1),as.vector(graph.frame$vertex2))	
		vertices_u <- unique(vertices)
		n0a <- graph.data.frame(graph.frame, directed=directed)
		n0a<-set.vertex.attribute(n0a,"id",index=V(n0a),vertices_u)
		n0a <- set.graph.attribute(n0a,"year",i)
		#sub=subset(V(n0a),!V(n0a)$id=="intentionality")
		#n0a=subgraph(n0a,sub)
		graphs[[counter <- counter + 1]] <- n0a
	}
	graphs
	},error=function(ex){print("error")})
}

read.file.to.table <- function(fileName,separator)
{
	t1 <- read.csv(fileName,sep=separator,header=TRUE)
	t1
}