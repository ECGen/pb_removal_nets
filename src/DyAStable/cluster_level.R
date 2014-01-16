library(plotrix)
library(igraph)
	
sort.data.frame <- function(x,key, ...)
{
tryCatch({
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
},error=function(ex){print("error")})
}

szeletel <- function(filename1="int02_id_CLEANED.csv",filename2="int02_py.csv",time.cut=5,cut=c("intentionality","INTENTIONALITY"))
{
tryCatch({
	t1=read.csv2(filename1,header=T,sep=";")
	t2=read.csv2(filename2,header=F,sep=";")
	evek=sort(unique(t2$V2))
	felbont=cut(evek,time.cut)
	period=as.data.frame(cbind(evek,as.vector(felbont)))
	colnames(period)=c("evek","period")
	t1m=merge(t1,t2,by.x="record",by.y="V1")
	colnames(t1m)=c("record","term","evek")
	t1p=merge(t1m,period,by.x="evek",by.y="evek")
	idok=split(t1p,t1p$period)
	stat=NULL
	graphs <- NULL
	k <- 0
	for (i in 1:length(idok)) 
	{
		if(length(idok[[i]]$term) > 0)
		{
			idok[[i]]$term=as.character(idok[[i]]$term)
			td=table(idok[[i]]$term,idok[[i]]$record)
			d=dist(td, method = "binary")
			d=as.matrix(d)
			w=1-d
			g=graph.adjacency(w, mode="undirected", weighted=T, diag=F, add.colnames=NULL)
			g$year <- idok[[i]]$evek[1]
			#g$year <- as.vector(idok[[i]]$period[1])
			g$last.year <- idok[[i]]$evek[length(idok[[i]]$evek)]
			g$period <- idok[[i]]$period[1]
			g <- set.vertex.attribute(g,"id",index=V(g),colnames(w))
			for(c in cut)
			{
				sub=subset(V(g),!V(g)$id==c)
				g <- subgraph(g,sub)
			}
			graphs[[k <- k + 1]] <- g
		}
	}	
	graphs
	},error=function(ex){print("error")})
}

sanyi.fuggveny <- function(filename1 ="au_career.csv",filename2="py_career.csv")
{
tryCatch({
	t1=read.csv2(filename1,header=T)
	t2=read.csv2(filename2,header=F)
	evek=unique(t2$V2)
	felbont=cut(evek,6)
	period=as.data.frame(cbind(evek,as.vector(felbont)))
	colnames(period)=c("evek","period")
	t1m=merge(t1,t2,by.x="record",by.y="V1")
	colnames(t1m)=c("record","au","evek")
	t1p=merge(t1m,period,by.x="evek",by.y="evek")
	idok=split(t1p,t1p$period)
	stat=NULL
	png(width=600,height=600,filename=paste("sorozat",".png"))
	par(mfrow=c(2,3),mar=c(3,0,3,0))
	for (i in 1:length(idok)) 
	{
		idok[[i]]$au=as.character(idok[[i]]$au)
		td=table(idok[[i]]$au,idok[[i]]$record)
		d=dist(td, method = "binary")
		d=as.matrix(d)
		w=1-d
		g=graph.adjacency(w, mode="undirected", weighted=T, diag=F, add.colnames=NULL)
		V(g)$komp=clusters(g,mode="weak")$membership
		sort(table(V(g)$komp),decreasing=T)->kompmeret
		summary(kompmeret[kompmeret>2])[5]->hatar
		kompmeret[kompmeret>=hatar]->kompmeret
		as.numeric(names(kompmeret))->kompmeret
		v=NULL
		for (m in 1:length(V(g))) 
		{
			v0=any(kompmeret==V(g)$komp[m])
			v=append(v,v0)
		}
		V(g)$kompmeret=v
		kell=subset(V(g),V(g)$kompmeret==T)
		cog=subgraph(g,kell)
		kompok=decompose.graph(cog,mode="weak",min.vertices = 3)
		dens=NULL
		bw=NULL
		for (k in 1:length(kompok)) 
		{
			dens0=graph.density(kompok[[k]])
			bw0=mean(betweenness(kompok[[k]]))
			dens=append(dens,dens0)
			bw=append(bw,bw0)
		}

		a=clusters(cog,mode="weak")$no
		#b=max(clusters(cog,mode="weak")$csize)
		#c=mean(clusters(cog,mode="weak")$csize)
		d=mean(E(cog)$weight)
		#e=mean(V(cog)$fokszam)
		f=mean(dens)
		h=mean(bw)
		stat0=c(a,d,f,h)
		stat=rbind(stat,stat0)
	
		plot(cog
			,layout=layout.fruchterman.reingold(cog,niter=200,area=vcount(cog)^3)
			,axes=T
			,vertex.label=""
			,edge.color="green3"
			,main=rownames(summary(idok))[i]
			,vertex.size=2
			,vertex.color=V(cog)$komp
			,vertex.shape="square"
			,edge.width=1.5
			,edge.label=""
			,asp=0.7
		)
	}
	dev.off()

	rownames(stat)=rownames(summary(idok))
	colnames(stat)=c("#clus","mean.weight","mean.dens","mean.bw")


	png(width=700,height=250,file="au_statistics.png")
	par(mar=c(3,6,3,2),mfrow=c(1,4),cex=1)
	for (t in c(1,2,3,4)) 
	{
		bp=barplot(stat[,t],main=colnames(stat)[[t]],names.arg=rownames(stat),col="white",border=T,horiz=T,las=2);
		lines(stat[,t],bp[,1],main="",type="o",col="red",lwd=2);
		grid()
	}
	dev.off()
	},error=function(ex){print("error")})
}

create.test <-function(graphs,year=2008,nodeDeleteNumber=20,nodeAddNumber=20,edgeDeleteNumber=50,edgeAddNumber=50)
{
tryCatch({
	graph <- graphs[[16]]
	new.graph.list <- list()
	for(i in 1:6)
	{
		nodes.to.delete <- sample(length(V(graph)),size=nodeDeleteNumber)
		graphn <- delete.vertices(graph,nodes.to.delete)
		edges.to.delete <- sample(length(E(graphn)),size=edgeDeleteNumber)
		graphn <- delete.edges(graphn,edges.to.delete)
		edges.to.add <- sample(length(V(graphn)),size=edgeAddNumber*2)
		graphn <- add.edges(graph,edges.to.add)
		new.graph.list[[i]] <- graphn 
	}
	new.graph.list
	},error=function(ex){print("error")})
}

test.broken.clusters <- function()
{
tryCatch({
	graphs <- NULL
	graphs[[1]] <- graph.full(64)
	size <- 32
	for(i in 1:6)
	{
		g <- NULL
		for(j in 1:64/size)
		{
			g[[j]] <- graph.full(size)
			
		}
	}
	},error=function(ex){print("error")})
}

bommarito <- function(graphs,years=NULL)
{
tryCatch({
	eigen.comm <- NULL
	betweenness.comm <- NULL
	fastgreedy.comm <- NULL
	walktrap.comm <- NULL
	average.comm.betweenness <- 0
	average.comm.fastgreedy <- 0
	average.comm.eigen <- 0
	stab.eigen <- 0
	stab.fastgreedy <- 0
	stab.betweenness <- 0
	y <- 0
	c <- 0
	first <- TRUE
	for(g in graphs)
	{
		y <- y + 1
		if(is.null(years) || length(years[years==g$year]) != 0)
		{
		c <- c + 1
		#print(paste("number: ", y))
		#print(paste("year: ",g$year))
		#print(paste("number of nodes: ",length(V(g))))
		eigen.matrix <- leading.eigenvector.community(g)
		betweenness.matrix <- edge.betweenness.community(g)
		fast.greedy.matrix <- fastgreedy.community(g)
		#walktrap.matrix  <- walktrap.community(g)
		ei <- NULL
		if(length(eigen.matrix$merges) !=0)
		{
			ei <- community.le.to.membership(eigen.matrix$merges, 0, eigen.matrix$membership)
			
		}else
		{	
			ei$membership <- rep(0,length(V(g)))
			ei$csize <- length(V(g))
		}
		bw <- community.to.membership(g,betweenness.matrix$merges,length(betweenness.matrix$merges[,1]),betweenness.matrix$membership)
		fg <- community.to.membership(g,fast.greedy.matrix$merges,length(fast.greedy.matrix$merges[,1]),fast.greedy.matrix$membership)
		
		
		eigen.comm[[y]] <- ei
		betweenness.comm[[y]] <- bw
		fastgreedy.comm[[y]] <- fg
		#print("community_detection_done")
		same.nodes.eigen <- data.frame()
		same.nodes.betweenness <- data.frame()
		same.nodes.fastgreedy <- data.frame()
		
		if(!first)
		{
			for(i in 0:(length(V(graphs[[y-1]]))-1))
			{
				for(j in i:(length(V(graphs[[y-1]]))-1))
				{
					if(i != j && eigen.comm[[y-1]]$membership[i+1] == eigen.comm[[y-1]]$membership[j+1])
					{
						same.nodes.eigen <- rbind(same.nodes.eigen,data.frame(node1=V(graphs[[y-1]])[i]$id,node2=V(graphs[[y-1]])[j]$id))
					}
					if(i != j && betweenness.comm[[y-1]]$membership[i + 1] == betweenness.comm[[y-1]]$membership[j + 1])
					{
						same.nodes.betweenness <- rbind(same.nodes.betweenness,data.frame(node1=V(graphs[[y-1]])[i]$id,node2=V(graphs[[y-1]])[j]$id))
					}
					if(i != j && fastgreedy.comm[[y-1]]$membership[ i + 1 ] == fastgreedy.comm[[y-1]]$membership[j + 1])
					{
						same.nodes.fastgreedy <- rbind(same.nodes.fastgreedy,data.frame(node1=V(graphs[[y-1]])[i]$id,node2=V(graphs[[y-1]])[j]$id))
					}
				}
			}
			
			#print("last year nodes done")
			same.cluster.eigen <- 0
			same.cluster.betweenness <- 0
			same.cluster.fastgreedy <- 0
			#print(same.nodes.fastgreedy)
			#print(same.nodes.eigen)
			#print(same.nodes.betweenness)
			all.cluster.eigen <- length(same.nodes.eigen$node1)
			all.cluster.betweenness <- length(same.nodes.betweenness$node1)
			all.cluster.fastgreedy <- length(same.nodes.fastgreedy$node1)
			for(i in 1:length(same.nodes.eigen$node1))
			{
				node1.this.time <- which(V(g)$id == toString(same.nodes.eigen$node1[i]))
				node2.this.time <- which(V(g)$id == toString(same.nodes.eigen$node2[i]))
				if(length(node1.this.time) != 0 && length(node2.this.time) != 0)
				{
					if(eigen.comm[[y]]$membership[node1.this.time] == eigen.comm[[y]]$membership[node2.this.time])
					{
						same.cluster.eigen <- same.cluster.eigen + 1
					}
				}
			}
			
			print("Eigen vector done")
			for(i in 1:length(same.nodes.betweenness$node1))
			{
				node1.this.time <- which(V(g)$id == toString(same.nodes.betweenness$node1[i]))
				node2.this.time <- which(V(g)$id == toString(same.nodes.betweenness$node2[i]))
				if(length(node1.this.time) != 0 && length(node2.this.time) != 0)
				{
					if(betweenness.comm[[y]]$membership[node1.this.time] == betweenness.comm[[y]]$membership[node2.this.time])
					{
						same.cluster.betweenness <- same.cluster.betweenness + 1
					}
				}
			}
			print("betweenness done")
			
			for(i in 1:length(same.nodes.fastgreedy$node1))
			{
				node1.this.time <- which(V(g)$id == toString(same.nodes.fastgreedy$node1[i]))
				node2.this.time <- which(V(g)$id == toString(same.nodes.fastgreedy$node2[i]))
				if(length(node1.this.time) != 0 && length(node2.this.time) != 0)
				{
					if(fastgreedy.comm[[y]]$membership[node1.this.time] == fastgreedy.comm[[y]]$membership[node2.this.time])
					{
						same.cluster.fastgreedy <- same.cluster.fastgreedy + 1
					}
				}
			}
			print("fast greedy done")
			print(g$year)
			print(same.cluster.eigen)
			print(all.cluster.eigen)
			print(stab.eigen)
			print(same.cluster.betweenness)
			print(all.cluster.betweenness)
			print(stab.betweenness)
			print(same.cluster.fastgreedy)
			print(all.cluster.fastgreedy)
			print(stab.fastgreedy)
			stab.eigen <- stab.eigen + same.cluster.eigen / all.cluster.eigen
			stab.betweenness <- stab.betweenness + same.cluster.betweenness / all.cluster.betweenness
			stab.fastgreedy <- stab.fastgreedy + same.cluster.fastgreedy / all.cluster.fastgreedy
		}
		first <- FALSE
		average.comm.betweenness <- average.comm.betweenness + length(bw$csize)
		average.comm.fastgreedy <- average.comm.fastgreedy + length(fg$csize)
		average.comm.eigen <- average.comm.eigen + length(ei$csize)
		}
	}
	stab.eigen <- stab.eigen / y
	stab.betweenness <- stab.betweenness/ y
	stab.fastgreedy <- stab.fastgreedy / y
	average.comm.betweenness <- average.comm.betweenness/y
	average.comm.fastgreedy <- average.comm.fastgreedy/y
	average.comm.eigen <- average.comm.eigen/y
	x <- NULL
	x[1] <- stab.eigen
	x[2] <- stab.betweenness
	x[3] <- stab.fastgreedy
	z <- NULL
	z[1] <- average.comm.eigen
	z[2] <- average.comm.betweenness
	z[3] <- average.comm.fastgreedy
	print(x)
	print(z)
	plot(x,z,col=c("green","blue","red"))
	print(average.comm.betweenness)
	print(average.comm.fastgreedy)
	print(average.comm.eigen)
	print(stab.betweenness)
	print(stab.fastgreedy)
	print(stab.eigen)
	},error=function(ex){print("error")})
}

statistic.based.dynamic.skeleton <- function(graphs,years=NULL,method="eigenvector",show.clusters=TRUE,cluster.limit = 5,th = 0.04,colors=TRUE)
{
tryCatch({
	shift <- 0.1
	if(is.null(years))ys <- length(graphs)
	else ys <- length(years)
	z <- 0
	spanw <- 1/ys
	results <- NULL
	for( g in 1:length(graphs))
	{
		if(is.null(years) || length(years[years == graphs[[g]]$year]) != 0)
		{
			memb <- rep(1,length(V(graphs[[g]])))
			clusno <- 1
			csize <- length(V(graphs[[g]]))
			merge.matrix <- NULL
			if(method == "eigenvector")
			{
				print(graphs[[g]])
				merge<- leading.eigenvector.community(graphs[[g]])
				merge.matrix <- merge$merge
				if(length(merge.matrix[,1]) > (cluster.limit -1))
					steps <- length(merge.matrix[,1]) - (cluster.limit - 1)
				else
					steps <- 1
				if(length(merge.matrix) != 0)
				communities <- community.le.to.membership(merge.matrix,steps,merge$membership)
			}
			if(length(merge.matrix) != 0)
			{
				memb <- communities$membership + 1
				clusno <- length(communities$csize)
				csize <- communities$csize
			}
			spanh <- 1/clusno
			xpos <- NULL
			ypos <- NULL
			rad <- NULL
			radius <- 0.001
			if(is.null(results))
			{
				for(c in 1:clusno)
				{
					xpos <- c(xpos,0)
					radnow <- 0.02 + 0.0005* csize[c]/10
					rad <- c(rad,radnow)
					ypos <- c(ypos,0 + radius +(c-1)*spanh)
					#ypos  <- c(ypos,0 + spanh)
				}
				res <- NULL
				res$memb <- memb
				res$csize <- csize
				res$year <- graphs[[g]]$year
				res$pos <- ypos + 0.02
				#res$colors <- lead.colors
				results[[z <- z+1]] <- res
				filename <- paste(paste("cluster_skeleton",unclass(Sys.time()),sep=""),".png",sep="")
				png(width=800,height=600,filename=filename)
				par(mar=c(1,1,1,1))
				plot.new()
				if(colors)
					col <- "green"
				else 
					col <- "grey"
				if(show.clusters)
					symbols(xpos + 0.02 + shift,ypos + 0.02,rad,inches=FALSE,add=TRUE,bg=col)
				text(shift,1,labels=graphs[[g]]$label,cex=1.2)
			}else
			{
				actualposes <- NULL
				d <- NULL
				spanh <- 1/clusno
				radius <- 	0.001
				if(colors)
					col <- rep(c("red","blue","yellow","grey","black","purple","orange","brown","white"),3)
				else 
					col <- rep("grey",50)
				for(c in 1:clusno)
				{
					rad <- 0.02 + 0.0005 * csize[c]/10
					if(c==1) d <- rad
					xpos <- 0+z*spanw
					ypos <- 0 + radius + (c-1)*spanh
					actualposes <- c(actualposes,ypos)
					if(show.clusters)
						symbols(xpos + 0.02 +shift,ypos + 0.02,rad,bg=col[g],inches=FALSE,add=TRUE)
					for(l in 1:length(results[[z]]$csize))
					{
						ind.curr <- which(memb == c)
						memb.this.year <- V(graphs[[g]])[ind.curr]$id
						ind.last.year <- which(results[[z]]$memb == l)
						memb.last.year <- V(graphs[[g-1]])[ind.last.year]$id
						inter <- intersect(memb.this.year,memb.last.year)
						uni <- max(csize[c],results[[z]]$csize[l])
						stat <- length(inter)/uni
						if(stat > th)
						{
							poses <- results[[z]]$pos[l]
							arrows(xpos-spanw + 0.02+shift,poses,xpos +0.02+shift,ypos + 0.02,code=2,lwd=3)
						}
					}
				}
				text(xpos +shift,1,labels=graphs[[g]]$label,cex=1.2)
				res <- NULL
				res$memb <- memb
				res$csize <- csize
				res$year <- graphs[[g]]$year
				res$pos <- actualposes + 0.02
				results[[z <- z+1]] <- res
			}
		}
	}
	dev.off()
	paste(paste(getwd(), "/",sep=""), filename,sep="")
	},error=function(ex){print("error")})
}

leadnode.based.cluster.skeleton <- function(graphs,years=NULL,method="eigenvector",lead.node="closeness",show.clusters=TRUE,cluster.limit=5,colors=TRUE)
{
tryCatch({
	shift <- 0.1
	filename <- paste(paste("cluster_skeleton",unclass(Sys.time()),sep=""),".png",sep="")
	png(width=800,height=600,filename=filename)
	result <- NULL
	k <- 0
	z <- 0
	g <- 0
	if(is.null(years))ys <- length(graphs)
	else ys <- length(years)
	spanw <- 1/ys
	if(is.null(show.clusters))stop("show.clusters must be boolean type")
	for(g in 1:length(graphs))
	{
		if(is.null(years) || length(years[years==graphs[[g]]$year]) != 0)
		{
			memb <- rep(1,length(V(graphs[[g]])))
			clusno <- 1
			csize <- length(V(graphs[[g]]))
			print(csize)
			print(g)
			if(method == "eigenvector")
			{
				merge<- leading.eigenvector.community(graphs[[g]])
				merge.matrix <- merge$merge
				if(length(merge.matrix[,1]) >cluster.limit-1)
					steps <- length(merge.matrix[,1]) - (cluster.limit-1)
				else
					steps <- 1
				if(length(merge.matrix) != 0)
				{
					communities <- community.le.to.membership(merge.matrix,steps,merge$membership)
				}
			}

			if(length(merge.matrix) != 0)
			{
				memb <- communities$membership + 1
				clusno <- length(communities$csize)
				csize <- communities$csize
			}
			if(lead.node == "betweenness")
			{
				graphs[[g]] <- set.vertex.attribute(graphs[[g]],"bw",value=betweenness(graphs[[g]],directed=FALSE))
			}else if(lead.node == "closeness")
			{
				graphs[[g]] <- set.vertex.attribute(graphs[[g]],"cl",value=closeness(graphs[[g]]))
			}
			if(!is.null(result))
			{
				leaders <- NULL
				actualposes <- NULL
				spanh <- 1/clusno
				if(colors)
					col <- rep(c("red","blue","yellow","grey","black","purple","orange","brown","white"),3)
				else
					col <- rep("grey",50)
				lead.col <- NULL
				nodes <- result[[length(result)]]$lead
				d <- NULL
				sp <- 0.001
				for(i in 1:clusno)
				{
					ind <- which(memb == i)
					n <- intersect(V(graphs[[g]])[ind]$id,V(graphs[[g-1]])[nodes]$id)
					oldCluster <- NULL
					for(j in n)
					{
						pos <- which(V(graphs[[g-1]])$id == j)
						l <- which(nodes == pos-1)
						oldCluster <- c(oldCluster,l)
					}
					rad <- 0.02 + 0.0005 * csize[i]/10
					if(i==1) d <- rad
					xpos <- 0+z*spanw
					ypos <- 0 + sp +(i-1)*spanh
					actualposes <- c(actualposes,ypos)
					if(show.clusters)
						symbols(xpos + 0.02+shift,ypos + 0.02,rad,bg=col[g],inches=FALSE,add=TRUE)
					
					v <- 0
					if(length(n) > 0)
					{
						arrowcol <- result[[length(result)]]$colors[oldCluster]
						poses <- result[[length(result)]]$pos[oldCluster]
						arrows(xpos-spanw + 0.02+shift,poses,xpos +0.02+shift,ypos + 0.02,col=arrowcol,code=2,lwd=3)
					}
				
					if(lead.node=="betweenness")
						lead <- sample(which(V(graphs[[g]])[ind]$bw == max(V(graphs[[g]])[ind]$bw)),1)
					else if(lead.node == "closeness")
					{
						c <- which(V(graphs[[g]])[ind]$cl == max(V(graphs[[g]])[ind]$cl))
						p <- c[1]
						lead <- ind[p]
					}
					u <- FALSE
					b <- 0
					for(j in n)
					{ 	
						b <- b + 1
						if(j == V(graphs[[g]])[lead]$id)
						{
							lead.col <- c(lead.col,arrowcol[b])
							u <- TRUE
						}
					}
					if(!u)
					{
						lead.col <- c(lead.col,sample(colors(),size=1))
					}
					leaders <- c(leaders,lead)
				}
				
				text(xpos+0.01 +shift,1,labels=graphs[[g]]$label,cex=1.2)
				res <- NULL
				res$memb <- memb
				res$csize <- csize
				res$lead <- leaders
				res$year <- graphs[[g]]$year
				res$pos <- actualposes + 0.02
				res$colors <- lead.col
				result[[k <- k+1]] <- res
			}
			else
			{
				leaders <- NULL
				spanh <- 1/clusno
				xpos <- NULL
				ypos <- NULL
				rad <- NULL
				sp <- 0.001
				for(i in 1:clusno)
				{
					ind <- which(memb == i)
					if(lead.node == "betweenness")
						lead <- sample(which(V(graphs[[g]])[ind]$bw == max(V(graphs[[g]])[ind]$bw)),1)
					else if(lead.node == "closeness")
					{
						#p <- sample(which(V(graphs[[g]])[ind]$cl == max(V(graphs[[g]])[ind]$cl)),1)
						c <- which(V(graphs[[g]])[ind]$cl == max(V(graphs[[g]])[ind]$cl))
						p <- c[1]
						lead <- ind[p]
						xpos <- c(xpos,0)
						radnow <- 0.02 + 0.0005* csize[i]/10
						rad <- c(rad,radnow)
						ypos <- c(ypos,0 + sp +(i-1)*spanh)
					}
					leaders <- c(leaders,lead)
				}
				if(colors)
					col <- "green"
				else
					col <- "grey"
				lead.colors <- sample(colors(),size=length(leaders))
				res <- NULL
				res$memb <- memb
				res$csize <- csize
				res$lead <- leaders
				res$year <- graphs[[g]]$year
				res$pos <- ypos + 0.02
				res$colors <- lead.colors
				result[[k <- k+1]] <- res
				plot.new()
				#col <- sample(colors(),size=1)
				#print(rad)
				#print("zutyu")
				if(show.clusters)
					symbols(xpos + 0.02+shift,ypos + 0.02,rad,inches=FALSE,add=TRUE,bg=col)
				text(shift,1,labels=graphs[[g]]$label,cex=1.2)
			}
			z <- z + 1	
		}	
		g <- g + 1
	}
	dev.off()
	paste(paste(getwd(), "/",sep=""), filename,sep="")
	},error=function(ex){print("error")})
}