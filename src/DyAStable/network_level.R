
plot.clique.numbers <- function(cliques)
{
	tryCatch({
	barplot(height=cliques$cl,names.arg=cliques$year,las=2,cex.names=1,cex.axis=0.5,main="number of complete subgraphs")
	},error=function(ex){print("error")})
}

vertex.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	vertex.numbers <- NULL
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			no <- length(V(n0a))
			vertex.numbers <- rbind(vertex.numbers,data.frame(number = no,year = n0a$year))
		}
	}
	vertex.numbers
},error=function(ex){print("error")})
}

plot.vertex.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	vertex.numbers <- vertex.numbers.time(graphs,years)
	barplot(height=vertex.numbers$number,names.arg=vertex.numbers$year,las=2,cex.names=1,cex.axis=0.5,main="vertex numbers")
	dev.off()
},error=function(ex){print("error")})
}

vertex.gaining.speed <- function(graphs,years=NULL)
{
tryCatch({
	speed <- 0
	no <- vertex.numbers.time(graphs,years)
	last <- no[1]
	for(n in no)
	{
		speed <- speed + abs(last - n)
		last <- n
	}
	speed
},error=function(ex){print("error")})
}



plot.edge.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	edge.numbers <- edge.numbers.time(graphs,years)
	barplot(height=edge.numbers$number,names.arg=edge.numbers$year,las=2,cex.names=1,cex.axis=0.5,main="edge numbers")
	dev.off()
},error=function(ex){print("error")})
}

edge.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	edge.numbers <- NULL
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			no <- length(E(n0a))
			edge.numbers <- rbind(edge.numbers,data.frame(number = no,year = n0a$year))
		}
	}
	edge.numbers
},error=function(ex){print("error")})
}

plot.degree.dist.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	par(mfrow=c(2,2))
	dist<- NULL
	i <- 0
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			dist <- degree.distribution(n0a)
			barplot(height=dist,las=2,cex.names=1,cex.axis=0.5,main=n0a$label)	
		}
	}
	dev.off()
},error=function(ex){print("error")})
}

graph.diameter.time <- function(graphs,years=NULL)
{
tryCatch({
	diameters <- NULL
	i <- 0
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			d <- diameter(n0a)
			diameters <- rbind(diameters,data.frame(d = d,year = n0a$year))
		}
	}
	diameters
},error=function(ex){print("error")})
}

plot.graph.diameter.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	diameters <- graph.diameter.time(graphs,years)
	barplot(height=diameters$d,names.arg=diameters$year,las=2,cex.names=1,cex.axis=0.5,main="graph diameters")
	dev.off()
},error=function(ex){print("error")})
}
average.path.lengths.time <- function(graphs,years=NULL)
{
tryCatch({
	averages <- NULL
	i <- 0
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			avg <- average.path.length(n0a,directed=FALSE)
			averages <- rbind(averages,data.frame(avg = avg,year = n0a$year))
		}
	}
	averages
},error=function(ex){print("error")})
}


plot.average.path.lengths <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	averages <- average.path.lengths.time(graphs,years)
	barplot(height=averages$avg,names.arg=averages$year,las=2,cex.names=1,cex.axis=0.5,main="average path length")
	dev.off()
},error=function(ex){print("error")})
}

clique.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	cliques <- NULL
	i <- 0
	for(n0a in graphs)
	{

		if(is.null(years) || length(years[years==n0a$year]) != 0){
			cliqueno <- clique.number(n0a)
			cliques <- rbind(cliques,data.frame(cl=cliqueno,year = n0a$year))
		}
	}
	cliques
},error=function(ex){print("error")})
}

plot.clique.numbers.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	cliques <- clique.numbers.time(graphs,years)
	barplot(height=cliques$cl,names.arg=cliques$year,las=2,cex.names=1,cex.axis=0.5,main="number of complete subgraphs")
},error=function(ex){print("error")})
}
graph.cohesion.time <- function(graphs,years=NULL)
{
tryCatch({
	coh <- NULL
	i <- 0
	for(n0a in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			c <- graph.cohesion(n0a)
			coh <- rbind(coh,data.frame(cohesion = c,year=n0a$year))
		}
	}
	coh
},error=function(ex){print("error")})
}

plot.graph.cohesion.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	coh <- graph.cohesion.time(graphs,years)
	barplot(height=coh$c,names.arg=coh$year,las=2,cex.names=1,cex.axis=0.5,main="cohesion of the graph")
	dev.off()
},error=function(ex){print("error")})
}

graph.density.time <- function(graphs,years)
{
tryCatch({
	densities <- NULL
	i <- 0
	for(graph in graphs)
	{
		if(is.null(years) || length(years[years==n0a$year]) != 0)
		{
			d <- graph.density(graph)
			densities <- rbind(densities,data.frame(den=d,year=graph$year))
		}
	}
},error=function(ex){print("error")})
}

plot.graph.density.time <- function(densities)
{
tryCatch({
	png(width=600,height=600)
	barplot(densities$den,names.arg=densities$year,las=2,cex.names=1,cex.axis=0.5,main="density of graphs")
},error=function(ex){print("error")})
}

plot.graph.density.time <- function(graphs,years=NULL)
{
tryCatch({
	png(width=600,height=600)
	densities <- graph.density.time(graphs,years)
	barplot(densities$den,names.arg=densities$year,las=2,cex.names=1,cex.axis=0.5,main="density of graphs")
},error=function(ex){print("error")})
}

plot.vnum.diameter.ratio<- function(graphs,years)
{
tryCatch({
	png(width=600,height=600)
	diameters <- NULL
	vnum <- NULL
	i <- 0
	for(n0a in graphs)
	{
		d <- diameter(n0a)
		diameters <- rbind(diameters,data.frame(d = d,year = n0a$year))
		vnum <- c(vnum,length(V(n0a)))
	}
	barplot(diameters,vnum)
},error=function(ex){print("error")})
}