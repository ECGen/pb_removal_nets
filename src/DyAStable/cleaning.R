simil.levenshtein <- function(string1, string2, case=TRUE, map=NULL) {
		
	if(!is.null(map)) {
		m <- matrix(map, ncol=2, byrow=TRUE)
		s <- c(ifelse(case, string1, tolower(string1)), ifelse(case, string2, tolower(string2)))
		for(i in 1:dim(m)[1]) s <- gsub(m[i,1], m[i,2], s)
		string1 <- s[1]
		string2 <- s[2]
	}
 
	if(ifelse(case, string1, tolower(string1)) == ifelse(case, string2, tolower(string2))) return(0)
 
	s1 <- strsplit(paste(" ", ifelse(case, string1, tolower(string1)), sep=""), NULL)[[1]]
	s2 <- strsplit(paste(" ", ifelse(case, string2, tolower(string2)), sep=""), NULL)[[1]]
	
	l1 <- length(s1)
	l2 <- length(s2)
	
	d <- matrix(nrow = l1, ncol = l2)
 
	for(i in 1:l1) d[i,1] <- i-1
	for(i in 1:l2) d[1,i] <- i-1
	for(i in 2:l1) for(j in 2:l2) d[i,j] <- min((d[i-1,j]+1) , (d[i,j-1]+1) , (d[i-1,j-1]+ifelse(s1[i] == s2[j], 0, 1)))
	
	d[l1,l2]
}

string.simil <- function(x,y)
{
	x <- toString(x)
	y <- toString(y)
	sizeX <- nchar(x)
	sizeY <- nchar(y)
	n <- min(sizeX,sizeY)
	k <- 0
	for(i in 1:n){
		if(strtrim(x,i) != strtrim(y,i))
			break
		k <- k + 1
	}
	nMax <- max(sizeX,sizeY)
	if(k<2){
		simil<- 0
	}else{
		simil <- k/nMax
	}
	simil
}

string.jaccard <- function(string1,string2)
{
	x <- toString(string1)
	y <- toString(string2)
	sizeX <- nchar(x)
	sizeY <- nchar(y)
	n <- min(sizeX,sizeY)
	intersect <- 0
	union <- 0
	for(i in 1:n){
		if(substr(x,i,i) == substr(y,i,i))
			intersect <- intersect + 1
	}
	nMax <- max(sizeX,sizeY)
#	if(intersect<3){
		#simil<- 0
#	}else{
		simil <- intersect/nMax
#	}
	simil
}

simil.authors <- function(string1,string2)
{
	simil <- NULL
	x <- toString(string1)
	y <- toString(string2)
	f <- FALSE
	words1 <- strsplit(x, " ")
	words2 <- strsplit(y, " ")

	if(words1[[1]][length(words1[[1]])] == words2[[1]][length(words2[[1]])])
	{
		l <- NULL
		k <- NULL
		for(i in words1[[1]])
		{
			l <- paste(l,i)
		}
		for(i in words2[[1]])
		{
			k <- paste(k,i)
		}
		if(simil.levenshtein(l,k) > 2)
			simil <- FALSE
		else 
			simil <- TRUE
	}else
	{
		simil <- FALSE
	}
	simil
}

simil.complex <- function(string1,string2)
{
	simil <- NULL
	x <- toString(string1)
	y <- toString(string2)
	f <- FALSE
	if(simil.levenshtein(x,y) <3){
		simil <- 1
	}else
	{
		words1 <- strsplit(x, " ")
		words2 <- strsplit(y, " ")
		sim <- 0
		for( i in words1[[1]])
		{
			for(j in words2[[1]])
			{
				if(nchar(i) >2 && nchar(j) > 2 || i == j)
				{
					if(simil.levenshtein(i,j)<=3)
					{
					sim <- sim + 1
					}
				}
				if(nchar(i) <= 2 || nchar(j) <=2)
					f <- TRUE
			}
		}
		t <- max(length(words1[[1]]),length(words2[[1]]))
		if(sim == 0)simil <-0
		if(sim<t && sim>0)
		{
			simil<- 2
		}else
		{
			simil <- 0
		}
		if(sim==t)simil <- 1
		simil
	}
}

check.co.authors <- function(string1,string2,bool1,bool2,description)
{

	simil <- 0
	
	if(bool1){
		a1 <- subset(description,author==string1,select=co.author)
		name1 <- a1$co.author
	}else
	{
		a1 <- subset(description,co.author==string1,select=author)
		name1 <- a1$author
	}
	if(bool2)
	{
		a2 <- subset(description,author==string2,select=co.author)
		name2 <- a2$co.author
	}else
	{
		a2 <- subset(description,co.author==string2,select=author)
		name2 <- a2$author
	}
	for(k in 1:length(name1)){
		for(l in 1:length(name2)){
			f <- FALSE
			if(toString(name1[k]) == "" || toString(name2[l]) == "")
			{
				simil <- check.titles(string1,string2,bool1,bool2,description)
			}else	
			{
				w1 <- strsplit(toString(name1[k])," ")
				w2 <- strsplit(toString(name2[l])," ")
				for(i in w1[[1]])
				{
					for( j in w2[[1]])
					{	
						if(simil.levenshtein(i,j)<=3 && nchar(i) >2 && nchar(j) > 2)
						{
							simil <- 1
							break
						}
						if(nchar(i) <=2 && nchar(j) <= 2 && j == i)
						{
							f <- TRUE
						}
					}
				}
				if(simil == 0 && f)
					simil <- check.titles(string1,string2,bool1,bool2,description)
				if(simil == 1) break
			}
		}
	}
	simil
}
check.titles <-function(string1,string2,bool1,bool2,description)
{
	simil <- NULL
	if(bool1){
		a1 <- subset(description,author==string1,select=Title)
	}else
	{
		a1 <- subset(description,co.author==string1,select=Title)
	}
	if(bool2){
		a2 <- subset(description,author==string2,select=Title)
	}else
	{
		a2 <- subset(description,co.author==string2,select=Title)
	}
	if(toString(a1$Title) == "" || toString(a2$Title) == "")
	{
		simil <- 0
	}else
	{
		if(simil.levenshtein(toString(a1$Title),toString(a2$Title)) <=4)
		{
			simil <- 1
		}else
		{
			simil <- 0
		}
	}
	simil
}
"~" <- function(...) UseMethod("~") 
"~.default" <- .Primitive("~") 
"~.character" <- function(...) paste(...,sep="") 

#png(width=600,height=600)
#par(mar=c(3,8,1,1),family="serif")
#setwd("D:/work")  #working dir
library(igraph)     
library(proxy)       
library(lattice) 

sort.data.frame <- function(x,key, ...)
{
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}

#levenshteint még megírjuk

levenshtein <- function(string1, string2, case=TRUE, map=NULL) {
		
	if(!is.null(map)) {
		m <- matrix(map, ncol=2, byrow=TRUE)
		s <- c(ifelse(case, string1, tolower(string1)), ifelse(case, string2, tolower(string2)))
		for(i in 1:dim(m)[1]) s <- gsub(m[i,1], m[i,2], s)
		string1 <- s[1]
		string2 <- s[2]
	}
 
	if(ifelse(case, string1, tolower(string1)) == ifelse(case, string2, tolower(string2))) return(0)
 
	s1 <- strsplit(paste(" ", ifelse(case, string1, tolower(string1)), sep=""), NULL)[[1]]
	s2 <- strsplit(paste(" ", ifelse(case, string2, tolower(string2)), sep=""), NULL)[[1]]
	
	l1 <- length(s1)
	l2 <- length(s2)
	
	d <- matrix(nrow = l1, ncol = l2)
 
	for(i in 1:l1) d[i,1] <- i-1
	for(i in 1:l2) d[1,i] <- i-1
	for(i in 2:l1) for(j in 2:l2) d[i,j] <- min((d[i-1,j]+1) , (d[i,j-1]+1) , (d[i-1,j-1]+ifelse(s1[i] == s2[j], 0, 1)))
	
	d[l1,l2]
}

string.simil <- function(x,y)
{
	x <- toString(x)
	y <- toString(y)
	sizeX <- nchar(x)
	sizeY <- nchar(y)
	n <- min(sizeX,sizeY)
	k <- 0
	for(i in 1:n){
		if(strtrim(x,i) != strtrim(y,i))
			break
		k <- k + 1
	}
	nMax <- max(sizeX,sizeY)
	if(k<2){
		simil<- 0
	}else{
		simil <- k/nMax
	}
	simil
}

string.jaccard <- function(string1,string2)
{
	x <- toString(string1)
	y <- toString(string2)
	sizeX <- nchar(x)
	sizeY <- nchar(y)
	n <- max(sizeX,sizeY)
	m <- min(sizeX,sizeY)
	i <- 1
	k <- 1
	union <- 0
	inter <- 0
	while(i <= sizeX && k <= sizeY)
	{
		if(i <= sizeX && k <= sizeY && substr(x,i,i) == substr(y,k,k))
		{
			inter <- inter + 1
			union <- union + 1
			i <- i+1
			k <- k + 1
		}
		else if(i <= sizeX && k <= sizeY && substr(x,i,i) < substr(y,k,k) || k>sizeX)
		{
			union <- union + 2
			i <- i + 1
			k <- k + 1
		}else if(i <= sizeX && k <= sizeY && substr(x,i,i) > substr(y,k,k) || i > sizeY)
		{
			union <- union + 2
			k <- k + 1
			i <- i + 1
		}
	}
	simil <- inter / union
	simil
}
png(width=600,height=600)
par(mar=c(3,8,1,1),family="serif")
#setwd("D:/work")  #working dir
library(igraph)     
library(proxy)       
library(lattice)     
words=read.csv("int02_id.csv",header=FALSE,sep="\t")

opt.trial  <- function(words,threshold = 0.5)
{
	terms <- words$V1
	l <- TRUE
	clusters <- NULL
	for(t in terms)
	{
		while(l)
		{
			w <- sample(terms,size=1,replace = FALSE)
			sim <- string.simil(words_sorted[i],words_sorted[j])
			clusters <- rbind(clusters,data.frame(word=prototype,no=k))
			if(sim > thershold)
			{
			}
		}
	}
}

sort.data.frame <- function(x,key, ...)
{
	if (missing(key)) {
		rn <- rownames(x)
		if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
		x[order(rn, ...), , drop=FALSE]
	} else {
		x[do.call("order", c(x[key], ...)), , drop=FALSE]
	}
}

clean.opt <- function(words,threshold = 0.6)
{
	words_data_frame <- data.frame(words)
	words_sorted <- sort.data.frame(words_data_frame,"V2")
	compared <- NULL
	k <-1
	terms<- tolower(words_sorted$V2)
	index <- words_sorted$V1
	i <- 1
	while( i <= length(terms))
	{
		compared <- rbind(compared, data.frame(ind=index[i],term=terms[i]))
		l <- TRUE
		k <- i+1
		while(l)
		{
			if(string.simil(terms[i],terms[k]) <threshold)
			{
				l <- FALSE
				i <- k
			}
			else
			{
				compared <- rbind(compared,data.frame(ind=index[k],term = terms[i]))
				k <- k + 1
			}
		}
	}
	compared_original <- sort.data.frame(compared,"ind")
	write.csv(compared_original,"int02_id_CLEANED.csv",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

clean <- function(words,threshold=0.5)
{
	compared <- list()
	for(i in 1:length(words$V2))
	{
		k <- 1
		l <- TRUE
		while(k <=length(compared) && l )
		{
			if(string.simil(compared[[k]],words$V2[i]) > threshold)
			{
				#l <- FALSE
				#z[[1]] <- toString(words$V2[i])
				#z[[2]] <- words$V1[i]
				#compared[[k]][[length(compared[[k]])+1]] <- z
			}
			k <- k + 1
		}
		if(l)
		{
			z <- list()
			z[[1]] <- toString(words$V2[i])
			z[[2]] <- words$V1[i]
			compared[[length(compared) +1]] <- list()
			compared[[length(compared)]][[1]] <- z
		}
	}
	compared
}
clean.words <-function(words,method="prefix",threshold=0.5)
{
	#words_sorted <- sort(words$V1)
	compared <- list()
	for(i in 1:length(words$V2))
	{
		k <- 1
		l <- TRUE
		while(k <=length(compared) && l )
		{
			#if(string.simil(compared[[k]],words$V2[i]) > threshold)
			{
				l <- FALSE
				z[[1]] <- words$V2[i]
				z[[2]] <- words$V1[i]
				compared[[k]][[length(compared[[k]])+1]] <- z
			}
			k <- k + 1
		}
		if(l)
		{
			z <- list()
			z[[1]] <- words$V2[i]
			z[[2]] <- words$V1[i]
			compared[[length(compared) +1]] <- list()
			compared[[length(compared)]][[1]] <- z
		}
	}

	for(i in 1:length(words_sorted))
	{
		for( j in 1:length(words_sorted))
		{
			if(method=="levenshtein"){
				similarities[i,j] <- levenshtein(words_sorted[i],words_sorted[j])
				if(similarities[i,j] >3)similarities[i,j] <- 0
			}
			#else if(method == "prefix")
				#similarities[i,j] <- string.simil(words_sorted[i],words_sorted[j])
			#else if(method == "jaccard")
			#	similarities[i,j] <- string.jaccard(words_sorted[i],words_sorted[j])
			#else
			#	stop("method not allowed")
		}
	}
	stop()
		diag(similarities) <- 0
		similarities[similarities<=threshold] <- 0
	clusters <- NULL
	k <- 1
	while(length(words_sorted!=0)){
		prototype <- sample(words_sorted,size=1,replace=FALSE)
		words_sorted <- words_sorted[words_sorted!=prototype]
		clusters <- rbind(clusters,data.frame(word=prototype,no=k))
		filt <-subset(similarities,select = prototype)
		sims<-subset(filt,filt>0)
		for(j in rownames(sims)){
			exists <- subset(clusters,word==j)
			if(length(exists$word)==0)
			{
				clusters <- rbind(clusters,data.frame(word=j,no=k))
				words_sorted <- words_sorted[words_sorted!=j]
			}
		}
		k <- k + 1
	}
	clusters
}

clean.prefix <- function(words)
{
	words_sorted <- sort(words$V1) #sort(words$terms)
	similarities <- matrix(nrow=length(words_sorted),ncol=length(words_sorted))
	colnames(similarities) <- words_sorted
	rownames(similarities) <- words_sorted
	for(i in 1:length(words_sorted))
	{
		for( j in 1:length(words_sorted)){
			similarities[i,j] <- string.simil(words_sorted[i],words_sorted[j])
		}
	}
	diag(similarities) <- 0
	similarities[similarities<=0.5] <- 0
#for(i in 1:length(words_sorted))
#{
#	index <- which(similarities[i,] == max(similarities[i,]))
#	maxVertex <- max(similarities[i,])
#	similarities[i,] <- 0
#	for(j in 1:length(index))
#		similarities[i,index[j]] <- maxVertex
#}
#for(i in 1:length(words_sorted))
#{
#	index <- which(similarities[i,] != 0)
#	for(j in 1:length(index)){
#		similarities[index[j],i] <- similarities[i,index[j]]
#	}
#}
#similGraph=graph.adjacency(similarities, mode="undirected", weighted=TRUE, diag=FALSE)
#similGraph=set.vertex.attribute(similGraph,"id",index=V(similGraph),words_sorted)
#plot(similGraph)
	clusters <- NULL
	k <- 1
	while(length(words_sorted!=0)){
		prototype <- sample(words_sorted,size=1,replace=FALSE)
		words_sorted <- words_sorted[words_sorted!=prototype]
		clusters <- rbind(clusters,data.frame(word=prototype,no=k))
		filt <-subset(similarities,select = prototype)
		sims<-subset(filt,filt>0)
		for(j in rownames(sims)){
			exists <- subset(clusters,word==j)
			if(length(exists$word)==0){
				clusters <- rbind(clusters,data.frame(word=j,no=k))
				words_sorted <- words_sorted[words_sorted!=j]
			}
		}
		k <- k + 1
	}
	clusters
}


clean.jaccard <- function(words)
{
	words_sorted <- sort(words$terms) #sort(words$terms)
	similarities <- matrix(nrow=length(words_sorted),ncol=length(words_sorted))
	colnames(similarities) <- words_sorted
	rownames(similarities) <- words_sorted
	for(i in 1:length(words_sorted))
	{
		for( j in 1:length(words_sorted)){
			similarities[i,j] <- string.jaccard(words_sorted[i],words_sorted[j])
		}
	}
	diag(similarities) <- 0
	similarities[similarities<=0.5] <- 0
	clusters <- NULL

	k <- 1
	while(length(words_sorted!=0)){
		prototype <- sample(words_sorted,size=1,replace=FALSE)
		words_sorted <- words_sorted[words_sorted!=prototype]
		clusters <- rbind(clusters,data.frame(word=prototype,no=k))
		filt <-subset(similarities,select = prototype)
		sims<-subset(filt,filt>0)
		for(j in rownames(sims)){
			exists <- subset(clusters,word==j)
			if(length(exists$word)==0){
				clusters <- rbind(clusters,data.frame(word=j,no=k))
				words_sorted <- words_sorted[words_sorted!=j]
			}
		}
		k <- k + 1
	}
	clusters
}

clean.levenshtein <- function(words)
{
	words_sorted <- as.vector(sort(words$terms)) #sort(words$terms)
	similarities <- matrix(nrow=length(words_sorted),ncol=length(words_sorted))
	colnames(similarities) <- words_sorted
	rownames(similarities) <- words_sorted
	for(i in 1:length(words_sorted))
	{
		for( j in 1:length(words_sorted)){
			similarities[i,j] <- levenshtein(words_sorted[i],words_sorted[j])
		}
	}
	diag(similarities) <- 0
	similarities[similarities>3] <- 0
	clusters <- NULL
	k <- 1
	while(length(words_sorted!=0)){
		prototype <- sample(words_sorted,size=1,replace=FALSE)
		words_sorted <- words_sorted[words_sorted!=prototype]
		clusters <- rbind(clusters,data.frame(word=prototype,no=k))
		filt <-subset(similarities,select = prototype)
		sims<-subset(filt,filt>0)
		for(j in rownames(sims))
		{
			exists <- subset(clusters,word==j)
			if(length(exists$word)==0)
			{
				clusters <- rbind(clusters,data.frame(word=j,no=k))
				words_sorted <- words_sorted[words_sorted!=j]
			}
		}
		k <- k + 1
	}
	clusters
}

clean.cit <- function(filename1="int02_citing.csv",filename2="int02_cited_InJournal.csv")
{
	desc_citing <- read.csv(filename1,header=FALSE)
	desc_cited <- read.csv(filename2,header=FALSE)
	authors_citing <- desc_citing$V2
	authors_cited <- desc_cited$V2
	authors <- c(toupper(desc_citing$V2),toupper(desc_cited$V2))
	index <- c(desc_citing$V1,desc_cited$V1)
	i <- 1
	compared <- NULL
	while(i<=length(authors))
	{
		l <- TRUE
		j <- 1
		while(l && j <= length(compared$author))
		{
			if(simil.levenshtein(authors[i],compared$author[j])<2)
			{
				compared <- rbind(compared,data.frame(ind=index[i],author=compared$author[j]))
				l <- FALSE
			}
			j <- j + 1
		}
		if(l)
		{
			compared <- rbind(compared,data.frame(ind=index[i],author=authors[i]))
		}
		i <- i + 1
	}
	outfile1 <- data.frame(ind=compared$ind[1:length(authors_citing)],author=compared$author[1:length(authors_citing)])
	outfile2 <- data.frame(ind=compared$ind[(length(authors_citing) +1):length(compared$author)],author=compared$author[(length(authors_citing)+1):length(compared$author)])
	write.csv(outfile1,"int02_citing_CLEANED.csv",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	write.csv(outfile2,"int02_cited_InJournal_CLEANED.csv",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

clean.authors <- function(fileName="int02_au.csv")
{
	desc <- read.csv(fileName,header=FALSE,sep="\t")
	desc_sorted <- sort.data.frame(desc,"V2")
	authors <- desc_sorted$V2
	authors <- toupper(authors)
	index <- desc_sorted$V1
	i <- 1
	k <- 1
	compared <- NULL
	while(i <=length(authors))
	{
		l <- TRUE
		j <- 1
		while(l && j <= length(compared$author))
		{
			if(simil.authors(authors[i],compared$author[j]))
			{
				compared <- rbind(compared,data.frame(ind=index[i],author=compared$author[j]))
				l <- FALSE
			}
			j <- j + 1
		}
		if(l)
		{
			compared <- rbind(compared,data.frame(ind=index[i],author=authors[i]))
		}
		i <- i + 1
	}
	comp_orig <- sort.data.frame(compared,"ind")
	write.csv(comp_orig,"int02_au_CLEANED.csv",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}
cluster.names <- function(fileName="int")
{  
	description=read.csv(fileName,header=TRUE,sep=";")
	authors <- description$author
	coauthors <- description$co.author
	names <- c(as.vector(authors),as.vector(coauthors))
	names <- names[names != ""]
	authorMatrix <- matrix(nrow=length(names),ncol=length(names))
	rownames(authorMatrix) <- names
	colnames(authorMatrix) <- names
	similMatrix <- matrix(nrow=length(names),ncol=length(names))
	rownames(similMatrix) <- names
	colnames(similMatrix) <- names
	for(i in 1:length(names))
	{
		for(j in 1:length(names))
		{
			sim <- simil.complex(names[i],names[j])
			if(sim == 2)
			{
				if(i>length(authors) && j > length(authors))
				{
					similMatrix[i,j] = check.co.authors(names[i],names[j],FALSE,FALSE,description)
				}
				if(i>length(authors) && j <= length(authors))
				{	
					similMatrix[i,j] = check.co.authors(names[i],names[j],FALSE,TRUE,description)
				}	
				if(i<=length(authors) && j <= length(authors))
				{
					similMatrix[i,j] = check.co.authors(names[i],names[j],TRUE,TRUE,description)
				}
				if(i<=length(authors) && j > length(authors))
				{
					similMatrix[i,j] = check.co.authors(names[i],names[j],TRUE,FALSE,description)
				}
			}else
			{
				similMatrix[i,j] <- sim
			}
		}
	}

	similGraph=graph.adjacency(similMatrix, mode="undirected", weighted=TRUE, diag=FALSE)
	components <- decompose.graph(similGraph)
	name_clusters <- NULL
	c <- 1
	for ( i in 1:length(components))
	{
		adj<-get.adjacency(components[[i]])
		authors <- rownames(adj)
		while(length(authors != 0))
		{
			name <- sample(authors,1,replace=FALSE)
			authors <- authors[authors != name]
			name_clusters <- rbind(name_clusters,data.frame(name=name,cluster=c))
			x <- subset(adj,select=toString(name))
			x <- subset(x,x==1)
			for(j in rownames(x))
			{
				exists<- subset(name_clusters,name==j)
				if(length(exists$name)==0)
				{
					name_clusters <- rbind(name_clusters,data.frame(name=j,cluster=c))
					authors <- authors[authors != j]
				}
			}
			c <- c + 1
		}
	}
	name_clusters
}

create.name.graph <- function(fileName="authors.csv",name_clusters)
{
	png(width=600,height=600)
description=read.csv(fileName,header=TRUE,sep=";")
mergedGraph <- data.frame()
s <- NULL
for(i in 1:length(description$author))
{
	r <- subset(name_clusters,name==toString(description$author[i]))
	r3 <- subset(name_clusters,name==toString(description$co.author[i]))
	if(length(mergedGraph)!=0){
		s <- subset(mergedGraph,c1==r$cluster & c2==r3$cluster)
	}
	if(length(s$author) == 0 && as.vector(description$co.author[i] != ""))
	{
			if(length(r3$cluster) == 0)
				mergedGraph <- rbind(mergedGraph,data.frame(author=description$author[i],co.author=description$co.author[i],c1=r$cluster,c2=0))
			else{
				if(length(mergedGraph) != 0){
					s <- subset(mergedGraph,c1 == r$cluster | c2 == r3$cluster)
					if(length(s$author)!= 0 && as.vector(description$co.author[i] != "")){
						if(s$c1 == r$cluster)
						{
							mergedGraph <- rbind(mergedGraph,data.frame(author=s$author,co.author=description$co.author[i],c1=s$c1,c2=r3$cluster))
						}
						if(s$c2 == r3$cluster)
						{
							mergedGraph <- rbind(mergedGraph,data.frame(author=description$author[i],co.author=s$co.author,c1=r$cluster,c2=s$c2))
						}
					}else
						mergedGraph <- rbind(mergedGraph,data.frame(author=description$author[i],co.author=description$co.author[i],c1=r$cluster,c2=r3$cluster))
				}else
					mergedGraph <- rbind(mergedGraph,data.frame(author=description$author[i],co.author=description$co.author[i],c1=r$cluster,c2=r3$cluster))
			}
	}
	
}
mg <- graph.data.frame(mergedGraph,directed = FALSE)
plot(mg,vertex.color="green",vertex.size=40,vertex.label=c(unique(as.vector(mergedGraph$author)),unique(as.vector(mergedGraph$co.author))))
dev.off()
mg
}