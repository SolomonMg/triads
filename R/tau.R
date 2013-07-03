# TODO: Add comment
# 
# Author: user01
###############################################################################


#R code for generating random graphs:
#requires packages ergm, intergraph

library(ergm)
#install.packages("intergraph")
library(intergraph)
library(igraph0)

#set up weighting vectors for clustering and hierarchy
clust.mask <- rep(0,16)
clust.mask[c(1,3,16)] <- 1
hier.mask <- rep(1,16)
hier.mask[c(6:8,10:11)] <- 0

#compute triad count and triad proportion for a given weighting vector
mask.stat <- function(my.graph, my.mask){
	n.nodes <- vcount(my.graph)
	n.edges <- ecount(my.graph)
	#set probability of edge formation in random graph to proportion of possible edges present in original
	p.edge <- n.edges/(n.nodes*(n.nodes +1)/2)
	r.graph <- as.network.numeric(n.nodes, density = p.edge)
	r.igraph <- as.igraph(r.graph)
	tc.graph <- triad.census(r.igraph)
	clust <- sum(tc.graph*my.mask)
	clust.norm <- clust/sum(tc.graph)
	return(c(clust,clust.norm))
}

#build 100 random graphs and compute their clustering and hierarchy measurements to create an empirical null distribution
emp.distro <- function(this.graph){
	clust <- matrix(rep(0,200), nrow=2)
	hier <- matrix(rep(0,200),nrow=2)
	for(i in c(1:100)){
		clust[,i] <- mask.stat(this.graph, clust.mask)
		hier[,i] <- mask.stat(this.graph, hier.mask)
		}
	my.mat <- rbind(clust, hier)
	rownames(my.mat) <- c("clust.ct", "clust.norm", "hier.ct", "hier.ct.norm")
	return(my.mat)
}
#find empirical p-value
get.p <- function(val, distro)
{
	distro.n <- sort(distro)
	distro.n <- distro.n - median(distro.n)
	val.n <- val - median(distro.n)
	p.val <- sum(abs(distro.n) > abs(val.n))/100
	return(p.val)
}

#load data
library(NetData)
data(kracknets, package = "NetData")

# Create sub-graphs based on edge attributes

# Reduce to non-zero edges and build a graph object

krack_full_nonzero_edges <- subset(krack_full_data_frame, (advice_tie > 0 | friendship_tie > 0 | reports_to_tie > 0))
head(krack_full_nonzero_edges)

krack_full <- graph.data.frame(krack_full_nonzero_edges)
summary(krack_full)

# Set vertex attributes
for (i in V(krack_full)) {
	for (j in names(attributes)) {
		krack_full <- set.vertex.attribute(krack_full, j, index=i, attributes[i+1,j])
	}
}
krack_advice <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "advice_tie")==0])
krack_friendship <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "friendship_tie")==0])
krack_reports_to <- delete.edges(krack_full, E(krack_full)[get.edge.attribute(krack_full,name = "reports_to_tie")==0])

#fix randomization if desired so results are replicable
#set.seed(3123)

start <- Sys.time()
g <- barabasi.game(10)
#plot(g, layout=layout.fruchterman.reingold)
distg <- emp.distro(g)
proc10 <- Sys.time()-start

start <- Sys.time()
g <- barabasi.game(100)
distg <- emp.distro(g)
proc100 <- Sys.time()-start

start <- Sys.time()
g <- barabasi.game(1000)
distg <- emp.distro(g)
proc1000 <- Sys.time()-start

start <- Sys.time()
g <- barabasi.game(5000)
distg <- emp.distro(g)
proc10000 <- Sys.time()-start

cbind(c(10,100,1000,5000), c(proc10, proc100, proc1000, proc10000))
#compute empirical distributions for each network
hc_advice <- emp.distro(krack_advice)
hc_friend <- emp.distro(krack_friendship)
hc_report <- emp.distro(krack_reports_to)

get.p(198, hc_full[1,])
get.p(194, hc_advice[1,])
get.p(525, hc_friend[1,])
get.p(1003, hc_report[1,])
get.p(979, hc_full[3,])
get.p(1047, hc_advice[3,])
get.p(1135, hc_friend[3,])
get.p(1314, hc_report[3,])

#generate 95% empirical confidence intervals for triad counts

#clustering
c(sort(hc_advice[1,])[5], sort(hc_advice[1,])[95])
c(sort(hc_friend[1,])[5], sort(hc_friend[1,])[95])
c(sort(hc_report[1,])[5], sort(hc_report[1,])[95])

#hierarchy
c(sort(hc_advice[3,])[5], sort(hc_advice[3,])[95])
c(sort(hc_friend[3,])[5], sort(hc_friend[3,])[95])
c(sort(hc_report[3,])[5], sort(hc_report[3,])[95])
