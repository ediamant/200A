###############################################
### Lab 1: estimating diversification rates ###
###############################################

## 1. Lineage through time plots
library(phytools)
setwd("/Users/EllieDiamant/Documents/Fall2017/200A/")
darter.tree<-read.tree("etheostoma_percina_chrono.tre")
plotTree(darter.tree,ftype="i",fsize=0.4,type="fan",lwd=1)
obj<-ltt(darter.tree,log.lineages=FALSE)
darter.tree
is.binary(darter.tree)
darter.tree<-multi2di(darter.tree)
darter.tree
is.binary(darter.tree)
#We have randomly resolved the internal nodes that were multifurcating….

obj<-ltt(darter.tree,plot=FALSE)
obj
plot(obj,log.lineages=FALSE,main="LTT plot for darters")
plot(obj,show.tree=TRUE,log.lineages=FALSE,main="LTT plot for darters")

plot(obj,log.lineages=FALSE,log="y",main="LTT plot for darters",
     ylim=c(2,Ntip(darter.tree)))
## we can overlay the pure-birth prediction:
h<-max(nodeHeights(darter.tree))
x<-seq(0,h,by=h/100)
b<-(log(Ntip(darter.tree))-log(2))/h
lines(x,2*exp(b*x),col="red",lty="dashed",lwd=2)

# We might like to compare the observed LTT, to simulated LTTs assuming a pure-birth process of the same duration & resulting in the same total number of species. We can do this using the phytools function pbtree:
trees<-pbtree(b=b,n=Ntip(darter.tree),t=h,nsim=100,method="direct",quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT of darters compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(darter.tree))),lty="dashed",lwd=2,col="red")
## now let's overlay our original tree
ltt(darter.tree,add=TRUE,lwd=2)
ltt95(trees,log=TRUE)
title(main="LTT of darters compared to simulated LTTs")
ltt(darter.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2)

#####################################
## 2: Fitting pure and birth-death models to trees

library(phytools)
fitbd <- birthdeath(darter.tree)
fitbd
#d/b (extinction fraction) and b-d (net diversification)
#b and d rates separately
bd(fitbd)
#
One problem with this function is that the birth and death rates are conditioned on full sampling of the tree. 

###################################
## 3: the gamma statistic

g<-sapply(trees,function(x) ltt(x,plot=FALSE)$gamma)
hist(g,main=expression(paste("Distribution of ",gamma," from simulation")))
mean(g)
var(g)
obj<-ltt(darter.tree,plot=FALSE)
print(obj)
coal.tree<-rcoal(n=100)
plotTree(coal.tree,ftype="off")
obj<-ltt(coal.tree,log.lineages=FALSE,log="y")

darter.gamma <- obj$gamma # ltt returns a gamma value as one of its elements
trees<-pbtree(n=100,nsim=200,scale=max(nodeHeights(coal.tree)))
ltt95(trees,log=TRUE)
title(main="Simulated coalescent trees compared to pure-birth LTTs")
ltt(coal.tree,add=TRUE,log.lineages=FALSE,col="red",lwd=2,lty="dashed")


###############
## 3 Incomplete Sampling
# simulation test
library(geiger)
age <- 25.91862
richness <- 216
darterbirth =  (log(richness) - log(2))/age
darterbirth
richness <- 216
missing <- 15
#this simulates gamma values when trees are undersampled.
#we will grow trees with n=34 and prune them down to 13 taxa
num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(darterbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
# create a histogram of the null distribution
hist(g1_null)

#arrow indicates where the observed gamma falls in the null you just generated
arrows(darter.gamma, 40, darter.gamma, 0, col="red", lwd=2) 

# Which of the null values are smaller (more negative) than the data?
smallerNull<-g1_null<=darter.gamma
# How many TRUEs are there?
count<-sum(smallerNull)

# finally, what is the p-value?
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval

#####################
### 4: Excercise
## 1
snake.tree <- read.tree("homalops.phy")
snake.tree
is.binary(snake.tree)
obj<-ltt(snake.tree,plot=FALSE)
snake.gamma <- obj$gamma
snake.gamma
# gamma = -3.241081

## 2
# Under the assumption that the clade is fully sampled, this negative gamma combined with the p-value of 0.0012 suggests that speciation occurred early in the lineage, akin to radiation.

## 3
#Because 13 taxa are missing from this phylogeny, the gamma would be expected to be more negatively biased
library(geiger)
age <- 22
richness <- 34
darterbirth =  (log(richness) - log(2))/age
darterbirth

richness <- 34
missing <- 13
num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(darterbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
par(mar = c(5, 4, 4, 2))
hist(g1_null)
smallerNull<-g1_null<=snake.gamma
count <- sum(smallerNull)
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
# the pvalue of the gamma is still significant (<0.05) when accounting for tree pruning, suggesting that the lineage went through diversification earlier on.

## 4
fitbd <- birthdeath(snake.tree)
bd(fitbd)
#birth rate: 0.06839495
#death rate: 0.00000 

## 5
setwd("/Users/EllieDiamant/Documents/Fall2017/200A/")
humm.tree <- read.tree("hummtree.tree")
humm.tree
#1
#This is a tree from McGuire et al. 2014 on hummingbird diversity, time-calibration, and biogeography. There are 294 tips and 293 internal nodes
#This tree was provided to me from McGuire. I put it in my repository for reference. McGuire, Jimmy A., et al. "Molecular phylogenetics and the diversification of hummingbirds." Current Biology 24.8 (2014): 910-916.
#2
is.binary(humm.tree)
obj<-ltt(humm.tree,plot=FALSE)
fitbd <- birthdeath(humm.tree)
bd(fitbd)
# birth rate = 0.1970653 and death rate = 0.0000000
#3
#Tree is from 2014 when hummingbirds were considered to be ~340 taxa (some taxa debated). I will use the 340 for this tree.They cite the age of the clade as 22my
humm.gamma <- obj$gamma
age <- 22
richness <- 340
hummbirth =  (log(richness) - log(2))/age
hummbirth

richness <- 340
missing <- 46
num_simulations<-200 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values
for(i in 1:num_simulations) {
  sim.bdtree(hummbirth, d=0, stop = "taxa", n=richness)->sim_tree 
  drop.random(sim_tree, missing)->prune # prune down to the # of taxa in the phylogeny
  gammaStat(prune)->g1_null[i]
}
par(mar = c(5, 4, 4, 2))
hist(g1_null)
smallerNull<-g1_null<=humm.gamma
count <- sum(smallerNull)
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval
#this p-value suggests that the gamma is not significant when taking into account tree pruning