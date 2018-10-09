####phylogenetic signal############################################
require(picante)
require(ape)
require(geiger)
library(dendextend)
require(ade4)
require(phytools)

setwd("C:/Users/silvia/Documents/Fission-Fusion")

load("FFsp.robj")

d1<-dsp
d1[,31:39]<-sapply(d1[,31:39], as.character)

d1$Species<-as.character(d1$Species)

#build phylo tree - read tree Bininda-Emons et al. 2007
mammal.tree<-read.nexus("MammalTree_nexus2.txt")

#verify d1 species in the tree
mammal.tree[[1]]$tip.label[match(d1[,1], mammal.tree[[1]]$tip.label)]

#remove tursiops sp.
d1<-d1[-72,]

#replace absent species by close species when possible
a<-c(1,5,18,30,49,57,65)
d1$Species[a]<-c("Antechinus_leo","Antechinus_bellus", "Capra_cylindricornis",
           "Chiropotes_israelita", "Orcaella_brevirostris", "Physeter_catodon", 
           "Sotalia_fluviatilis")

rownames(d1)<-d1[,1]

tree<-prune.sample(t(d1), mammal.tree[[1]]) #extract only sp of interest from tree
plot(tree, no.margin=T, label.offset=1, cex=0.9)
tree.dendro<-as.hclust(reorder(multi2di(tree))) #convert to hclust
tree.dendro<-as.dendrogram(tree.dendro)         #convert to dendrogram


#solve polytomies
treedi<-multi2di(tree)

##sort dataset based on phylogeny
d1.sort<-d1[tree$tip.label,]
d1.sort<-d1.sort[,which(apply(d1.sort, 2, function(x) sum(!is.na(x)))>20)]

##separate categorical and continous vars

a<-sapply(d1.sort, is.character)
d1.fact<-d1.sort[,a]
d1.fact<-data.frame(unclass(d1.fact))
rownames(d1.fact)<-d1.fact$Species

d1.cont<-d1.sort[,-which(a)]

##test phylogenetic signal of categorical variables
###1. transform tree topology to one that has all internal branch lengths 
###   multiplied by 0 (i.e. lambda=0) creating one giant basal polytomy
lambda0<-rescale(treedi,"lambda", 0)#transforms the tree topology 
par(mfrow=c(1,2))#remember from day1 session2 that this sets the graphical parameters so that the plotting device has 1 row and 2 colums, so we can now plot two trees next to each other.
plot(treedi)
plot(lambda0)



#Now find the maximum likelihood estimate of lambda for each categorical variable

getNamedVector <- function(colname, df = d1.fact){
  a <- d1.fact[,colname]
  names(a) <- rownames(df)
  return(a)
}

lambda.fact<-list()
fitlambda0<-list()
diff.aic<-list()
for (i in 2:ncol(d1.fact))
{
  trait<-na.omit(getNamedVector(names(d1.fact[i])))
  
  lambda.fact[[i]]<-fitDiscrete(treedi, dat=trait, model="ER", transform="lambda", 
                                control=list(method="subplex",niter=200,
                                             hessian=F))
  fitlambda0[[i]]<-fitDiscrete(lambda0, trait,control=list(method="subplex",
                                                           niter=200))
  
  diff.aic[[i]]<-fitlambda0[[i]]$opt$aic-lambda.fact[[i]]$opt$aic
}

save(lambda.fact, fitlambda0, diff.aic, file="PhyloSignalCategorical2017.robj")

##Use maximum likelihood to estimate lambda for continuous variables
getNamedVectorCont <- function(colname, df = d1.cont){
  a <- d1.cont[,colname]
  names(a) <- rownames(df)
  return(a)
}

lambda.cont<-list()
fitlambda0.cont<-list()
diff.aic.cont<-list()
for (i in 1:ncol(d1.cont)){
  trait<-c(na.omit(getNamedVectorCont(names(d1.cont[i]))))
  lambda.cont[[i]]<-fitContinuous(treedi, trait, model="lambda",control=list(
    method="Nelder-Mead",niter=200, hessian=T))
  fitlambda0.cont[[i]]<-fitContinuous(lambda0, trait, model="lambda",
                                      control=list(method="Nelder-Mead", niter=200))
  
  diff.aic.cont[[i]]<-fitlambda0.cont[[i]]$opt$aic-lambda.cont[[i]]$opt$aic
}
save(lambda.cont, fitlambda0.cont, diff.aic.cont, file="PhyloSignalCont.robj")


# alternatively, use phylotools
tree.phylo<-as.phylo(treedi)
lambda.phylo.cont<-apply(d1.cont, 2, function(x) phylosig(tree.phylo, x,
                                                method="lambda", test=T))

#plot tanglegram
treedi.d<-as.dendrogram(treedi.hclust)
d1.phylo$tip.label<-sort(d1.phylo$tip.label)
treedi$tip.label<-sort(treedi$tip.label)

tanglegram(d1.phylo, treedi)

d1.tree<-d1[match(tree$tip.label, rownames(d1)),]
d1.tree[order(match(d1.tree[,1], tree$tip.label)),]
