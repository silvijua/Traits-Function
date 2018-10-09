###############HIERARCHICAL CLUSTERING, PLOTS AND PCA ######################
require(cluster)
require(pvclust)
require(plyr)
require(gplots)
require(ade4)
require(plotrix)
require(phylobase)
require(dendextend)
setwd("/Users/silvijua/Documents/Fission-Fusion")
setwd("C:/Users/silvia/Favorites/Documents/Fission-Fusion")

# HIERARCHICAL CLUSTERING OF SPECIES--------

##If imputated data:
load(file="ImputedSpecies2017.robj") ##load sp. data
load(file="FFsp.robj")

d1<-d1scaled #for raw data
d1<-arrange(d1, Species)

d1[,2:15]<-apply(d1[,2:15], 2, function(x) scale(x, center=min(x, na.rm=T),
                  scale=max(x, na.rm=T)-min(x, na.rm=T)))

rownames(d1)<-d1[,1]
b<-c(2, 4, 22,46,54,61)
sp<-c("Antechinus_subtropicus", "Antechinus_agilis", "Capra_aegagrus",
      "Orcaella_heinsohni", "Physeter_macrocephalus",
      "Sotalia_guianensis")
d1[,1]<-as.character(d1[,1])
d1[,1][b]<-sp
d1<-as.data.frame(cbind(Sp.tree=rownames(d1), d1))
rownames(d1)<-d1[,2]

for (i in 17:25){
  d1[,i]<- factor(d1[,i])
}

#dsp<-d1 
#load("SpeciesData.robj")
a<-c(5,6, 8, 10, 12:15,17:22,25 ) #variables to include are the ones identified as important in the FAMD
#weight<-apply(d1[,c(8:20, 22:32)],2, function(x) sum(!is.na(x)))
#d1<-dsp
w.cont<-1/meanCV[c(7,8,10,15,18,20,25,26)]
w.disc<-rep(mean(w.cont), 7)
for (i in 17:25){
  d1[,i]<- factor(d1[,i])
}

popdist<-as.matrix(daisy(d1[,a], metric="gower"), weights=c(w.cont, w.disc)) #create distance matrix

# if raw data:

weight<-apply(d1[,a],2, function(x) sum(!is.na(x)))/10 #para datos con NAs
popdist<-as.matrix(daisy(d1[,a], metric="gower", weights=weight))

###To be able to run pvclust with distance gower, generate a different data set
###based on this matrix and using cmdscale
# Compute the eigenvalues
x <- cmdscale(popdist,1,eig=T)

# Plot the eigenvalues and choose the correct number of dimensions (eigenvalues close to 0)
plot(x$eig, 
     type="h", lwd=5, las=1, 
     xlab="Number of dimensions", 
     ylab="Eigenvalues")

# Recover the coordinates that give the same distance matrix with the correct number of dimensions    
x <- cmdscale(popdist,30)

# pvclust() clusters columns
set.seed(1407)
d1.pvclust <- parPvclust(data=t(x), method.hclust="average", method.dist="euclidean", 
                    nboot=1000, iseed=1407)
plot(d1.pvclust, cex=0.8)
pvrect(d1.pvclust, max.only=T)

save(d1.pvclust, file="Pvclust_VARShighvariation.robj")
load(file="Pvclust_imp_Eucl_average_weighted.robj")


hclaverunct.selvars<-d1.pvclust
hclawardeucl.selvars<-d1.pvclust
hclavereucl.selvars<-d1.pvclust
hclward2corr<-d1.pvclust
hclaverunct<-d1.pvclust

save(hclwardEucl,hclward2Eucl, hclAverEucl, hclward2cor,hclavercor,
     hclavercor.selvars, hclavereucl.selvars,hclwardeucl.selvars,
     hclcompeucl.selvars, hclcompcor.selvars, hclwardeucl.selvars, 
     hclward2eucl.selvars, hclwardcor.selvars, file="Hcl.robj")

hcl<-list(hclwardEucl,hclward2Eucl, hclAverEucl, hclward2cor,hclavercor,
          hclavercor.selvars, hclavereucl.selvars,hclwardeucl.selvars,
          hclcompeucl.selvars, hclcompcor.selvars, hclwardeucl.selvars,
          hclward2eucl.selvars,hclwardcor.selvars)

##plot and evaluate pvclust output
se<-list()
for(i in 1:length(hcl)){
  plot(hcl[[i]], cex=0.8)
  pvrect(hcl[[i]], alpha=0.95)
  seplot(hcl[[i]], ylim=c(0,0.2))
}

d1.pvclust<-hclaverunct  #select best hierarchical clustering model

# PLOT CIRCULAR CLUSTER OF SOCIAL SYSTEMS ---------------------------------------------------
load(file="Pvclust_VARShighvariation.robj")
pv95max<-pvpick(d1.pvclust, max.only=T) #identify top significant clusters


##define colors for clusters in plot
#color clusters
cluster.pop<-c(rep(NA, times=nrow(d1) ))
d.st<-cbind(d1, cluster=cluster.pop)
summ<-list()

colors<-brewer.pal(6, "Dark2") #rainbow(length(pv95max$clusters),alpha=0.9) #

for (i in 1:length(pv95max$clusters)){
  d.st[d.st[,2] %in% pv95max$clusters[[i]],length(d.st)]<-colors[i]
}

colors_to_use<-d.st$cluster
colors_to_use <- colors_to_use[order.hclust(d1.pvclust$hclust)]


d1.pvclust[[1]]$labels<-as.character(d.st[d1.pvclust[[1]]$labels, 1])
#change labels when names are not complete

b<-c(2, 4, 22,46,54,61)
sp<-c("Antechinus_subtropicus", "Antechinus_agilis", "Capra_aegagrus",
      "Orcaella_heinsohni", "Physeter_macrocephalus",
      "Sotalia_guianensis")
d1.pvclust[[1]]$labels<-as.character(d1.pvclust[[1]]$labels)
d1.pvclust[[1]]$labels[b]<-sp

superradial.phylog <- function (phylog, circle = 1, cleaves = 1, cnodes = 0, labels.leaves = names(phylog$leaves), 
                                clabel.leaves = 1, labels.nodes = names(phylog$nodes), clabel.nodes = 0, 
                                draw.box = FALSE, colors.leaves = 1) 
{
  if (!inherits(phylog, "phylog")) 
    stop("Non convenient data")
  
  leaves.number <- length(phylog$leaves)
  leaves.names <- names(phylog$leaves)
  nodes.number <- length(phylog$nodes)
  nodes.names <- names(phylog$nodes)
  
  if (length(labels.leaves) != leaves.number) 
    labels.leaves <- names(phylog$leaves)
  if (length(labels.nodes) != nodes.number) 
    labels.nodes <- names(phylog$nodes)
  if (circle < 0) 
    stop("'circle': non convenient value")
  leaves.car <- gsub("[_]", " ", labels.leaves)
  nodes.car <- gsub("[_]", " ", labels.nodes)
  opar <- par(mar = par("mar"), srt = par("srt"))
  on.exit(par(opar))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  dis <- phylog$droot
  dis <- dis/max(dis)
  rayon <- circle
  dis <- dis * rayon
  dist.leaves <- dis[leaves.names]
  dist.nodes <- dis[nodes.names]
  plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
               xaxt = "n", yaxt = "n", xlim = c(-2, 2), ylim = c(-2, 
                                                                 2), xaxs = "i", yaxs = "i", frame.plot = FALSE)
  d.rayon <- rayon/(nodes.number - 1)
  alpha <- 2 * pi * (1:leaves.number)/leaves.number
  names(alpha) <- leaves.names
  x <- dist.leaves * cos(alpha)
  y <- dist.leaves * sin(alpha)
  xcar <- (rayon + d.rayon) * cos(alpha)
  ycar <- (rayon + d.rayon) * sin(alpha)
  
  q2 <- alpha > pi/2 & alpha <= pi
  q3 <- alpha > pi & alpha <= 3*pi/2
  
  
  srt = alpha * 180/pi
  srt[q2] <- srt[q2] - 180
  srt[q3] <- srt[q3] - 180
  
  pos <- rep(4, length(alpha))
  pos[q2 | q3] <- 2
  
  
  if (clabel.leaves > 0) {
    for (i in 1:leaves.number) {
      segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
    }
    for (i in 1:leaves.number) {
      par(srt = srt[i])
      text(xcar[i], ycar[i], leaves.car[i], adj = 0, cex = par("cex") * 
             clabel.leaves, pos=pos[i], offset=0, col=colors.leaves[i])
      segments(xcar[i], ycar[i], x[i], y[i], col = grey(0.7))
    }
  }
  
  if (cleaves > 0) {
    for (i in 1:leaves.number) points(x[i], y[i], pch = 21, 
                                      bg = "black", cex = par("cex") * cleaves)
  }
  ang <- rep(0, length(dist.nodes))
  names(ang) <- names(dist.nodes)
  ang <- c(alpha, ang)
  for (i in 1:length(phylog$parts)) {
    w <- phylog$parts[[i]]
    but <- names(phylog$parts)[i]
    ang[but] <- mean(ang[w])
    b <- range(ang[w])
    a.seq <- c(seq(b[1], b[2], by = pi/180), b[2])
    lines(dis[but] * cos(a.seq), dis[but] * sin(a.seq))
    x1 <- dis[w] * cos(ang[w])
    y1 <- dis[w] * sin(ang[w])
    x2 <- dis[but] * cos(ang[w])
    y2 <- dis[but] * sin(ang[w])
    segments(x1, y1, x2, y2)
  }
  if (cnodes > 0) {
    for (i in 1:length(phylog$parts)) {
      w <- phylog$parts[[i]]
      but <- names(phylog$parts)[i]
      ang[but] <- mean(ang[w])
      points(dis[but] * cos(ang[but]), dis[but] * sin(ang[but]), 
             pch = 21, bg = "white", cex = par("cex") * cnodes)
    }
  }
  points(0, 0, pch = 21, cex = par("cex") * 2, bg = "red")
  if (clabel.nodes > 0) {
    delta <- strwidth(as.character(length(dist.nodes)), cex = par("cex") * 
                        clabel.nodes)
    for (j in 1:length(dist.nodes)) {
      i <- names(dist.nodes)[j]
      par(srt = (ang[i] * 360/2/pi + 90))
      x1 <- dis[i] * cos(ang[i])
      y1 <- dis[i] * sin(ang[i])
      symbols(x1, y1, delta, bg = "white", add = TRUE, 
              inches = FALSE)
      text(x1, y1, nodes.car[j], adj = 0.5, cex = par("cex") * 
             clabel.nodes)
    }
  }
  if (draw.box) 
    box()
  return(invisible())
}

hcl.phylog<-hclust2phylog(d1.pvclust$hclust) #transform hclust to phylog
png(filename="HCL_MaxLevel.png", width = 2000, height = 2000,
    units="px", res=300)#, compression="none")
superradial.phylog(hcl.phylog, circle=0.7,cleaves=0.3,labels.leaves=names(hcl.phylog$leaves),
                   clabel.leaves=0.8, colors.leaves=colors_to_use)

dev.off()


#plot smaller trees of the large clusters
d1.pvclust.phylo<-as.phylo(d1.pvclust$hclust)
d1.pvclus.phylo1<-prune(d1.pvclust.phylo, c(pv95max$clusters[[2]], 
                    "Syncerus caffer", "Brachyteles_hypoxanthus"))

pv95max<-pvpick(d1.pvclust, max.only=F) #identify top significant clusters

#plot cluster 1
colors_to_use<-c(rep("black",2), "orangered3", "black", "orangered3","black",
                 rep("orangered3",2), rep("black",2),rep("darkgoldenrod2", 3),
                 "white")

png(filename="ClusterSp_Cluster1.png", width = 1600, height = 1000, 
    res=250)
plot(d1.pvclus.phylo1, label.offset=0.004, tip.color=colors_to_use,
     root.edge=F)
dev.off()

#plot cluster 2
d1.pvclus.phylo2<-prune(d1.pvclust.phylo, c(pv95max$clusters[[1]], 
                "Syncerus caffer", "Brachyteles_hypoxanthus"))

colors_to_use<-c(rep("cadetblue",5), rep("navy", 5), "black", "black",
                 "lightseagreen", "navy", "black", "black", "navy",
                 "black","lightseagreen", "gold4", rep("black", 2),
                 "cadetblue", rep("black", 3), rep("darkgreen", 2), "black",
                 "gold4", rep("black",2), "lightseagreen",
                 rep("green3",2), rep("black",2),"gold4", 
                 rep("navy", 4), rep("gold4", 2), rep("black",5), 
                 "white", "gold4", rep("black",5) )

png(filename="ClusterSp_Cluster2.png", width = 1300, height = 2300, 
    res=250)
plot(d1.pvclus.phylo2, label.offset=0.004, tip.color=colors_to_use,
     root.edge=F, cex=0.8)
dev.off()

##COMPARE PHYLOGENY AND SOCIAL SYSTEMS -----------------------------------

###distance matrix for social data for imputated data
d1<-d1[match(tree$tip.label, d1[,1]),]
rownames(d1)<-d1[,1]

d1.phylo<-as.phylo(d1.pvclust$hclust)

pvclust.dist<-as.matrix(cophenetic(d1.phylo)) ###distance matrix of pvclust

#distance matrix for phylogeny
tree.dendro<-as.hclust(reorder(multi2di(tree)))
phylo.dist<-as.matrix(cophenetic(tree))

#run mantel test
m.test<-mantel.test(pvclust.dist, phylo.dist, nperm=1000, graph=T)
m.test.r<-mantel(pvclust.dist, phylo.dist, method="pearson", permutations=1000)


##################################################################


###If using Raw data:
###distance matrix for social data for raw data
d1<-d1[match(tree$tip.label, d1[,1]),]
rownames(d1)<-d1[,1]

weight<-apply(d1[,-1],2, function(x) sum(!is.na(x)))/100
weight[c(1,5:6,8:9,12,18,23:26,29,32)]<-weight[c(1,5:6,8:9,12,18,23:26,29,32)]*2
popdist<-as.matrix(daisy(d1[,-1], metric="gower", weights=weight))
popdist[which(is.na(popdist))]<-mean(popdist, na.rm=T)

d1.phylo<-as.phylo(d1.pvclust$hclust)
a<-c(2,5,18,29, 48,56,64)

#replace absent species by close species when possible
d1.phylo$tip.label[a]<-c("Antechinus_leo", "Antechinus_bellus", "Capra_walie", "Chiropotes_israelita",
           "Orcaella_brevirostris", "Physeter_catodon","Sotalia_fluviatilis")

pvclust.dist<-as.matrix(cophenetic(d1.phylo)) ###distance matrix of pvclust

#distance matrix for phylogeny
tree.dendro<-as.hclust(reorder(multi2di(tree)))
phylo.dist<-as.matrix(cophenetic(tree))

#run mantel test
m.test<-mantel.test(pvclust.dist, phylo.dist, nperm=1000, graph=T)
m.test.r<-mantel(pvclust.dist, phylo.dist, method="pearson", permutations=1000)


# ANALIZE CLUSTER CHARACTERISTICS -----------------------------------------

######Summarize variables per cluster and plot#####################

##select max clusters of the best output of pvclust
pv95max<-pvpick(d1.pvclust, max.only=T)

sel.clusters<-c(1:6)
pv95max.sel<-pv95max$clusters[sel.clusters]
#upload pop data
load("FFpops.robj")

#calculate summary stats per cluster for all variables and assign cluster id to
#each pop

cluster<-c(rep(NA, times=nrow(dpop) ))
dpop<-cbind(dpop, cluster=cluster) #add cluster column 
summ<-list()
for (i in 1:length(pv95max$clusters)){
  summ[[i]]<-summary(dpop[dpop[,1] %in% pv95max$cluster[[i]],])
  dpop[dpop[,1] %in% pv95max$cluster[[i]],length(dpop)]<-i
}

#scale continuous variables
scaled.dpop<-apply(dpop[,6:34], 2, function(x) scale(x, center=T, scale=T))
dpop<-cbind(Species=dpop[,1], scaled.dpop, dpop[,35:44])

save(dpop, file="PopScaled.robj")

##select contiuous variables
dpop.cont<-dpop[,2:30]

####plot each cluster by variable (cont. vars)######
##1. load standardized vars
load("PopScaled.robj")

##2.EStimate median and standard error for each variables per cluster

cols<-colnames(dpop.cont)

apply(dpop[,cols],2,function(x) wilcox.test(x[which(dpop$cluster==2)], 
                                            x[which(dpop$cluster==3)]))
dpop.median<-ddply(dpop, .(cluster), function(x) apply(x[,cols],2,median,na.rm=TRUE))
dpop.SD<-ddply(dpop, .(cluster), function(x) apply(x[,cols],2, sd,
                                                 na.rm=TRUE))

n<-ddply(dpop,.(cluster), nrow)
##3.plot 
      
vars<-c("Male:Female", "Density", "Mean.AR", "Mean.G","Mean.SG",
         "Max.SG","Rel.SG","Mean.FF")


###all vars together per cluster (continuous variables)
png(filename="BarplotAllVarsAllCluster.png", width = 2000, height = 2000, 
    res=250)

par(mfrow=c(length(pv95max$clusters),1), mgp=c(0,0.25,0), tck=0.1,  mar=c(0.8,1,0.2,1.2),
    cex.lab=1, cex.axis=0.8, oma=c(6.3,0.8,0.2,3), las=2)

sel.vars<-c(8,9,11,16,19,21,26,27)

for (i in 1:length(pv95max$clusters)) {
  bar<-barplot2(as.matrix(dpop.median[i,sel.vars]), width=1.4, space=c(0,0.1), 
                names.arg=NULL,axes=T,
                beside=T, col=colors[i],xpd=T, border=colors[i],
                cex.axis=0.8, ylim=c(-1,1),
                plot.ci=T, ci.l=as.matrix(dpop.median[i,sel.vars]-dpop.SD[i,sel.vars]), 
                ci.u=as.matrix(dpop.median[i,sel.vars]+dpop.SD[i,sel.vars]))
  
  axis(4, at=c(0, -0.5), label=c(paste0("Cluster", " ", i), paste0(
    "n=",n[i,2])), tick=F, cex.axis=1)
  
  if(i > (length(pv95max$clusters)-1)){
    par(mgp=c(0,0,0), mar=c(0.6,0.2,0,0),las=2 )
    axis(1, cex.axis=1.2, at = bar+0.2, tick=F, outer=T,
         label=vars, padj=0)
  }
  
}

dev.off

dpop.scaled<-dpop
save(dpop.scaled, dpop.median, file="PopScaled.robj")

####plot each cluster by variable (disc. vars)######

##1. load standardized vars (see code above)
####2.EStimate freq. distrib. each variables per cluster
a<-c(24:33)
dpop.fact<-as.data.frame(sapply(dpop[,a], as.factor))
levels(dpop.fact$cluster)<-c(1:length(levels(dpop[,29])))

freq<-list()
for (i in 1:9) {
  freq[[i]]<-ddply(dpop.fact, .(cluster), function(x) as.data.frame(table(x[,i], 
                                                                   useNA="no")))
  total<-c()
  total<-ddply(freq[[i]], .(cluster), function(x) sum(x[,3], na.rm=T))
  total<-rep(total[,2], each=length(levels(dpop.fact[,i])))
    
  freq[[i]]<-cbind(freq[[i]], rel=freq[[i]][,3]/total)
}

###3. plot
png(filename="BarplotCategorical53.png", width = 2000, height = 2000, 
    res=250)

par(mfcol=c(2,3), 
    mgp=c(0.2,0.3,0), tck=0.05,  mar=c(0.5,0.5,0.5,0.5),
    cex.lab=1, cex.axis=1, oma=c(1.5,1,1.5,10), las=2)

for (i in c(5,7:8)) {
  for (j in c(1:2,4,6:9)){
    if(i==5){
    bar.fact<-barplot2(as.matrix(freq[[i]][freq[[i]][,1]==j,4]), width=0.3, space=c(0,0.1), 
                     names.arg=NULL,ylim=c(0,1),
                     axes=T,beside=T, col=colors[j],xpd=F, border=colors[j],
                     cex.axis=0.8,
                     plot.ci=F)
    
    }
  else
    bar.fact<-barplot2(as.matrix(freq[[i]][freq[[i]][,1]==j,4]), width=0.3, space=c(0,0.1), 
                       names.arg=NULL,ylim=c(0,1),
                       axes=F,beside=T, col=colors[j],xpd=F, border=colors[j],
                       cex.axis=0.8,
                       plot.ci=F)
  if(i>7){
    axis(4, at=c(0.6, 0.45), label=c(paste0("Cluster", " ", j), paste0(
      "n=",n[j,2])), tick=F, cex.axis=1.2)
  }
  
  if(j > 8){
    par(mgp=c(0,0.3,0), mar=c(0,0.2,0.2,0),las=1 )
    axis(1, cex.axis=1.2, at = bar.fact[1:length(unique(freq[[i]][,2])),1], tick=F, outer=T,
         label=as.character(c(0:(length(unique(freq[[i]][,2]))-1))), padj=0)
    par (mar=c(0.5,0.5,0.5,0.5), tck=0, lty=0)
  }
  
  if(j==1){
  par(las=1)
  axis(3, at=bar.fact[5], label=colnames(d1.fact)[i], tick=F, cex.axis=1.2 )
  }
}
}

dev.off()


###run pca analysis of standardized median
dpop.median[is.na(dpop.median)]<-0

pca<-prcomp(~ PropFemale+PropMale+PropImmat+MaleFemaleRatio+PopDensity+PopSize+  
              AI_mean+AI_min+AI_max+social.budget+GrpSize_mean+GrpSize_min+
              GrpSize_max+SubGrpSize_mean+SubGrpSize_min+SubGrpSize_max+
              sg_size_sd+Sgsize_range+PropSG+Ffmean+FF_max, data=dpop.median[-12,],
            scale=F, na.action=na.exclude)

barplot(100*summary(pca)$importance[2,], ylab="% variance")

print(pca)  #load of each variable
summary(pca) #importance of components

png(filename="PCAclusterPop_median.png", width = 2000, height = 2000, 
       res=250)
par(mfrow=c(1,1),mgp=c(2, 0.6,0),mar=c(3,0.2,3,0.2), cex=0.8, cex.axis=0.9,
    cex.lab=0.9)
biplot(pca, choice=c(1,2), cex=c(1.5,1.2), col=c("black", "darkcyan"))
ndev.off()

plot3d(pca.sp$x[,1:3]) ##plot in 3d

#because clusters 5 and 9 are to close, we can pull them together

###run pca analysis of standardized continuous variables
d.st[is.na(d.st)]<-0
pca.sp<-prcomp(~ AI_mean+AI_max+GrpSize_mean+GrpSize_min+GrpSize_max+
                 SubGrpSize_mean+SubGrpSize_min+SubGrpSize_max+sg_size_sd
               +Sgsize_range+PropSG+Ffmean+FF_max+IndSubgrpDistance_mean, 
               data=d.st[!is.na(d.st[,25]), c(-1,-25)],
               scale=F, na.action=na.omit)
barplot(100*summary(pca.sp)$importance[2,], ylab="% variance")
biplot(pca.sp)
