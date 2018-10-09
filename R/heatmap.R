########## HEATPLOT CLUSTERS FROM PVCLUST#################
require(gplots)
require(Heatplus)

#load medians of scaled variables
dpop.median

#load data without imputation as an alternative and prepare it for analysis
dpop<-read.csv("FFpopulationsImp.csv")
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(dpop[,1], brain.eco[,3]),]
dpop$Species<-paste0(brain.eco[,1],"_", brain.eco[,2])
dpop[c(56,103,366,367),1]<-c("Chiropotes_albinasus", "Chiropotes_sagulatus",
                             "Tursiops_sp.", "Tursiops_sp.")
rownames(dpop)<-paste0(dpop$Species, "_", dpop$Pop)
dpop<-dpop[,c(1,16,18,25,28:32,34,36,39, 21:24, 35, 45,46)]
dpop[dpop==""]<-NA

#calculate summary stats per cluster for all variables and assign cluster id to
#each pop

pv95max<-pvpick(d1.pvclust, max.only=T)
sel.clusters<-c(1,2)#c(16,7,15,2,17,3, 18, 19)
pv95max.sel<-pv95max$clusters[sel.clusters]
#upload pop data
load("dpop.robj")

#calculate summary stats per cluster for all variables and assign cluster id to
#each pop
dpop<-data.frame(dpop, CVSubgroup=dpop$sg_size_sd/dpop$SubGrpSize_mean)

cluster<-c(rep(NA, times=nrow(dpop) ))
dpop<-data.frame(dpop, cluster=cluster) #add cluster column 
summ<-list()
for (i in 1:length(pv95max.sel)){
  summ[[i]]<-summary(dpop[dpop[,1] %in% pv95max.sel[[i]],])
  dpop[dpop[,1] %in% pv95max.sel[[i]],length(dpop)]<-i
}

#scale continuous variables
scaled.dpop<-apply(dpop[,11:28], 2, function(x) scale(x, center=T, scale=T))
dpop<-cbind(Species=dpop[,1], scaled.dpop, dpop[,29:38])

##select contiuous variables
cont<-c(2:30)
dpop.cont<-dpop[,cont]
##2.EStimate mean and standard error for each variables per cluster

dpop[,31:40]<-lapply(dpop[,31:40], as.factor)
dpop.mean<-ddply(dpop, .(cluster, Species), function(x) apply(x[,cont],2,
                            function(y) mean(y,na.rm=TRUE)*sum(!is.na(y))))

NAs<-ddply(dpop,.(cluster), function(x) apply(x[,cont], 2, 
                                              function(y) sum(!is.na(y))))

dpop.mean<-dpop.mean[,-2]


means<-data.frame(dpop.mean)
  
for (i in 2:length(dpop.mean)){
  for (j in 1:length(sel.clusters)){
   means[dpop.mean$cluster %in% NAs[j,1],i]<- 
     dpop.mean[dpop.mean$cluster %in% NAs[j,1],i]/NAs[j,i]
  }
  
}                             

dpop.mean<-ddply(means,.(cluster), function(x) apply(x[,2:30], 2, function(y) sum(y,na.rm=T)))


dpop.SD<-ddply(dpop, .(cluster, Species), function(x) apply(x[,cont],2,
                                function(y) sd(y,na.rm=TRUE)*sum(!is.na(y))))

dpop.SD<-dpop.mean[,-2]

SD<-data.frame(dpop.SD)

for (i in 2:length(dpop.SD)){
  for (j in 1:length(sel.clusters)){
  SD[dpop.SD$cluster %in% NAs[j,1],i]<- 
      dpop.SD[dpop.SD$cluster %in% NAs[j,1],i]/NAs[j,i]
  }
  
}               

n<-ddply(dpop,.(cluster), nrow)

#create dendrogram for clusters
dfordendro<-as.matrix(cbind(Cluster=c(1:6), VarX=c(6,7,9,10,14,16)))
rownames(dfordendro)<-dfordendro[,1]
dendro1<-as.dendrogram(hclust(dist(dfordendro[,2])))
plot(dendro1)



labvars<-c("Male:Female", "Density", "Mean.AR", "Mean.G","Mean.SG",
           "Max.SG","Rel.SG","Mean.FF")

medians<-as.matrix(dpop.mean[-7,c(8,9,11,16,19,21,26,27)])
rownames(medians)<-paste0("Cluster", " ", c(1:6))
colnames(medians)<-labvars

##with heatplus
###prepare binomial and categorical data
####EStimate freq. distrib. each variables per cluster
dpop.fact<-dpop.scaled[,c(31:40)]
dpop.fact<-as.data.frame(sapply(dpop.fact, as.factor))

freq<-list()
for (i in 1:(ncol(dpop.fact)-1)) {
  freq[[i]]<-ddply(dpop.fact, .(cluster), function(x) as.data.frame(
                                          table(x[,i],useNA="no")))
  total<-c()
  total<-ddply(freq[[i]], .(cluster), function(x) sum(x$Freq, na.rm=T))
  total<-rep(total$V1, each=length(levels(dpop.fact[,i])))
  
  freq[[i]]<-cbind(freq[[i]], rel=freq[[i]]$Freq/total)
}

binom<-data.frame(KinAssoc=c(1,1,1,1,NA,1),
                  SexSeg=c(NA,NA,1,1,NA,1), 
        Composition_A=c(0,0,0,1,0,0),
        Composition_B=c(0,0,1,0,0,0), 
        Composition_C=c(0,1,1,0,1,1), 
        Composition_D=c(0,0,1,1,1,0),
        Composition_E=c(0,0,1,1,1,1), 
        Composition_F=c(1,1,1,0,1,0), 
        Composition_G=c(0,0,1,0,1,0)
        )


png(filename="Heatplot2.png", width = 2000, height = 2000, 
    res=300)
par(oma=c(1,1,1,10), mar=c(3,1,7,10), font.lab=1)

plot(annHeatmap2(medians, legend=1, dendrogram=list(Col=list(status="no"), 
                  Row=list(status="no")), 
                 labels=list(Row=list(nrow=6, col=colors), 
                  Col=list(nrow=1, side=3)),
                 ann=list(Row=list(data=binom[6:1,], inclRef=F,
                 control=list(boxw=0.4,boxh=0.2, hbuff=0.005, vbuff=0.005,
                  cex.label=5, nacol=gray(0.7)))), scale="none",col=
                  bluered))


dev.off()
