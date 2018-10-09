###########PCA SOCIAL TRAITS WITH MISSING DATA##############
library(ks)
 
load("FFpops.robj")
dpop<-scaled.d1
dpop[dpop==""]<-NA
dpop$Species<-as.character(dpop$Species)
#dpop<-dpop[-c(33, 339,340),]

dpop.c<-dpop[,c(1, 3:13)]
dpop.d<-dpop[,c(1, 14:19)]

##PCA with mean to replace missing values: not a good idea

rep.mean.PCA<-PCA(dpop.c[,c(1,3:15)], quali.sup=1)

#PCA with multiple replacement of missing values (the PCA is based only on the observed values)
###estimate number of dimensions needed for imputation
nd<-estim_ncpPCA(dpop.c[,2:12], ncp.min=1, ncp.max=6, scale=FALSE)

###impute and pca with 3 dimension for imputation
imp<-imputePCA(dpop.c[,2:12], 1)

pca.imp<-PCA(data.frame(Species=dpop$Species, imp$co), quali.sup=1) #without weights

w<-apply(dpop.c[,2:12], 2, function(y) sum(!is.na(y))/335)
pca.imp.w<-PCA(data.frame(Species=dpop$Species, imp$co), quali.sup=1, col.w=w) #with weight

#plot PCAS
###without weights - INDIV.

plot(pca.imp, choix="ind", xlim=c(-11, 18),habillage=1, col.var="black", invisible="quali",
     title="")
plot(pca.imp, choix="var", new.plot=T, cex=0.8)

###with weights - INDIV.

plot(pca.imp.w, choix="ind", xlim=c(-5, 10),ylim=c(-7,7), habillage=1, col.var="black", invisible="quali",
     title="")
plot(pca.imp.w, choix="var", new.plot=T, cex=0.8)
             

#######Factor Analysis for mixed data ###########################################

#improving input? check correlation
# correl<-cor(dpop[,-c(1:4, 29:37)], use="pairwise.complete.obs", method="pearson")
# 
# high.corr<-which(correl>=0.8, arr.ind=T) #remove PropM,soc.budget, FF_sd, AI_min, AI_sd,
#                                          #soc.budget, GrpMin and Max, FF_min & max
# high.corr2<-which(correl<=-0.8, arr.ind=T)
# 
# #test high cors for statistical significance
# cor.p<-list()
# cont.d<-dpop[,-c(1:4, 29:37)]
# for (i in 1:nrow(high.corr)) {
#   cor.p[[i]]<-cor.test(cont.d[,high.corr[i,2]], cont.d[,high.corr[i,1]], 
#                        alternative="two.sided", method="pearson")
# }

# #dpop.clean<-dpop[,-c(2:4,6,12,14,15,17,18,22,26:28)] #remove correlated
# 
# #remove variables that might include too much noise from paper sampling
# dpop.clean<-dpop.clean[,-c(3,6:8)]
# #sin.na<-apply(dpop.clean, 2, function(y) sum(!is.na(y))) # remove vars with less than 36 obs
# #dpop.clean<-dpop.clean[,!sin.na<36]
# sin.na.pop<-apply(dpop.clean, 1, function(y) sum(!is.na(y)))
# dpop.clean<-dpop.clean[!sin.na.pop<2,] #remove obs with no vars after removing vars.
# data.FA<-dpop.clean

data.FA<-dclean
#convert discrete variables to factors
data.FA$PrefferedAssoc<-as.integer(data.FA$PrefferedAssoc)
data.FA$PrefferedAssoc<-as.factor(data.FA$PrefferedAssoc)

data.FA$KinAssoc<-as.integer(data.FA$KinAssoc)
data.FA$KinAssoc<-as.factor(data.FA$KinAssoc)

data.FA$SexualSeg<-as.integer(data.FA$SexualSeg)
data.FA$SexualSeg<-as.factor(data.FA$SexualSeg)

data.FA$Subgroup_composition[data.FA$Subgroup_composition=="A"]<-1
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="B"]<-2
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="C"]<-3
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="D"]<-4
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="E"]<-5
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="F"]<-6
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="G"]<-7
data.FA$Subgroup_composition[data.FA$Subgroup_composition=="H"]<-8
data.FA$Subgroup_composition<-as.integer(data.FA$Subgroup_composition)
data.FA$Subgroup_composition<-as.factor(data.FA$Subgroup_composition)

data.FA$FFSexDiff<-as.integer(data.FA$FFSexDiff)
data.FA$FFSexDiff<-as.factor(data.FA$FFSexDiff)

data.FA$FF_tempvar<-as.integer(data.FA$FF_tempvar)
data.FA$FF_tempvar<-as.factor(data.FA$FF_tempvar)

data.FA$Species<-as.factor(data.FA$Species)

data.FA<-data.FA[,-c(2:5)]

save(data.FA, file="data.FA")

#data.FA[,2:11]<-apply(data.FA[,2:11], 2, function(x) scale(x, center=T, scale=T))

nd<-estim_ncpFAMD(data.FA, ncp.min=0, ncp.max=6, method.cv="Kfold")

#try removing 15, 9,4,3...AI_mean removed leads to best result
imp.FA<-imputeFAMD(data.FA[,-1], ncp=3, threshold = 1e-06, maxiter=1000)

#FAMD using the imputed data with probability of being present in categorical variables
fa<-FAMD(data.FA[,-1], tab.com=imp.FA$tab.disj)
save(fa,imp.FA, file="FAMD.robj")

#FAMD using the imputed data with imputed categories in categorical variables
#fa<-FAMD(imp.FA$completeObs[,-1], ncp=5)
load("data.FA")
load("FAMD.robj")
summary(fa)
dimdesc(fa)

n<-apply(data.FA, 2, function(x) length(x[-which(is.na(x))]))

##Plots
plot(fa, choix="ind", xlim=c(-10,10), ylim=c(-5,10), lab.var = F,
     cex=0.8)

##plot with ggplot
# extract population coords in PC for plotting
PC1 <- fa$ind$coord[which(fa$ind$coord[,1]>-17),1]
PC2 <- fa$ind$coord[which(fa$ind$coord[,1]>-17),2]
PC3 <- fa$ind$coord[which(fa$ind$coord[,1]>-17),3]
labs <- rownames(data.FA)[which(fa$ind$coord[,1]<30)]
PCs <- data.frame(cbind(PC1,PC2, PC3))
rownames(PCs) <- labs

ggplot(PCs, aes(PC1,PC2, label=rownames(PCs))) + 
  geom_text() 


# Now extract quantitative variables
#
sig12<-c(1:4, 7,9:11)
sig13<-c(2:6, 8,9,11)
vPC1 <- fa$quanti.var$contri[,1] * c(-1,1,-1,1,0,0,-1,0, -1,-1, 1)
vPC2 <- fa$quanti.var$contri[,2] * c(1,1, -1,-1,0,0,1,0,1,1,-1)
VPC3 <- fa$quanti.var$contri[,3] * c(0,1,1,1,-1,-1,0,1,1,0,1)
vlabs <- rownames(fa$quanti.var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2, VPC3))
rownames(vPCs) <- vlabs
colnames(vPCs) <- colnames(PCs)
#
# and plot them
#
pv <- ggplot() + theme(aspect.ratio=1) + theme_bw(base_size = 20) 
# no data so there's nothing to plot
# put a faint circle there, as is customary
# angle <- seq(-pi, pi, length = 50) 
# df <- data.frame(x = sin(angle), y = cos(angle)) 
# pv <- pv + geom_path(aes(x, y), data = df, colour="grey70") 
# #
##plot pops and variables together
pv <- pv + geom_point(data=PCs, aes(PC1,PC2, label=rownames(PCs), colour="grey")) 
ppv<- pv + geom_text(data=vPCs, aes(x=vPC1,y=vPC2,label=rownames(vPCs)), 
                     size=4) + xlab("PC1") + ylab("PC2")
ppv <- pv + geom_segment(data=vPCs, aes(x = 0, y = 0, xend = vPC1*0.9, 
                yend = vPC2*0.9), arrow = arrow(length = unit(1/2, 'picas')),
                color = "black")
ppv # show plot 

##plot pops kernels and variables
pkv<- ppv + geom_density_2d(data=PCs, aes(x=PC1, y=PC2)) + 
      geom_point(data=PCs, aes(PC1,PC2, label=rownames(PCs), colour="grey")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
pkv


################ KERNEL DENSITY ESTIMATION DIAZ(2016) Nature##############################
H <- Hpi(x=PCs[,1:2])      # optimal bandwidth estimation
est<- kde(x=PCs[,c(1,2)], H=H[c(1,2),c(1,2)], compute.cont=TRUE)     # kernel density estimation

# set contour probabilities for drawing contour levels
cl<-contourLevels(est, prob=c(0.5, 0.25, 0.05), approx=TRUE)

##adjust coords of variables to plot
a<-sig12
fa.var1<-fa$var$contrib[a,1] * c(-1,1,-1,1,-1, -1,-1, 1) *.4
fa.var2<-fa$var$contrib[a,2] * c(1,1, -1,-1,1,1,1,-1) * .4
fit<-data.frame(cbind(fa.var1, fa.var2))

pdf("FAMD_Pops.pdf", width=8, height=6)
par(mar=c(4,4,2,2))
plot(est, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, ylab="", xlab="", 
    cex.axis=0.75, ylim=c(-10, 10), xlim=c(-10, 10),las=1) 
plot(est,abs.cont=cl[1], labels=c(0.5),labcex=0.7, add=TRUE, lwd=0.75, col="grey30")
plot(est,abs.cont=cl[2], labels=c(0.95),labcex=0.7, add=TRUE, lwd=0.5, col="grey60")
plot(est,abs.cont=cl[3], labels=c(0.99),labcex=0.7, add=TRUE, lwd=0.5, col="grey60")
points( PCs[,1:2], pch=16, cex=0.25, col="black") 
text(fit, cex=0.75, 
     labels=c("Females", "PropMales", "PropYoung", "M:F", "GS-Mean", "SGS-Max", "SG-CV",
              "Rel-SGS"), 
     pos=c(3,4,2,4,2,2,3,1), offset=0.1)
arrows(0,0, fit[,1], fit[,2], code=2, col=1, lty=1, lwd=1.3, length=0.1)
mtext("PC1", cex=0.8, side=1, line=2)
mtext("PC2", cex=0.8, side=2, line=2) #, las=2)
dev.off()

###plot FAMD pops by mammal order
library(colorspace)
pal <- choose_palette()

order<-read.csv("FFSpeciesOrder.csv")
order$Order<-as.character(order$Order)
data.FA$Species<-as.character(data.FA$Species)
order$Spcomplete<-as.character(paste0(order[,1], "_", order[,2]))
order.data.FA<-order[match(data.FA$Species, order$Spcomplete),4]
order.data.FA<-as.factor(order.data.FA)
data.FA$Order<-order.data.FA
#colors<-c('#e41a1c','#377eb8','#ff7f00','#984ea3','#00441b','#ffff33','#b3de69','#fb9a99','#999999')
colors<-c( rep("white",8),'#e41a1c')#, rep("white", 0))
col.order<-data.frame(cbind(levels(data.FA$Order), colors))
colors<-col.order[match(data.FA$Order, col.order[,1]), 2]

#pdf("FAMD_Pops_Individuals2.pdf", width=8, height=6) 
par(mar=c(4,4,2,2))
plot(est,abs.cont=cl[1], xlab="", ylab="",labels=c(0.5),labcex=0.7, 
     lwd=0.75, col="grey30", cex.axis=0.75, ylim=c(-10, 10), xlim=c(-10, 10))
plot(est,abs.cont=cl[2], labels=c(0.95),labcex=0.7, add=TRUE, lwd=0.5, col="grey60")
plot(est,abs.cont=cl[3], labels=c(0.99),labcex=0.7, add=TRUE, lwd=0.5, col="grey60")

points(PCs[,1:2], pch=21, cex=0.8, col="black", bg=as.character(colors), xlab="", ylab= "", cex.axis=0.75,
     ylim=c(-10, 10), xlim=c(-11, 11))
legend("topright", legend=levels(data.FA$Order), col="black", pt.bg=as.character(col.order$colors), pch=21, box.lwd=0.6, cex=0.8, 
       y.intersp = 0.8)
mtext("PC1", cex=0.8, side=1, line=2)
mtext("PC2", cex=0.8, side=2, line=2)
dev.off()
