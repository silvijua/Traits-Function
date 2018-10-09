#######################################################################################
#### Imputing missing data and prepare data for analysis -- SPECIES ########
#######################################################################################
require(phytools)
require(ape)
require(picante)
require(caper)
require(plyr)
require(mice)
require(lattice)
require(magrittr)

setwd("C:/Users/silvia/Documents/Fission-Fusion")

##load data at species level and change names
load("FFsp.robj")

#brain.eco<-read.csv("BodyEcoTraits.csv")
#brain.eco<-brain.eco[match(dsp[,1], brain.eco[,3]),]
#dsp[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])
rownames(d1)<-d1[,1]
dsp<-d1

# Dealing with phylogenetic signal ----------------------------------------

##estimate variance-covariance matrix for all the vars with phylogenetic signal using method of 
## Revell (2009)
mammal.tree<-read.nexus("MammalTree_nexus2.txt")

#verify species in the tree
a<-which(is.na(mammal.tree[[1]]$tip.label[match(dsp[,1], mammal.tree[[1]]$tip.label)]))

#replace absent species by close species when possible
replace.sp<-c("Antechinus_leo", "Antechinus_bellus", "Capra_walie",
              "Orcaella_brevirostris", "Physeter_catodon", 
              "Sotalia_fluviatilis")

#replace names by those used in the phylogeny
dsp[,1]<-as.character(dsp[,1])
dsp[a,1]<-replace.sp
rownames(dsp)<-dsp[,1]

tree<-prune.sample(t(dsp), mammal.tree[[1]]) #prune tree to include only species of interest
treedi<-multi2di(tree)

#check phylosignal and use phylopars to impute continuous variables with high
#phylosignal and mice for those without signal. Impute categorical with high
#phylosignal by using values from same genus.

###Prepare data for phylopars on the website software at species level####
phylo.vars<-c("PropImmat", "GrpSize_mean", "GrpSize_max")

dsp<-dsp[treedi$tip.label,]
phylopars.d<-dsp[,phylo.vars]

phylopars.d[is.na(phylopars.d)]<-""
phylopars.d<-sapply(phylopars.d, function(x) as.numeric(x))
rownames(phylopars.d)<-dsp[,1]

write.table(phylopars.d,  file="data.phylopars.txt", sep = "\t") ##export data
write.tree(treedi, file="NewickTree.txt") #export tree in newick format

#Run phylopars online

###phylopars in R####
library(Rphylopars)

dsp<-dsp[treedi$tip.label,]

d1.cont<-data.frame(dsp[,c(1,3,4,6,11,13,15,16,18)])
d1.cont[,1]<-as.factor(d1.cont[,1])

#fit phylopars: estimate covariance parameters
PPE<-phylopars(d1.cont, treedi, pheno_error=FALSE,
               phylo_correlated=T,species_identifier="Species", optim_limit=100)


#predict (impute) values               
d1.cont.pred<-phylopars.predict(PPE)
d1.cont.pred$predicted

#Cross-validation of imputed data
PPE_cv <- phylopars.crossvalidate(PPE) 
plot(PPE_cv)

#Replace NAs with Phylopars predicted values on species dataset
phylopars.sp<-read.delim("PhyloparsSp.txt", sep="\t")
abs<-apply(phylopars.sp[,-1],2,abs)
phylopars.sp <-data.frame(Species=as.character(phylopars.sp[,1]), abs)
phylopars.sp$Species <- as.character(phylopars.sp$Species)

phylopars.sp[phylopars.sp=="" ]<-NA

for (i in 2:nrow(phylopars.sp)){
  if (is.na(phylopars.sp$Species[i]) )
    phylopars.sp$Species[i] <- phylopars.sp$Species[i-1]
}

phylopars.sp<-ddply(phylopars.sp, .(Species), function(x) apply(x[,2:length(x)], 2,
                                                                mean))

phylopars.sp$Species<-gsub(" ", "_", phylopars.sp$Species)

##based on the ones that had better predictions after cross-validation
for ( i in c(5,24)){
  SpNA<-dsp[which(is.na(dsp[,i])),1]
  dsp[which(is.na(dsp[,i])),i] <- phylopars.sp[match(SpNA, phylopars.sp$Species), colnames(dsp)[i]]
}

write.csv(dsp, file="SpeciesImpute.csv")


###Impute data with multivariate imputation by chained equations (MICE) --------

PhyloVars<-read.csv("PhyloparsOutput2017.csv")
PhyloVars[,1]<-as.character(PhyloVars[,1])#add variables imputed with phylopars
rownames(PhyloVars)<-PhyloVars[,1]


dsp.fact<-read.csv("FFSpDiscImputed.csv")
rownames(dsp.fact)<-dsp.fact$Species
dsp.fact<-dsp.fact[as.character(dsp$Species),]
#verify species in the tree
a<-which(is.na(mammal.tree[[1]]$tip.label[match(dsp.fact[,1], mammal.tree[[1]]$tip.label)]))

#replace absent species by close species when possible
replace.sp<-c("Antechinus_leo", "Antechinus_bellus", "Capra_walie",
              "Orcaella_brevirostris", "Physeter_catodon", 
              "Sotalia_fluviatilis")

#replace names by those used in the phylogeny
dsp.fact[,1]<-as.character(dsp.fact[,1])
dsp.fact[a,1]<-replace.sp
rownames(dsp.fact)<-dsp.fact[,1]
dsp.fact<-dsp.fact[treedi$tip.label,]
dsp.fact[dsp.fact==""]  <- NA

dsp$HierarcStr<-dsp.fact$HierarcStr
dsp$PropImmat[which(is.na(dsp$PropImmat))]<-PhyloVars$PropImmat[which(
  is.na(dsp$PropImmat))]


#examine the proportion of usable cases given the distribution of NAs
p<-md.pairs(dsp)
p.summ<-round(p$mr/(p$mr + p$mm), 2)  #remove predictors with less than 30% usable cases

write.csv(p.summ, file="Psumm.csv")
write.csv(correl, file="correl.csv")


pred.variables<-list()

for (i in 2:ncol(p.summ)) {
  pred.variables[[i]]<-which(p.summ[colnames(p.summ)[i],] >= 0.3)
}

##add sample size/species as a variables that explains NA responses
cases<-read.csv("FFCases_countsamplesize.csv")
cases$Species<-as.character(cases$Species)
cases<-arrange(cases, Species)

#verify species in the tree
a<-which(is.na(mammal.tree[[1]]$tip.label[match(cases$Species, mammal.tree[[1]]$tip.label)]))

#replace absent species by close species when possible
replace.sp<-c("Antechinus_leo", "Antechinus_bellus", "Capra_walie",
              "Orcaella_brevirostris", "Physeter_catodon", 
              "Sotalia_fluviatilis")

#replace names by those used in the phylogeny
cases[a,1]<-as.character(c(rep(replace.sp[1],7), replace.sp[2], rep(replace.sp[3], 3),
                         "Chiropotes_sagulatus", replace.sp[4], 
                         rep(replace.sp[5],31), rep(replace.sp[6], 16)))
cases<-cases[!is.na(match(as.character(cases$Species),as.character(dsp$Species))),]
Sample<-data.frame(table(as.character(cases$Species)))

dsp<-arrange(dsp, Species)
Sample<-arrange(Sample, Var1)
dsp$Sample<-Sample[,2]

##check for correlations
correl<-cor(dsp[,2:15], use="pairwise.complete.obs", method="pearson")

#load prediction matrix
mat<-read.csv("PredictiveMat10.csv")

dsp.cut<-dsp[,-1]   ##remove variables that won't be imputed and don't predict
colnames(mat[,2:25])<-colnames(dsp.cut)
mat<-as.matrix(mat[,2:25])
rownames(mat)<-colnames(mat)
dsp.cut$PrefferedAssoc<-as.factor(dsp.cut$PrefferedAssoc)
dsp.cut$KinAssoc<-as.factor(dsp.cut$KinAssoc)
dsp.cut$SexualSeg<-as.factor(dsp.cut$SexualSeg)
dsp.cut$HierarcStr<-as.factor(dsp.cut$HierarcStr)
dsp.cut$Subgroup_composition<-as.factor(dsp.cut$Subgroup_composition)
dsp.cut$FFSexDiff<-as.factor(dsp.cut$FFSexDiff)
dsp.cut$FusionBehav<-as.factor(dsp.cut$FusionBehav)
dsp.cut$FissionBehav<-as.factor(dsp.cut$FissionBehav)
dsp.cut$FF_tempvar<-as.factor(dsp.cut$FF_tempvar)

#define boundaries for values for posterior distribution
post<-c("squeeze(dsp.cut[,1])", "squeeze(dsp.cut[,2])", "squeeze(dsp.cut[,3])",
        "squeeze(dsp.cut[,4])", "squeeze(dsp.cut[,5])", "squeeze(dsp.cut[,6])",
        "squeeze(dsp.cut[,7])", "squeeze(dsp.cut[,8])", "squeeze(dsp.cut[,9])",
        "squeeze(dsp.cut[,10])", "squeeze(dsp.cut[,11])", "squeeze(dsp.cut[,12])",
        "squeeze(dsp.cut[,13])", "squeeze(dsp.cut[,14])", rep('', 9))

#create 20 imputation values per variable
mice.sp<-mice(dsp.cut[,-24], m=20, method=c("norm", '' , rep("norm",3), 
                                      rep("pmm",2),
                                      rep("norm", 4), "pmm", 
                                      rep("norm",2), "polyreg", "logreg",
                                      rep("polyreg",7)), 
              predictorMatrix=mat[-24,-24], visitSequence="monotone", maxit=50, 
              print=TRUE)

save(mice.sp, file="ImputedSp2017Truncated.robj")

#view imputation values
mice.sp$imp

#plot to inspect results
stripplot(mice.sp)

#squeeze to truncate
for(i in c(1, 3:14)){
mice.sp$imp[[i]]<-apply(mice.sp$imp[[i]], 2, function(x) squeeze(x, 
                                                  range(dsp.cut[,i], na.rm=T)))
}

#check for convergence
plot(mice.sp)

#check imputated values
densityplot(mice.sp, data = ~ PopDensity + PopSize,
            layout = c(2,1))
#to review: AI_mean, SG_mean, prefferedAssoc,KinAssoc, SexualSeg, SubgroupCompo,
#FFSexDiff, FissionBehav, FF_tempvar


#fill out the original data set with the all imputation values for each variable
comp.pops<-complete(mice.sp, "long")

save(comp.pops, file="MiceComplete2017.robj")

comp.sp.cont<-ddply(comp.pops, .(.id), function(x) apply(x[,3:16], 
                                                         2, median))

comp.sp.disc<-ddply(comp.pops, .(.id), function(x) apply(x[,17:25], 2, function(y)
  names(sort(table(y), decreasing=T)[1])))

d.sp<-data.frame(Species=comp.sp.cont$.id, comp.sp.cont[,-1], comp.sp.disc[,-1])

d1<-d.sp


###Go back to preparing species and check for correlations

save(d1, file="ImputedSpecies2017.robj")


# Imputating population data ----------------------------------------------
####
dpop<-read.csv("FFpopulationsImp.csv")   #social data at pop level
rownames(dpop)<-paste0(dpop$Species, "_", dpop$Pop)

###Prepare data for phylopars at population level
dpop<-read.csv("FFpopulationsImp.csv")
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(dpop[,1], brain.eco[,3]),]
dpop[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])
dpop[c(56,103,366,367),1]<-c("Chiropotes_albinasus", "Chiropotes_sagulatus",
                              "Tursiops_sp.", "Tursiops_sp.")
rownames(dpop)<-paste0(dpop$Species, "_", dpop$Pop)

dpop.cont<-dpop[,c(1,13,14,25:28,32)]
dpop.cont$Species<-as.factor(dpop.cont$Species)
dpop.cont[is.na(dpop.cont)]<-""

##pool together all observations by species by var
dpop.phylopars<-ddply(dpop.cont, .(Species), function(x) apply(x[,-1], 2, 
                                                           function(y) paste(y, collapse=";")))
rownames(dpop.phylopars)<-dpop.phylopars[,1]
write.table(dpop.phylopars[,-1],  file="dpopphylopars.txt", sep = "\t")

#run phylopars and Replace NAs with Phylopars on population dataset

##If replacing values with phylo pars values in species data set:
phylopars.pop<-read.delim("PhyloparsPop.txt", sep="\t")
abs<-apply(phylopars.pop[,-1],2,abs)
phylopars.pop <-data.frame(Species=as.character(phylopars.pop[,1]), abs)
phylopars.pop$Species <- as.character(phylopars.pop$Species)

phylopars.pop[phylopars.pop=="" ]<-NA

for (i in 2:nrow(phylopars.pop)){
  if (is.na(phylopars.pop$Species[i]) )
    phylopars.pop$Species[i] <- phylopars.pop$Species[i-1]
}

phylopars.pop<-ddply(phylopars.pop, .(Species), function(x) apply(x[,2:length(x)], 2,
                                                                    mean))

phylopars.pop$Species<-gsub(" ", "_", phylopars.pop$Species)

d1.cont<-d1[,c(1:9,13:17,19)]

for ( i in c(4,8:9,11:12,14)){
  SpNA<-d1.cont[which(is.na(d1.cont[,i])),1]
  d1.cont[which(is.na(d1.cont[,i])),i] <- phylopars.sp[match(SpNA, phylopars.sp$Species), i]
}

for ( i in c(3,5,13)){
  SpNA<-d1.cont[which(is.na(d1.cont[,i])),1]
  d1.cont[which(is.na(d1.cont[,i])),i] <- phylopars.pop[match(SpNA, phylopars.pop$Species), i]
}

write.csv(d1.cont, file="PopsContImpute.csv")

#load imputed data for populations
pop.cont<-read.csv("PopsContImpute.csv")
pop.disc<-read.csv("PopDiscImpute.csv")
popdata<-data.frame(pop.cont, pop.disc[,-c(1:2)])
colnames(popdata)[1]<-"Pop"
popdata[popdata==""]<-NA

save(popdata, file="PopImputed.robj")


# #impute dpop using mice -------------------------------------------------

##change names for full version
dpop<-read.csv("FFpopulationsImp.csv")
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(dpop[,1], brain.eco[,3]),]
dpop[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])
dpop[c(103,366,367),1]<-c("Chiropotes_sagulatus", 
                             "Tursiops_sp.", "Tursiops_sp.")
rownames(dpop)<-paste0(dpop$Species, "_", dpop$Pop)

mat<-read.csv("Pred.csv")
dpop<-dpop[!is.na(match(dpop$Species,d1$Species)),colnames(mat)[-1]] #match to columns of mat
rownames(mat)<-mat$X

for (i in c(25:33)){
  dpop[,i]<-factor(as.numeric(dpop[,i]))
}
dpop<-as.data.frame(cbind(dpop[,c(1:10,15:25,27:29)], dpop[,c(11:14, 26, 30:33)]))

mat<-rbind(mat[c(1:10,15:25,27:29, 11:14, 26, 30:33),])
mat<-as.matrix(cbind(mat[,c(2:11,16:26,28:30, 12:15, 27, 31:34)]))

mice.pop<-mice(dpop, m=20, method=c(rep("pmm", 24),rep("logreg",4),
              "polyreg", "logreg", "polyreg","polyreg", "logreg"),
              predictorMatrix=mat, visitSequence="monotone", maxit=100, 
              print=TRUE)
save(mice.pop, file="mice.pop.robj")

#plot to inspect results
stripplot(mice.pop)

#check for convergence
plot(mice.pop)

#check imputated values
densityplot(mice.sp, scales = list(x = list(relation = "free")),
            layout = c(4,1))

#fill out the original data set with the all imputation values for each variable
comp.pops<-complete(mice.pop, "long")

comp.pop.cont<-ddply(comp.pops, .(.id), function(x) apply(x[,3:26], 
                                                         2, median))

comp.pop.disc<-ddply(comp.pops, .(.id), function(x) apply(x[,27:35], 2, function(y)
  names(sort(table(y), decreasing=T)[1])))

##change names for full version
datapop<-read.csv("FFpopulationsImp.csv")
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(datapop[,1], brain.eco[,3]),]
datapop[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])
datapop[c(103,366,367),1]<-c("Chiropotes_sagulatus", 
                         "Tursiops_sp.", "Tursiops_sp.")

datapop<-datapop[!is.na(match(datapop$Species,d1$Species)),] #match to columns of mat


d.pop<-data.frame(Species=datapop$Species, comp.pop.cont[,-1], 
                  comp.pop.disc[,-1])


dpop<-d.pop

save(dpop, file="DatapopImputed.robj")
write.csv(dpop, file="FFPopImputed.csv")








#impute dpop using the output of mice at the species level
nreps<-ddply(dpop,.(Species), function(x) apply(x[2:ncol(x)], 2, function(y) sum(is.na(y))))
b<-comp.pops[,c(".id", colnames(dpop)[-1])]

for (i in 2:ncol(dpop)){
  for (j in 1:nrow(nreps)){
    
    dpop[which(is.na(dpop[dpop$Species==nreps[j,1],i])),i]<-sample(b[which(!is.na(match(b$.id,nreps[j,1]))),i],
                                                                   nreps[j,2], replace=T)
  }
 
}
comp.pops




##load eco and brain size data, change species names of pops
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(dpop[,1], brain.eco[,3]),]
dpop[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])
dpop[56,1]<-"Chiropotes_albinasus"

#match names of pops with species
d1<-arrange(d1,Species)
d1$Species[c(2,4,22,29,48,56,64)]<-c("Antechinus_stuartii", "Antechinus_agilis",  
                                  "Capra_aegagrus","Chiropotes_sagulatus","Orcaella_heinsohni",
                                  "Physeter_macrocephalus", "Sotalia_guianensis")
dpop<-dpop[!is.na(match(dpop$Species,d1$Species)), ]


#verify pop species in the tree
b<-which(is.na(mammal.tree[[1]]$tip.label[match(d1[,1], mammal.tree[[1]]$tip.label)]))

##to replace pops species names
replace<-c(rep("Antechinus_leo",4), "Antechinus_bellus", 
           "Capra_caucasica","Capra_caucasica","Chiropotes_sagulatus", 
           "Orcaella_brevirostris", rep("Physeter_catodon",4), 
           rep("Sotalia_fluviatilis",6), rep("Tursiops_sp",2)) 
d1[b,1]<-replace
d2<-data.frame(species=unique(d1[,1]))
rownames(d2)<-d2[,1]


# Imputating with GLS - nice try :-0 --------------------------------------

phylovars<-d1[,c(1, 3,6,13, 16,18,21)] #select vars with phylo signal as tested by lambda.fit
lambda<-c(0,0.495,1,0.44,1,1,0.957)

dsp<-dsp[!is.na(dsp[,1]),colnames(d1)]
dsp.disc<-as.data.frame(dsp[,c(1,10:12,18,20:22)])
dsp<-as.data.frame(dsp[, c(1:9,13:17,19)])
write.csv(dsp.disc, "SpeciesDiscImpute.csv")

##impute continuous variables based on phylo distance using pGLM, only those with a high phylogenetic signal

tree<-prune.sample(t(dsp), mammal.tree[[1]])
big.fake<-rlnorm(nrow(dsp))
dsp$fake<-big.fake

c.data<-list()
pglm<-list()
pglm.null<-list()
aic<-list()

##fit gls to the data using a fake random variable and the phylogenetic  vcv
for (i in 2:(length(dsp)-1)){
  c.data[[i]]<-comparative.data(phy=tree, data=as.data.frame(dsp[,c(1,i,16)]), names.col='Species',
                                vcv=T, na.omit=F)
  pglm.null[[i]]<-pgls(log(c.data[[i]]$data[,1]) ~ c.data[[i]]$data[,2], data=c.data[[i]], lambda=0.000001)
  pglm[[i]]<-pgls(log(c.data[[i]]$data[,1]) ~ c.data[[i]]$data[,2], data=c.data[[i]], lambda="ML")
  aic[[i]]<-AIC(pglm[[i]], pglm.null[[i]])
}

predictNAs<-list() ##predict for the variables with high significant phylogenetic distance
fake2<-rlnorm(length(c.data[[2]]$phy$tip.label))
fake2<-as.data.frame(cbind(c.data[[2]]$phy$tip.label, fake=fake2))
rownames(fake2)<-fake2[,1]

for (i in c(2,6,10:14)){
  predictNAs[[i]]<-predict(pglm[[i]], data=fake2)
}

predicted<-data.frame(matrix(unlist(predictNAs), nrow=69, byrow=T))
colnames(predicted)<-colnames(dsp[,c(2,6,10:14)])
predicted<-cbind(Species=c.data[[2]]$phy$tip.label, predicted)
rownames(predicted)<-predicted[,1]

##replace missing values with predicted values
d1.pred<-ddply(d1[,colnames(predicted)], .(Species), function(x) apply(x[,-1], 2, function(y)
  sum(is.na(y))))

##generate replacement values
all.rep.predict<-list()
for (i in 2:length(predicted)){
  mod<-pglm[[i]]
  rep.predict<-as.data.frame(matrix(,nrow=nrow(fake2),ncol=max(d1.pred[,i])))
  for (j in 1:max(d1.pred[,i])){
    mod$model$coef[1,1]<-mean(rnorm(100,pglm[[i]]$model$coef[1,1], pglm[[i]]$sterr[1]*sqrt(pglm[[i]]$n)))
    rep.predict[,j]<-predict(mod, data=fake2)
  }
  rownames(rep.predict)<-c.data[[2]]$phy$tip.label
  all.rep.predict[[i]]<-rep.predict
}

#replace values
cols<-d1[,colnames(predicted)]

NAs<-list()
for (i in 2:length(cols)){
  NAs[[i]]<-cols[is.na(cols[,i]),c(1,i)]
  NA.freq<-ddply(NAs[[i]], .(Species), function(x) apply(x, 2, function(y) sum(is.na(y))))
  NA.freq[,1]<-sort(unique(NAs[[i]][,1]))
  
  for (j in 1:nrow(NA.freq)){
    print(j)
    NAs[[i]][NAs[[i]][,1]==NA.freq[j,1], 2]<-t(all.rep.predict[[i]][NA.freq[j,1],c(1:NA.freq[j,2])])
  }
  cols[rownames(NAs[[i]]),i]<-NAs[[i]][,2]
}

#import variable imputed with phylopars
phylo.imp<-read.table("PhyloparsSGrpSD.txt", sep="\t", header=T)
phylo.imp$X[phylo.imp$X==""]<-NA

##replace names changed because of phylogeny
sg_size_sd<-phylo.imp[!is.na(phylo.imp$X),c(1,3)]
sg_size_sd$X<-gsub(" ", "_", sg_size_sd$X)
sg_size_sd$X[ sg_size_sd$X %in% replace.sp[-8]] <-c("Capra_aegagrus","Orcaella_heinsohni","Sotalia_guianensis",
                                                    "Physeter_macrocephalus","Antechinus_subtropicus",
                                                    "Antechinus_agilis")                                              

sg_size_sd<-sg_size_sd[-1,]
C_sagulatus<-c(0, mean(c(13.9,10,14.4881)))
sg_size_sd<-rbind(sg_size_sd, C_sagulatus) 
sg_size_sd[72,1]<-"Chiropotes_sagulatus"
sg_size_sd<-arrange(sg_size_sd,X)
rownames(sg_size_sd)<-sg_size_sd$X

##replace sg sd in imputed d1 object
d1<-d.sp
brain.eco<-read.csv("BodyEcoTraits.csv")
brain.eco<-brain.eco[match(d1[,1], brain.eco[,3]),]
d1[,1]<-paste0(brain.eco[,1],"_", brain.eco[,2])   #change names of d1 species for full name
d1<-arrange(d1, Species)
d1$sg_size_sd<-sg_size_sd$sg_size_sd

save(d1, file="ImputedSpecies.robj")
