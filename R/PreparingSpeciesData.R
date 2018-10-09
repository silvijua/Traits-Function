#######################################################################################
           ####Preparing Data Set Social Systems Mammals-- SPECIES ########
#######################################################################################

#estimating medians from pop dataset
d<-read.csv("FFpopulations.csv")

dsp<-ddply(d, .(Species), function(x) apply(x[,c(7:35)],2,
                                                function(y) median(y,na.rm=T))) %>%
  arrange(Species)

dsp.disc<-read.csv("FFSpDisc.csv")

dsp<-data.frame(dsp, dsp.disc[,c(3:11)])

dsp$Subgroup_composition<-as.character(dsp$Subgroup_composition)
dsp$FissionBehav<-as.character(dsp$FissionBehav)
dsp$FusionBehav<-as.character(dsp$FusionBehav)
  
dsp[dsp==""]  <- NA
save(dsp, file="dsp.robj")

d1<-dsp   #social data at species level

# Prepare data at species level -------------------------------------------

#count non-NAs by trait ((demo, interactions, community, stability) for each species
pop.var<-data.frame(filled=apply(d1[,2:ncol(d1)], 2, function(y) sum(!is.na(y))/75))

##eliminate vars with very few sp (<30%)

pop.var[pop.var[,1]<0.3,]

d1<-subset(d1, ,-c(social.budget,FF_sd, FF_min))


#count non-NAs by species 

d.nonNA<-apply(d1[,2:ncol(d1)], 1, function(y) sum(!is.na(y))/35)

summary.d<-data.frame(d1[,1], d.nonNA)

##eliminate species with =<30% vars
d1.temp<-d1[-which(summary.d[,2]<0.3), ]
d1<-d1.temp[-which(d1.temp$Species=="Tursiops_sp."),]
save(d1, file="FFSp.robj")

##Impute data (SKIP for raw data) - Code ImputatingData
#load(file="ImputedSpecies2.robj")

#improving input? select continuous variables and check correlation
nominal<-c(1,28:36)


d1.cont<-d1[,2:27]
correl<-cor(d1.cont[,1:ncol(d1.cont)], use="pairwise.complete.obs", method="pearson")
high.corr<-which(correl>=0.7, arr.ind=T) 
high.corr2<-which(correl<=-0.7, arr.ind=T)
high.corr<-cbind(high.corr, id=c(1:nrow(high.corr)))


#test high cors for statistical significance
cor.p<-list()
cont.d1<-d1.cont[,1:ncol(d1.cont)]
for (i in 1:nrow(high.corr)) {
  cor.p[[i]]<-cor.test(cont.d1[,high.corr[i,2]], cont.d1[,high.corr[i,1]], 
                       alternative="two.sided", method="pearson")
}

high.p<-list()
for (i in 1:length(cor.p)) {
  high.p[i]<-cor.p[[i]]$p.value<0.05  #select the non significant correlations
}

#count non-NAs by trait ((demo, interactions, community, stability) for each species
pop.var<-data.frame(filled=apply(d1[,2:ncol(d1)], 2, function(y) sum(!is.na(y))))


#exclude highly correlated variables (r>.7) and less filled ones (<30%)
cor.vars<-c(3,4,6,12,14,16,19,21:24)
d1<-d1[,-cor.vars]

#check emptyness for each species again
d.nonNA<-apply(d1[,2:ncol(d1)], 1, function(y) sum(!is.na(y))/24)
summary.d<-data.frame(d1[,1], d.nonNA)
d1.temp<-d1[-which(summary.d[,2]<0.3), ]

##remove variables not so informative because replicated in meaning (e.g. propfemale)
d1<-d1[,-3]
save(d1, file="FFSp.robj")


###center and scale variables 
#load("ImputedSpecies2.robj") ##load sp. data
#load(file="PopImputed.robj") ##load pop data

scaled.d1<-apply(d1[,2:15], 2, function(x) scale(x, center=T, scale=T))

d1scaled<-data.frame(Species=d1[,1],  scaled.d1, d1[,16:24])
save(dsp, d1, d1scaled, file="FFsp.robj")

### to reduce noise in the data, select for each type of variables, those with #################
### highest species representation                                             ##################

pop.var<-as.data.frame(cbind(Variable=row.names(pop.var), frequency=pop.var[,1]))
pop.var[pop.var==""]<-NA
pop.var.demo<-arrange(pop.var[1:7,], desc(frequency))
demo<-pop.var.demo[1:5,]

pop.var.interact<-arrange(pop.var[8:14,], desc(frequency))
interact<-pop.var.interact[c(1:2,4:6),]

pop.var.comm<-arrange(pop.var[15:25,], desc(frequency))
comm<-pop.var.comm[c(1:3,7,9,11),]

pop.var.stab<-arrange(pop.var[26:32,], desc(frequency))
stab<-pop.var.stab[1:5,]

#join all top 5 variables of each type in a single data frame and extract this columns from d1

vartop<-rbind(demo, interact, comm, stab)
d1f<-subset(d1f, select=as.character(vartop[,1]))

save(d1f, file="d1f.robj")
