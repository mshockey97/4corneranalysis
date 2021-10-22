#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Multievarient GLM and forth corner analysis (https://heather-grab.github.io/Entom-4940/rql.html)
#Coder: Matt shockey (mshockey97@gmail.com)
#Date: 10/6/2021
#Purpose: Determine why trees are where they are
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Table of Contents~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Setup workspace
# Step 2: Pull data
# Step 3: Split by transect/plot



#Log Notes~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Setup workspace-------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Clear workspace 
remove(list=ls())

#Gather libraries of interest
library(tidyverse)
library(vegan)
library(readxl)
library(forestmangr)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(mvabund)
library(lattice)


#read in data
#species code for later
spec_code=read.csv("data/species_code.csv")
spec_code=spec_code %>% mutate(AccSpeciesName=scientific.name)


#load R matrix R: the site by environment matrix that describes the environmental conditions at each site.
envi <- read.csv("data/masterdata.csv", row.names = 1)
envi=envi[,c(1,2,3,4,13,14)]



#load the Q matrix the species by traits matrix that describes the traits of each species (pull from data base)
my_data <- read.csv("data//trait.csv")
my_data=left_join(my_data,spec_code)
q=my_data %>% 
  select(code,TraitName,OrigValueStr) %>%  
  group_by(code,TraitName) %>% 
  #pivot q wider
  pivot_wider(names_from = TraitName, values_from=OrigValueStr,values_fn=mean)
#write to CSV
write_csv(q,"data//traits.csv")
traits=read.csv("data/traits.csv")
#Filter out traits we have data for
traits=traits %>% 
  select(code, Species.tolerance.to.waterlogging,Species.habitat.characterisation...Plant.requirement..soil.pH,
         Leaf.nitrogen.phosphorus..N.P..ratio,Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...petiole.excluded,
         Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded,Leaf.phosphorus..P..content.per.leaf.dry.mass,
         Leaf.nitrogen..N..content.per.leaf.dry.mass) %>% 
  na.omit() %>% 
  #clean up
  mutate(code=code,'waterlog tolerance'=Species.tolerance.to.waterlogging,soilpH=Species.habitat.characterisation...Plant.requirement..soil.pH,
         'NP Ratio g/cm3'=Leaf.nitrogen.phosphorus..N.P..ratio, 'LMA petiol undefined mm2 mg-1'=Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded,
         'LMA petiol excluded mm2 mg-1'=Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...petiole.excluded,
         'leaf N mg/g'= Leaf.nitrogen..N..content.per.leaf.dry.mass, 'leaf P mg/g'= Leaf.phosphorus..P..content.per.leaf.dry.mass, 
        .keep="unused")
write_csv(traits,"data//traitsfinal.csv")
traits=read.csv("data/traitsfinal.csv", row.names = 1)
#load L matrix: the site by species matrix that describes which species occur at each site
abd <- read_csv("data/nmds.csv",col_names = T)
rownames(abd)<-c(0:49)


#exploring the data---------------------------------------------
#look at abundance across the floodplain
Tree_spp=mvabund(abd)
rownames(Tree_spp)<-c(0:49)
plot(Tree_spp)


#manyglm Plot to see if there is a pattern/cone 
mod1 <- manyglm(Tree_spp ~ envi$Hand, family="poisson")
plot(mod1)
mod2= manyglm(Tree_spp ~ envi$Hand, family="negative_binomial")
plot(mod2)
mod3=manyglm(Tree_spp~envi$Hand+envi$`Dist.from.river`+envi$plotsd, family = "negative_binomial")
plot(mod3)




#hypothesis testing with Anova
anova(mod3, nBoot=99)
anova(mod2, p.uni="adjusted", nBoot=99)



#multivarient SDM (https://cran.r-project.org/web/packages/mvabund/mvabund.pdf)
sdm_fit=traitglm(abd,envi)
sdm_fit$fourth
#plot
a        = max( abs(sdm_fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.spp = levelplot(t(as.matrix(sdm_fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.spp)


#Fourth corner analysis--------------------------------

#forth corner analysis
traits=traits
abd4th=abd %>% 
  mutate(QUMI=QUMI.1) %>% 
  select(ACRU,CACA18,CACO15,DIVI5,LIST2,QULA3,QUNI,QUMI,ULAM) 
fit=traitglm(abd4th,envi,traits)
fit$fourth #print fourth corner terms
#plot residuals
plot(fit)
#Check anova
anova(fit, nBoot = 10)
#summary can give you the effects of indivudual traits by environment interactions
#this is slow use nBoot=10 for test but 1000 for final analysis
summary(fit, nBoot=10)
#create heatmap
a        = max( abs(fit$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(fit$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)


#choosing best model using LASSO (Least Absolute Shrinkage and Selection Operator)
ft1=traitglm(abd4th,envi,traits,method="glm1path")
ft1$fourth #notice LASSO penalty has shrunk many interactions to zero
#plot residuals
plot(ft1)
#plot heat map
a        = max( abs(ft1$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(ft1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)

