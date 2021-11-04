#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Multievarient GLM and forth corner analysis (https://heather-grab.github.io/Entom-4940/rql.html)
#Coder: Matt shockey (mshockey97@gmail.com)
#Date: 10/6/2021
#Purpose: Determine why trees are where they are
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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
envi=envi[,c(2,3,4,13,14)]



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
plot(Tree_spp)

#manyglm Plot to see if there is a pattern/cone. Its best to use negative binomial
#family for count data. Binomial for presence absence and poisson for presence only. 
#because we have count data I will be using a negative binomial distribution 
#source: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00190.x
#too crazy
mod1=manyglm(Tree_spp~envi$Hand+envi$`Dist.from.river`+envi$plotsd+envi$plotele+envi$Inundation.duration, family = "negative_binomial")
plot(mod1)
#second best model
mod2= manyglm(Tree_spp ~ envi$Hand, family="negative_binomial")
plot(mod2)

mod3=manyglm(Tree_spp~envi$Hand+envi$`Dist.from.river`+envi$plotsd, family = "negative_binomial")
plot(mod3)

mod4=manyglm(Tree_spp~envi$Hand*envi$`Dist.from.river`*envi$plotsd, family = "negative_binomial")
plot(mod4)
#best model
mod5= manyglm(Tree_spp ~ envi$Inundation.duration, family="negative_binomial")
plot(mod5)

mod6=manyglm(Tree_spp ~ envi$Hand+envi$Inundation.duration, family="negative_binomial")
plot(mod6)

mod7=manyglm(Tree_spp ~ envi$Inundation.duration+envi$Hand, family="negative_binomial")
plot(mod7)
#not significant
mod8=manyglm(Tree_spp ~ envi$Dist.from.river, family="negative_binomial")
plot(mod8)
#not significant
mod9=manyglm(Tree_spp ~ envi$plotsd, family="negative_binomial")
plot(mod9)
#significant
mod10=manyglm(Tree_spp ~ envi$plotele, family="negative_binomial")
plot(mod10)
#plot ele significant HAND 0.09
mod11=manyglm(Tree_spp ~ envi$plotele+envi$Hand, family="negative_binomial")
plot(mod11)
#hand yes plot ele 0.32
mod12=manyglm(Tree_spp ~ envi$Hand+envi$plotele, family="negative_binomial")
plot(mod12)
#Third best model
mod13=manyglm(Tree_spp ~ envi$plotele+envi$Inundation.duration, family="negative_binomial")
plot(mod13)

mod14=manyglm(Tree_spp ~ envi$plotele+envi$Hand+envi$Inundation.duration, family="negative_binomial")
plot(mod14)
#checking how much one model explains responces more than another
#to compare the fits of two models, you can use the anova() 
#function with the regression objects as two separate arguments. The anova() function will take the model objects as 
#arguments, and return an ANOVA testing whether the more complex model is significantly better at capturing the data than #
#the simpler model. If the resulting p-value is sufficiently low (usually less than 0.05), we conclude that the more complex 
#model is significantly better than the simpler model, and thus favor the more complex model. If the p-value is not sufficiently
#low (usually greater than 0.05), we should favor the simpler model.
anova(mod2,mod6)
anova(mod2,mod7)
anova(mod2,mod4)
anova(mod5,mod6)
anova(mod5,mod7)
anova(mod5,mod4)

#check models using predict function
modcheck1=predict(mod1, type = "response")
check1=abd-modcheck1
check1=round(check1,digits = 0)
summary(check1)
sum(AIC(mod1))
#HAND (2nd best model)
modcheck2=predict(mod2, type = "response")
check2=abd-modcheck2
check2=round(check2,digits = 0)
summary(check2)
sum(AIC(mod2))


modcheck3=predict(mod3, type = "response")
check3=abd-modcheck3
check3=round(check3,digits = 0)
summary(check3)
sum(AIC(mod3))


modcheck4=predict(mod4, type = "response")
check4=abd-modcheck4
check4=round(check4,digits = 0)
summary(check4)
sum(AIC(mod4))
#inundation duration (best model)
modcheck5=predict(mod5, type = "response")
check5=abd-modcheck5
check5=round(check5,digits = 0)
summary(check5)
sum(AIC(mod5))

modcheck6=predict(mod6, type = "response")
check6=abd-modcheck6
check6=round(check6,digits = 0)
summary(check6)
sum(AIC(mod6))
#plot ele (third best model)
modcheck10=predict(mod10, type = "response")
check10=abd-modcheck10
check10=round(check10,digits = 0)
summary(check10)
sum(AIC(mod10))
#hand and plot ele both of which are significant (until adjusted)
modcheck11=predict(mod11, type = "response")
check11=abd-modcheck11
check11=round(check11,digits = 0)
summary(check11)
sum(AIC(mod11))
#plot ele and inundation duration better than HAND and both are signifigant even when adjusted
modcheck13=predict(mod13, type = "response")
check13=abd-modcheck13
check13=round(check13,digits = 0)
summary(check13)
sum(AIC(mod13))
anova(mod5,mod13)

modcheck14=predict(mod14, type = "response")
check14=abd-modcheck14
check14=round(check14,digits = 0)
summary(check14)
sum(AIC(mod14))
anova(mod5,mod14)
#hypothesis testing with Anova
anova(mod1,p.uni="adjusted", nBoot=99)
summary(mod1)


anova(mod2, p.uni="adjusted", nBoot=99)
summary(mod2)

anova(mod3,p.uni="adjusted", nBoot=99)
summary(mod3)

anova(mod5, p.uni="adjusted", nBoot=99)
summary(mod5)

anova(mod6,p.uni="adjusted", nBoot=99)
summary(mod6)

anova(mod10,p.uni="adjusted", nBoot=99)
summary(mod10)

anova(mod11,p.uni="adjusted", nBoot=99)
summary(mod11)

anova(mod12,p.uni="adjusted", nBoot=99)
summary(mod12)

anova(mod13,p.uni="adjusted", nBoot=99)
summary(mod13)

anova(mod14,p.uni="adjusted", nBoot=99)
summary(mod14)
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
anova.traitglm(sdm_fit)
#multivarient sdm adding lasso penalty
sdm_fitl=traitglm(abd,envi,method = "glm1path")
sdm_fitl$fourth
#plot
a        = max( abs(sdm_fitl$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.spp = levelplot(t(as.matrix(sdm_fitl$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.spp)
anova.traitglm(sdm_fitl)
#Fourth corner analysis--------------------------------

#forth corner analysis
traits=traits
abd4th=abd %>% 
  mutate(QUMI=QUMI.1) %>% 
  select(ACRU,CACA18,CACO15,DIVI5,LIST2,QULA3,QUNI,QUMI,ULAM) 
fit=traitglm(abd4th,envi,traits)
fit$fourth 
#print fourth corner terms
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

anova.traitglm(fit)

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

