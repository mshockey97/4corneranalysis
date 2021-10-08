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



#load R matrix R: the site by envoronment matrix that describes the envonmental conditions at each site.
envi <- read.csv("data/masterdata.csv", row.names = 1)
envi=envi[,c(1,2,3,4,13,14)]



#load the Q matrix the species by traits matrix that describes the traits of each species (pull from data base)

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
mod3=manyglm(Tree_spp~envi$Hand+envi$`Dist from river`+envi$plotsd, family = "negative_binomial")
plot(mod3)
mod4= manyglm(Tree_spp ~ envi$Hand*1000, family="gamma")
plot(mod4)



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

