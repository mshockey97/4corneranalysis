#PERMANOVA
library(vegan)
library(tidyverse)
library(dplyr)
library(readr)
library(vegan)
library(readxl)
#load plot data
envi <- read.csv("data/masterdata.csv")
spec=read_xlsx("data/sipseyspreadsheet2_frequency.xlsx")
#remove first col and add correct row names
spec=spec[,c(-1)]
rownames(spec)<-c(0:49)
#run permanova using Adonis
adonis(spec~Inundation.duration+Hand+BA+plotsd+plotele, data = envi)

#Making pariwise comparison

