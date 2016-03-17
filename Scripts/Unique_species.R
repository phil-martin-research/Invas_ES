#script to produce a list of unique species names for invasive work

library(plyr)
library(dplyr)
library(reshape)

Species<-read.csv("Data/Species2.csv")
as.vector(Species)

Species_melt<-melt(Species,id.vars = "ID")

write.csv(unique(Species_melt$value),file = "Data/Species2_check.csv")
