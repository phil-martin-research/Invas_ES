#this is a script to map the locations of study sites and get an idea of the climatic variables at the sites

rm(list=ls())

#load packages
library(ggplot2)
library(raster)
library(plyr)

#load in data
loc<-read.csv("Data/Invasive_locations.csv",stringsAsFactors=FALSE)
loc2<-unique(loc)
str(loc2)
loc2$Latitude<-as.numeric(loc2$Latitude)
loc2$Longitude<-as.numeric(loc2$Longitude)
loc3<-as.matrix(cbind(loc2[3],loc2[2]))

#load in climate data
#first mean temperature
Temp<-raster("Data/Climate/Bioclim/bio1.bil")
#now precipitation
Prec<-raster("Data/Climate/Bioclim/bio12.bil")
#now cwd
CWD<-raster("Data/Climate/CWD/CWD.bil")

#now extract data for each of these layers
loc2$Temp<-extract(Temp,loc3,method="bilinear")
loc2$Precip<-extract(Prec,loc3,method="bilinear")
loc2$CWD<-extract(CWD,loc3,method="bilinear")




ggplot(loc2,aes(x=Temp/10))+geom_histogram()
ggplot(loc2,aes(x=Precip))+geom_histogram()
ggplot(loc2,aes(x=CWD))+geom_histogram()

