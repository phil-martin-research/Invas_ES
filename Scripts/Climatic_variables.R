#this is a script to map the locations of study sites and get an idea of the climatic variables at the sites

rm(list=ls())

#load packages
library(ggplot2)
library(raster)
library(plyr)

#load in data
loc<-read.csv("Data/ES_invasives_new.csv",stringsAsFactors=FALSE)
head(loc)
loc$Latitude<-as.numeric(loc$Latitude)
loc$Longitude<-as.numeric(loc$Longitude)
loc2<-as.matrix(cbind(loc$Longitude,loc$Latitude))

#load in climate data
#first mean temperature
Temp<-raster("Data/Climate/Bioclim/bio1.bil")
#now precipitation
Prec<-raster("Data/Climate/Bioclim/bio12.bil")
#now cwd
CWD<-raster("Data/Climate/CWD/CWD.bil")

#now extract data for each of these layers
loc$Temp<-extract(Temp,loc2,method="bilinear",buffer=1000,fun=mean,na.rm=T)
loc$Precip<-extract(Prec,loc2,method="bilinear",buffer=1000,fun=mean,na.rm=T)
loc$CWD<-extract(CWD,loc2,method="bilinear",buffer=1000,fun=mean,na.rm=T)

is.na(loc$CWD)

ggplot(loc,aes(x=Temp/10))+geom_histogram()+facet_wrap(~EF_type)
ggplot(loc,aes(x=Precip))+geom_histogram()+facet_wrap(~EF_type)
ggplot(loc,aes(x=CWD))+geom_histogram()+facet_wrap(~EF_type)

write.csv(loc,"Data/Inv_trait_climate.csv",row.names = FALSE)

AGB<-subset(loc,EF_type=="Aboveground biomass")
write.csv(AGB,"Data/AGB_trait_climate.csv",row.names = FALSE)

Nit<-subset(loc,EF_type=="Soil nitrogen")
write.csv(Nit,"Data/Nitrogen_trait_climate.csv",row.names = FALSE)

Moist<-subset(loc,EF_type=="Soil moisture")
write.csv(Moist,"Data/Moisture_trait_climate.csv",row.names = FALSE)

Water<-subset(loc,EF_type=="Water supply")
write.csv(Water,"Data/Moisture_trait_climate.csv",row.names = FALSE)

