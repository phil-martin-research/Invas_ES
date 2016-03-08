#script to sort out taxonomy for invasive species work

library(taxize)
library(devtools)
library(traits)

#load in species data
Species<-read.csv("Data/Species_sorted2.csv")


Species_sub<-as.vector(Species)

Resolve<-NULL
for (i in 1:nrow(Species_sub)){
  Species_subset<-Species[i,1]
  Resolve_sub<-gnr_resolve(names = Species_subset,best_match_only=T)
  Resolve<-rbind(Resolve,Resolve_sub)
  print(i)
  }


write.csv(data.frame(Resolve),"Data/Species_resolved.csv")

#read .csv back in after 

Species_fixed<-read.csv("Data/Species_resolved_fixed.csv")
Americas<-data.frame(Region="america",where=c( "Continental US", "Alaska", "Canada", "Caribbean Territories", "Central Pacific Territories", "Hawaii", "Mexico"))
Europe<-data.frame(Region="europe",where=c(  "Albania", "Austria", "Azores", "Belgium", "Islas_Baleares", "Britain", "Bulgaria", "Corse", "Kriti", "Czechoslovakia", "Denmark", "Faroer", "Finland", "France", "Germany", "Greece", "Ireland", "Switzerland", "Netherlands", "Spain", "Hungary", "Iceland", "Italy", "Jugoslavia", "Portugal", "Norway", "Poland", "Romania", "USSR", "Sardegna", "Svalbard", "Sicilia", "Sweden", "Turkey", "USSR_Northern_Division", "USSR_Baltic_Division", "USSR_Central_Division", "USSR_South_western", "USSR_Krym", "USSRSouth_eastern_Division"))
Regions<-rbind(Americas,Europe)
Species_origin<-NULL

#loop through for each country in list and each species
for (i in 1:nrow(Species_fixed)){
  Sp_sub<-as.character(Species_fixed[i,3])
 for(y in 1:nrow(Regions)){
   Sp_origin<-is_native(Sp_sub, where=as.character(Regions[y,2]), region = as.character(Regions[y,1]))
   if (length(Sp_origin)!=1){
   Species_or_sub<-data.frame(Species=Sp_sub,Region=as.character(Regions[y,1]),Country=as.character(Regions[y,2]),Origin="Unknown")
   Species_origin<-rbind(Species_origin,Species_or_sub)
   }else{
   Species_or_sub<-data.frame(Species=Sp_sub,Region=as.character(Regions[y,1]),Country=as.character(Regions[y,2]),Origin=Sp_origin)
  Species_origin<-rbind(Species_origin,Species_or_sub) 
   }
 }
 }



write.csv(Species_origin,"Data/Species_origin.csv")




