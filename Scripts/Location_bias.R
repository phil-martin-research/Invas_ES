library(rvest)
library(ggmap)
library(xml2)
library(ggplot2)
library(plyr)
library(reshape)

ISC_records<-read.csv("Data/ISC_records.csv",stringsAsFactors = F)

Location_summary<-NULL
for (i in 1:nrow(ISC_records)){
  isc <- read_html(ISC_records[i,4])
  ns <- xml_ns(isc)
  Location<-xml_text(xml_find_all(isc, xpath="//div[@id='toDistributionTable']/table/tbody/tr/td[1]", ns))
  Nat_Inv<-xml_text(xml_find_all(isc, xpath="//div[@id='toDistributionTable']/table/tbody/tr/td[4]", ns))
  Invasive<-xml_text(xml_find_all(isc, xpath="//div[@id='toDistributionTable']/table/tbody/tr/td[6]", ns))
  Location_sub<-subset(Location,Location!="ASIA"&Location!="AFRICA"&Location!="NORTH AMERICA"&
                         Location!="CENTRAL AMERICA AND CARIBBEAN"&Location!="SOUTH AMERICA"&
                         Location!="OCEANIA"&Location!="EUROPE"&Location!="ANTARCTICA"&Location!="SEA AREAS")
  Location_sub<-ifelse(Location_sub=="Georgia (Republic of)","Republic of Georgia",Location_sub)
  Location_sub<-ifelse(Location_sub=="-Nei Menggu","Inner Mongolia",Location_sub)
  Location_sub<-ifelse(Location_sub=="Turkey-in-Asia","Turkey",Location_sub)
  if (length(Location_sub)>0){
  Location2<-data.frame(Name=ISC_records[i,1],Location=Location_sub,Nat_Inv=Nat_Inv,Invasive=Invasive,record=i)
  Location2$Location<-gsub("-", "", Location2$Location, fixed = TRUE)
  Location3<-subset(Location2,Invasive=="Invasive")
  Sub_loc<-NULL
  }
  if (nrow(Location3)>0){
  for (j in 1:nrow(Location3)){
  Sub_loc3<-Location3[j,]
  Lat_longs<-geocode(Sub_loc3$Location,output = "more",source = "google")
  Sub_loc3$long<-Lat_longs$lon
  Sub_loc3$lat<-Lat_longs$lat
  Sub_loc3$north<-Lat_longs$north
  Sub_loc3$south<-Lat_longs$south
  Sub_loc3$east<-Lat_longs$east
  Sub_loc3$west<-Lat_longs$west
  Sub_loc3$country<-Lat_longs$country
  Sub_loc<-rbind(Sub_loc3,Sub_loc)
  }
  Location_summary<-rbind(Sub_loc,Location_summary)
  }
  print(paste(round((i/nrow(ISC_records))*100,2),"% finished"))
}


write.csv(Location_summary,file = "Data/Bias_analyses/Invasive_location.csv")


###############################################################################
#this section allows mapping of where invasive species occur###################
#with data taken from the Invasive Species Compendium##########################
###############################################################################

Geog_bias<-read.csv("Data/Bias_analyses/Invasive_location.csv")
head(Geog_bias)

Geom_summ<-ddply(Geog_bias,.(Long,Lat),summarise,Inv_count=length(Invasive))

#now map the locations of these
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
P_Bias1<- ggplot() +   mapWorld
P_Bias2 <- P_Bias1+ geom_point(data=Geom_summ,aes(x=Long, y=Lat,size=Inv_count) ,color="blue") 
P_Bias2

ggplot()+geom_point(data=Geom_summ,aes(x=Long, y=Lat,size=Inv_count) ,color="blue") 

