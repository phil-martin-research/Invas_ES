library(rvest)
library(ggmap)
library(xml2)
library(ggplot2)
library(plyr)
library(reshape)

ISC_records<-read.csv("Data/ISC_records.csv",stringsAsFactors = F)

#create a loop to gather information on the location of species from ISC
#but NOT their latitude and longitude
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
  if (length(Location_sub)>0){
  Location2<-data.frame(Name=ISC_records[i,1],Location=Location_sub,Nat_Inv=Nat_Inv,Invasive=Invasive,record=i)
  Location2$Location<-gsub("-", "", Location2$Location, fixed = TRUE)
  Location3<-subset(Location2,Invasive=="Invasive")
  }
  Location_summary<-rbind(Location3,Location_summary)
  print(paste(round((i/nrow(ISC_records))*100,2),"% finished"))
}


head(Location_summary)

write.csv(Location_summary,file = "Data/Bias_analyses/Invasive_location.csv")


#for this section I need to produce a vector with all the unique locations collected from ISc
#and use this to geocode the location to return max and min lat and long values for the 
#region/country

Geog_bias<-read.csv("Data/Bias_analyses/Invasive_location.csv",stringsAsFactors = F)
head(Geog_bias)
#replace problematic place names
Geog_bias$Location<-ifelse(Geog_bias$Location=="Czechoslovakia (former)","Czech republic",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Yugoslavia (former)","Yugoslavia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Georgia (Republic of)","Rebublic of Georgia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Western Siberia","Siberia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Congo Democratic Republic","Democratic Republic of Congo",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="BosniaHercegovina","Bosnia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Micronesia, Federated states of","Federated states of Micronesia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Micronesia, Federated states of","Federated states of Micronesia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Sulawesi","Celebes",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Nei Menggu","Inner Mongolia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Former USSR","Russia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="GuineaBissau","Guinea-Bissau",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="Eastern Siberia","Siberia",Geog_bias$Location)
Geog_bias$Location<-ifelse(Geog_bias$Location=="TurkeyinAsia","Turkey",Geog_bias$Location)


un_loc<-as.character(unique(Geog_bias$Location))
geocode_summary<-NULL
for (i in 1:length(un_loc)){
  sub_loc<-geocode(un_loc[i],output = "more")
  sub_loc2<-data.frame(Location=un_loc[i],type=sub_loc$type,lon=sub_loc$lon,lat=sub_loc$lat,
             north=sub_loc$north,south=sub_loc$south,east=sub_loc$east,west=sub_loc$west)
  data.frame(Location=un_loc[i],type=sub_loc$type,lon=sub_loc$lon,lat=sub_loc$lat,
             lon_poly=rep(c(sub_loc$east,sub_loc$west),each=2),
             lat_poly=rep(c(sub_loc$north,sub_loc$south),2))
  geocode_summary<-rbind(sub_loc2,geocode_summary)
}

geocode("Guinea-Bissau",output = "more")

ggplot(geocode_summary,aes())

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

