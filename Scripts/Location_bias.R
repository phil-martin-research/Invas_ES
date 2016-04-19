library(rvest)
library(ggmap)
library(xml2)

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
  if (length(Location_sub)>0){
  Location2<-data.frame(Name=ISC_records[i,1],Location=Location_sub,Nat_Inv=Nat_Inv,Invasive=Invasive,record=i)
  Location2$Location<-gsub("-", "", Location2$Location, fixed = TRUE)
  Location3<-subset(Location2,Invasive=="Invasive")
  }
  if (nrow(Location3)>0){
  Lat_longs<-geocode(Location3$Location)
  Location3$Long<-Lat_longs$lon
  Location3$Lat<-Lat_longs$lat
  #row.names(Location3)<-make.names(row.names(Location3),unique=TRUE)
  Location_summary<-rbind(Location3,Location_summary)
  }
}


write.csv(Location_summary,file = "Data/Bias_analyses/Invasive_location.csv")
