#script to scrape data on invasive plant species records
#from invasive species compendium and use this to see if there are any geographic or climate
#biases in the dataset we have used to assess impacts of invasive species on carbon pools

#author:Phil Martin
#date last edited:2016-04-20


#load packages
library(rvest)
library(ggmap)
library(xml2)
library(ggplot2)
library(plyr)
library(reshape)
library(raster)
library(gridExtra)
library(grid)
library(taxize)

#import list of all plant species that have datasheets in the ISC
ISC_records<-read.csv("Data/ISC_records.csv",stringsAsFactors = F)

#gather information on the location of species, in text format, from ISC
#but NOT their latitude and longitude, as this is not available
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

#loop through unique place names and return longitude and latitude data
un_loc<-as.character(unique(Geog_bias$Location))
geocode_summary<-NULL
for (i in 1:length(un_loc)){
  sub_loc<-geocode(un_loc[i],output = "more")
  sub_loc2<-data.frame(Location=un_loc[i],type=sub_loc$type,lon=sub_loc$lon,lat=sub_loc$lat)
  geocode_summary<-rbind(sub_loc2,geocode_summary)
}

#merge lat and longs with species data
Geom_bias_loc<-merge(Geog_bias,geocode_summary,by="Location")
#get count of number of invasive species for each unique long and lat
Geom_bias_summary<-ddply(Geom_bias_loc,.(Location,lon,lat),summarise,inv_count=length(Name))

#now map the locations of these
theme_set(theme_bw(base_size=12))
world<-map_data("world")
P_Bias1<- ggplot()+geom_polygon( data=world, aes(x=long, y=lat, group = group),fill="grey")
P_Bias2 <- P_Bias1+ geom_point(data=Geom_bias_summary,aes(x=lon,y=lat,size=inv_count),alpha=0.4,colour="blue")
P_Bias3<-P_Bias2+theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size=1.5,colour="black",fill=NA),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank())
P_Bias4<-P_Bias3+scale_size("Number of invasive plant \nspecies recorded")+coord_equal()+ylim(c(-55,90))+theme(legend.position="bottom")

#map locations of studies
studies<-read.csv("Data/Studies_climate.csv")
Studies_summary<-ddply(studies,.(long,lat),summarise,inv_count=length(EF_type))

P_Studies1<- ggplot()+geom_polygon( data=world, aes(x=long, y=lat, group = group),fill="grey")
P_Studies2 <- P_Studies1+ geom_point(data=Studies_summary,aes(x=long,y=lat,size=inv_count),alpha=0.4,colour="blue")
P_Studies3<-P_Studies2+theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_rect(size=1.5,colour="black",fill=NA),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank())
P_Studies4<-P_Studies3+scale_size("Number of datapoints \nused in study")+coord_equal()+ylim(c(-55,90))+theme(legend.position="bottom")

pdf("Figures/Bias_map.pdf")
grid.arrange(P_Bias4,P_Studies4)
dev.off()

#now look at climate space occupied by these species
CWD_raster<-raster("Data/Climate/CWD/cwd.bil")
Geom_bias_loc$CWD<-extract(CWD_raster,cbind(Geom_bias_loc$lon,Geom_bias_loc$lat),method='bilinear')
Inv_CWD<-ggplot(Geom_bias_loc,aes(x=CWD))+geom_histogram(aes(y=(..count../sum(..count..))*100))+xlim(c(-2500,0))
Inv_CWD2<-Inv_CWD+ylab("Percentage of locations with\n invasive plant records")+xlab("Climatic water deficit (mm)")

Study_CWD<-ggplot(studies,aes(x=CWD))+geom_histogram(aes(y=(..count../sum(..count..))*100))+xlim(c(-2500,0))
Study_CWD2<-Study_CWD+ylab("Percentage of datapoints \nused in study")+xlab("Climatic water deficit (mm)")

pdf("Figures/Bias_CWD.pdf")
grid.arrange(Inv_CWD2,Study_CWD2)
dev.off()

Geom_bias_loc[,10]
data.frame(CWD=studies[,4]
  


#now try putting all these figures together into one
pdf("Figures/Biases.pdf")
grid.arrange(P_Bias4,Inv_CWD2,P_Studies4,Study_CWD2,ncol=2)
dev.off()

#now a figure of just the difference in percentage of studies for different CWD values
Geom_bias_loc$CWD_bin<-cut(Geom_bias_loc$CWD,breaks = seq(-2500,0,by = 250),labels= seq(-2500,-250,by = 250),include.lowest=T)
Geom_bias_loc2<-Geom_bias_loc[complete.cases(Geom_bias_loc),]
head(Geom_bias_loc2)
Geom_bias_CWD_summary<-ddply(Geom_bias_loc2,.(CWD_bin),summarise,perc=(length(Name)/nrow(Geom_bias_loc2))*100)

studies$CWD_bin<-cut(studies$CWD,breaks = seq(-2500,0,by = 250),labels= seq(-2500,-250,by = 250),include.lowest=T)
studies2<-studies[complete.cases(studies),]
head(studies2)
studies_bias_CWD_summary<-ddply(studies2,.(CWD_bin),summarise,perc=(length(lat)/nrow(studies2))*100)

CWD_comp<-merge(studies_bias_CWD_summary,Geom_bias_CWD_summary,by="CWD_bin",all.y=T)
CWD_comp[is.na(CWD_comp)]<-0
CWD_comp$diff<-CWD_comp$perc.x-CWD_comp$perc.y

#plot pf differences in percentage of sites in different climatic zones
theme_set(theme_bw(base_size=12))
CWD_P1<-ggplot(CWD_comp,aes(x=CWD_bin,y=diff))+geom_point()+geom_hline(yintercept=0,lty=2)
CWD_P1+ylab("Difference in percentage of sites \nused in our study vs global records")+xlab("Climatic water deficit (mm)")
ggsave("Figures/CWD_diff_bias.pdf",width = 6,height=4,units = "in",dpi=400)
ggsave("Figures/CWD_diff_bias.png",width = 6,height=4,units = "in",dpi=400)


##################################################################
#use this section to look at taxonomic biases#####################
##################################################################

Tax_bias_summary<-ddply(Geom_bias_loc,.(Name),summarise,total=length(Location))


Corr_names<-gnr_resolve(names = Tax_bias_summary$Name,best_match_only=T)

Inv_species<-data.frame(Inv_sp=Corr_names$matched_name,total=Tax_bias_summary$total,Type="Global")

Studies<-read.csv("Data/Studies_combined.csv")
Study_species<-ddply(Studies,.(Inv_sp),summarise,total=length(Study))
Study_species$Type<-"Study"

unique(Study_species$Inv_sp)

Combined_species<-rbind(Study_species,Inv_species)

Try_gf<-read.csv("Data/TRY_growth_form.csv")

Species_gf<-merge(Combined_species,Try_gf,by.x="Inv_sp",by.y="AccSpeciesName",all.x=T)
Species_gf2<-subset(Species_gf,PlantGrowthForm!=""&PlantGrowthForm!="<NA>")
head(Species_gf2)
ddply(Species_gf2,.(Type),summarise,tot_sum=sum(total))
Species_gf2$perc<-ifelse(Species_gf2$Type=="Global",(Species_gf2$total/8602)*100,(Species_gf2$total/116)*100)
Species_gf3<-ddply(Species_gf2,.(Type,PlantGrowthForm),summarise,total=sum(perc))
Species_gf4<-rbind(Species_gf3,data.frame(Type="Study",PlantGrowthForm="fern",total=0))


theme_set(theme_bw(base_size=12))
GF_plot1<-ggplot(Species_gf4,aes(x=PlantGrowthForm,y=total,colour=Type,fill=Type))+geom_bar(stat="identity",width=0.6,position = position_dodge(width = 0.8))
GF_plot2<-GF_plot1+ylab("Percentage of locations (%)")+xlab("Plant growth form")+scale_color_brewer("Data set",palette="Set1")+scale_fill_brewer("Data set",palette="Set1")
GF_plot2+coord_cartesian(ylim=c(0,45))
ggsave("Figures/Growthform_diff_bias.pdf",width = 6,height=4,units = "in",dpi=400)
ggsave("Figures/Growthform_diff_bias.png",width = 6,height=4,units = "in",dpi=400)
