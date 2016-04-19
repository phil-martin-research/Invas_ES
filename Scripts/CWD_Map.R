#this is a script to create a map of climatic water deficit

rm(list=ls())

#load packages
library(ggplot2)
library(raster)
library(plyr)

#now load in data
CWD<-raster("Data/Climate/CWD/CWD.bil")

#convert the raster to points for plotting
map.p <- rasterToPoints(CWD)

#Make the points a dataframe for ggplot
df <- data.frame(map.p)
#Make appropriate column headings


#Now make the map
P1<-ggplot(data=df, aes(y=y, x=x))+geom_raster(aes(fill=CWD))+scale_fill_continuous("CWD (mm per year)",low="red",high="light grey")
P1+theme_bw()+coord_equal()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
                                  axis.title.x=element_blank(),panel.grid.major = element_blank(),axis.title.y=element_blank(),
                                  panel.grid.minor = element_blank(),legend.position = "right",legend.key = element_blank())
ggsave("Figures/CWD_map.png",dpi=600)
