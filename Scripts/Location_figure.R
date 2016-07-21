#this script produces a figure of the study locations used for our meta-analysis

library(ggplot2)

#import data
studies<-read.csv("Data/Studies_climate.csv")
studies2<-subset(studies,!is.na(CWD))
mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp<-ggplot()+mapWorld

loc_map<-mp+geom_point(data=studies2,aes(x=long,y=lat,colour=EF_type),shape=1)+facet_wrap(~EF_type,ncol=1)+scale_colour_brewer("Carbon pool type",palette = "Set1")
loc_map+coord_equal()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),panel.grid.major = element_blank(),axis.title.y=element_blank(),
                            panel.grid.minor = element_blank(),legend.position = "right",legend.key = element_blank())
ggsave("Figures/Study_location_map.png",dpi=600)

sum(ifelse(studies2$lat>0,1,0))/nrow(studies2)
