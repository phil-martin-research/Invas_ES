#this script is to analyse how the change in plant height affects aboveground biomass storage

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)

#read in data
AGB<-read.csv("Data/AGB_trait_climate.csv")
AGB<-AGB[complete.cases(AGB),]

AGB$Height_diff<-AGB$Inv_height-AGB$Native_height
AGB$Height_RR<-log(AGB$Inv_height)-log(AGB$Native_height)
AGB$Height_prop<-(AGB$Inv_height-AGB$Native_height)/AGB$Native_height
AGB$Woody_diff<-AGB$Inv_woodiness-AGB$Native_woodiness
AGB$Diff_RR<-log(AGB$EF_I)-log(AGB$EF_UI)
AGB$CWD2<-(AGB$CWD-mean(AGB$CWD))/sd(AGB$CWD)

#correct the variance for all measurements so that they are equal to SD
AGB$SE_UI<-ifelse(AGB$Var=="SE",AGB$SE_UI/sqrt(AGB$SS_UI),AGB$SE_UI)
AGB$SE_I<-ifelse(AGB$Var=="SE",AGB$SE_I/sqrt(AGB$SS_I),AGB$SE_I)


#do some data exploration
ggplot(AGB,aes(x=Height_diff,y=Diff_RR))+geom_point()+geom_smooth()
ggplot(AGB,aes(x=Height_RR,y=Diff_RR,label=Study))+geom_point()+geom_smooth(method="lm")
ggplot(AGB,aes(x=Height_prop,y=exp(Diff_RR)-1))+geom_point()
ggplot(AGB,aes(x=CWD,y=Diff_RR))+geom_point()
ggplot(AGB,aes(x=CWD2,y=Diff_RR))+geom_point()
ggplot(AGB,aes(x=Precip,y=Diff_RR))+geom_point()
ggplot(AGB,aes(x=Temp/10,y=Diff_RR))+geom_point()

#now calculate effect sizes in metafor
AGB2<-subset(AGB,!is.na(Height_RR)&Terr_aqu=="Terrestrial")
AGB_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=AGB2)


M0<-rma.mv(yi~1,vi,random=list(~1|Study),data=AGB_ES,method="ML")
M1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=AGB_ES,method="ML")
M2<-rma.mv(yi~CWD2,vi,random=list(~1|Study),data=AGB_ES,method="ML")
M3<-rma.mv(yi~Height_RR+CWD2,vi,random=list(~1|Study),data=AGB_ES,method="ML")
M4<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=AGB_ES,method="ML")
AICc(M0,M1,M2,M3,M4)

summary(M1)

#get r squared value

1-(deviance(M4)/deviance(M0))

#cerate new dataset for predictions
new.data<-expand.grid(Height_RR=seq(min(AGB_ES$Height_RR),max(AGB_ES$Height_RR),length.out = 100),
                      CWD2=seq(min(AGB_ES$CWD2,na.rm = T),max(AGB_ES$CWD2,na.rm = T),length.out = 100))
#create new dataframe to produce convex hull
AGB3<-AGB2[,c("Height_RR","CWD2")]

#get convex hull
ch <- chull(AGB3$Height_RR, AGB3$CWD2)
poly.df <- AGB3[c(ch, ch[1]),]
poly <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(poly.df[,1:2]))),1)))


# Create a SpatialPointsDataFrame with new.data:

new.data.poly <- new.data
coordinates(new.data.poly) <- ~Height_RR+CWD2

# Extract the points in new.data which are covered by the polygon:

new.data$inp <- over(new.data.poly, poly)
new.data <- new.data[complete.cases(new.data),]

# Calculate yi as you did:

new.data$yi<-(new.data$Height_RR*0.7048) + 0.5151 + (0.4241*new.data$CWD2) + ((new.data$CWD2*new.data$Height_RR)*-0.1593)

# Plot
theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=CWD2,fill=yi)) +
  geom_raster() + 
  scale_fill_gradient2("Change in community biomass \n(log response ratio)",low="blue",mid="light grey",high="red")+geom_point(shape=1,size=2,data=AGB2,fill=NA)
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Climatic water deficit (mm)")
ggsave("Figures/AGB_climate.pdf",width = 8,height = 4,units = "in",dpi = 400)
