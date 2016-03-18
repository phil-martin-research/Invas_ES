#this script is to analyse how the change in plant height affects soil carbon

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)

#read in data
Carb<-read.csv("Data/SC_trait_climate.csv")
head(Carb)
Carb<-subset(Carb, select = -c(Native_C3,Inv_C3) )
Carb<-subset(Carb, Terr_aqu=="Terrestrial" )
Carb<-Carb[complete.cases(Carb),]

Carb$Height_diff<-Carb$Inv_height-Carb$Native_height
Carb$Height_RR<-log(Carb$Inv_height)-log(Carb$Native_height)
Carb$Height_prop<-(Carb$Inv_height-Carb$Native_height)/Carb$Native_height
Carb$Woody_diff<-Carb$Inv_woodiness-Carb$Native_woodiness
Carb$Diff_RR<-log(Carb$EF_I)-log(Carb$EF_UI)
Carb$CWD2<-(Carb$CWD-mean(Carb$CWD))/sd(Carb$CWD)

#correct the variance for all measurements so that they are equal to SD
Carb$SE_UI<-ifelse(Carb$Var=="SE",Carb$SE_UI/sqrt(Carb$SS_UI),Carb$SE_UI)
Carb$SE_I<-ifelse(Carb$Var=="SE",Carb$SE_I/sqrt(Carb$SS_I),Carb$SE_I)


#do some data exploration
ggplot(Carb,aes(x=Height_RR,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Carb,aes(x=CWD,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Carb,aes(x=Precip,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Carb,aes(x=Temp/10,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Carb,aes(x=Woody_diff,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Carb,aes(x=Height_RR,y=CWD,size=Diff_RR))+geom_point()


#now calculate effect sizes in metafor
Carb_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=Carb)


M0<-rma.mv(yi~1,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M2<-rma.mv(yi~CWD2,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M3<-rma.mv(yi~Precip,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M4<-rma.mv(yi~Temp,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M5<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M6<-rma.mv(yi~Height_RR*Precip,vi,random=list(~1|Study),data=Carb_ES,method="ML")
M7<-rma.mv(yi~Height_RR*Temp,vi,random=list(~1|Study),data=Carb_ES,method="ML")
AICc(M0,M1,M2,M3,M4,M5,M6,M7)

#get r squared value

1-(deviance(M5)/deviance(M0))

#create new dataset for predictions
new.data<-expand.grid(Height_RR=seq(min(Carb_ES$Height_RR),max(Carb_ES$Height_RR),length.out = 1000),
                      CWD2=seq(min(Carb_ES$CWD2,na.rm = T),max(Carb_ES$CWD2,na.rm = T),length.out = 1000))
#create new dataframe to produce convex hull
Carb3<-Carb[,c("Height_RR","CWD2")]

#get convex hull
ch <- chull(Carb3$Height_RR, Carb3$CWD2)
poly.df <- Carb3[c(ch, ch[1]),]
poly <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(poly.df[,1:2]))),1)))


# Create a SpatialPointsDataFrame with new.data:
new.data.poly <- new.data
coordinates(new.data.poly) <- ~Height_RR+CWD2

# Extract the points in new.data which are covered by the polygon:
new.data$inp <- over(new.data.poly, poly)
new.data <- new.data[complete.cases(new.data),]

# Calculate yi as you did:
new.data$yi<-(new.data$Height_RR*0.1052) + 0.1943 + (1.2867*new.data$CWD2) + ((new.data$CWD2*new.data$Height_RR)*0.0341)

# Plot
sd(Carb$CWD)
mean(Carb$CWD)

theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=(CWD2*374.6319)+-551.739,fill=yi)) +
  geom_raster() + 
  scale_fill_gradient2("Change in soil carbon storage \n(log response ratio)",low="blue",mid="light grey",high="red")+geom_point(shape=1,size=2,data=Carb,fill=NA)
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Climatic water deficit (mm)")
ggsave("Figures/Carb_climate.pdf",width = 8,height = 4,units = "in",dpi = 400)
