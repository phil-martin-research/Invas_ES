#this script is to analyse how the change in plant height affects aboveground biomass storage

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)

#read in data
Nit<-read.csv("Data/Nitrogen_trait_climate.csv")
Nit<-Nit[complete.cases(Nit),]

Nit$Height_diff<-Nit$Inv_height-Nit$Native_height
Nit$Height_RR<-log(Nit$Inv_height)-log(Nit$Native_height)
Nit$Height_prop<-(Nit$Inv_height-Nit$Native_height)/Nit$Native_height
Nit$Woody_diff<-Nit$Inv_woodiness-Nit$Native_woodiness
Nit$Fix_diff<-Nit$Inv_Nit_fix-Nit$Nat_Nit_fix

Nit$Diff_RR<-log(Nit$EF_I)-log(Nit$EF_UI)
#correct the variance for all measurements so that they are equal to SD
Nit$SE_UI<-ifelse(Nit$Var=="SE",Nit$SE_UI/sqrt(Nit$SS_UI),Nit$SE_UI)
Nit$SE_I<-ifelse(Nit$Var=="SE",Nit$SE_I/sqrt(Nit$SS_I),Nit$SE_I)


#do some data exploration
ggplot(Nit,aes(x=Height_diff,y=Diff_RR))+geom_point()+geom_smooth()
ggplot(Nit,aes(x=Height_RR,y=Diff_RR,label=Study))+geom_point()+geom_smooth(method="lm")
ggplot(Nit,aes(x=Height_prop,y=exp(Diff_RR)-1))+geom_point()
ggplot(Nit,aes(x=CWD,y=Diff_RR))+geom_point()
ggplot(Nit,aes(x=Precip,y=Diff_RR))+geom_point()
ggplot(Nit,aes(x=Temp/10,y=Diff_RR))+geom_point()
ggplot(Nit,aes(x=Woody_diff,y=Diff_RR))+geom_jitter()
ggplot(Nit,aes(x=Fix_diff,y=Diff_RR))+geom_point()

#now calculate effect sizes in metafor
Nit2<-subset(Nit,!is.na(Height_RR)&Terr_aqu=="Terrestrial"&Fix_diff>-1)
Nit_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=Nit2)


M0<-rma.mv(yi~1,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M2<-rma.mv(yi~CWD,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M3<-rma.mv(yi~Precip,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M4<-rma.mv(yi~Temp,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M5<-rma.mv(yi~Height_RR*CWD,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M6<-rma.mv(yi~Height_RR*Precip,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M7<-rma.mv(yi~Height_RR*Temp,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M8<-rma.mv(yi~Height_RR*Fix_diff,vi,random=list(~1|Study),data=Nit_ES,method="ML")

AICc(M0,M1,M2,M3,M4,M5,M6,M7,M8)

#get r squared value

1-(deviance(M8)/deviance(M0))

#cerate new dataset for predictions
new.data<-expand.grid(Height_RR=seq(min(Nit_ES$Height_RR),max(Nit_ES$Height_RR),length.out = 1000),
                      Fix_diff=seq(min(Nit_ES$Fix_diff,na.rm = T),max(Nit_ES$Fix_diff,na.rm = T),length.out = 1000))
#create new dataframe to produce convex hull
Nit3<-Nit2[,c("Height_RR","Fix_diff")]

#get convex hull
ch <- chull(Nit3$Height_RR, Nit3$Fix_diff)
poly.df <- Nit3[c(ch, ch[1]),]
poly <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(poly.df[,1:2]))),1)))


# Create a SpatialPointsDataFrame with new.data:

new.data.poly <- new.data
coordinates(new.data.poly) <- ~Height_RR+Fix_diff

# Extract the points in new.data which are covered by the polygon:

new.data$inp <- over(new.data.poly, poly)
new.data <- new.data[complete.cases(new.data),]

# Calculate yi as you did:

new.data$yi<-(new.data$Height_RR*0.0001) + 0.3616 + (0.3200*new.data$Fix_diff) + ((new.data$Fix_diff*new.data$Height_RR)*-2.1919)

# Plot
theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=Fix_diff,fill=yi)) +
  geom_raster() + 
  scale_fill_gradient2("Change in soil nitrogen \n(log response ratio)",low="blue",mid="light grey",high="red")+geom_point(shape=1,size=2,data=Nit2,fill=NA)
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Difference in native and invasive nitrogen fixing ability")
ggsave("Figures/Nitrogen_traits.pdf",width = 8,height = 6,units = "in",dpi = 400)
