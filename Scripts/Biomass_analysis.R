#this script is to analyse how the change in plant height affects aboveground biomass storage

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)

#read in data
AGB<-read.csv("Data/AGB_trait_climate.csv")
head(AGB)

AGB$Height_diff<-AGB$Inv_height-AGB$Native_height
AGB$Height_RR<-log(AGB$Inv_height)-log(AGB$Native_height)
AGB$Height_prop<-(AGB$Inv_height-AGB$Native_height)/AGB$Native_height
AGB$Woody_diff<-AGB$Inv_woodiness-AGB$Native_woodiness
AGB$Diff_RR<-log(AGB$EF_I)-log(AGB$EF_UI)
head(AGB)


#do some data exploration
ggplot(AGB,aes(x=Height_diff,y=Diff_RR))+geom_point()+geom_smooth()
ggplot(AGB,aes(x=Height_RR,y=Diff_RR,label=Study))+geom_point()+geom_smooth(method="lm")
ggplot(AGB,aes(x=Height_prop,y=exp(Diff_RR)-1))+geom_point()
ggplot(AGB,aes(x=CWD,y=Diff_RR))+geom_point()
ggplot(AGB,aes(x=Precip,y=Diff_RR))+geom_point()
ggplot(AGB,aes(x=Temp/10,y=Diff_RR))+geom_point()



#correct the variance for all measurements so that they are equal to SD
AGB$SE_UI<-ifelse(AGB$Var=="SE",AGB$SE_UI/sqrt(AGB$SS_UI),AGB$SE_UI)
AGB$SE_I<-ifelse(AGB$Var=="SE",AGB$SE_I/sqrt(AGB$SS_I),AGB$SE_I)
names(AGB)

#now calculate effect sizes in metafor
AGB_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=AGB)
AGB_ES2<-subset(AGB_ES,!is.na(Height_RR))


M0<-rma.mv(yi~1,vi,random=list(~1|Study),data=AGB_ES2,method="ML")
M1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=AGB_ES2,method="ML")
M2<-rma.mv(yi~Height_RR*CWD,vi,random=list(~1|Study),data=AGB_ES2,method="ML")
M3<-rma.mv(yi~Height_RR+Temp,vi,random=list(~1|Study),data=AGB_ES2,method="ML")
AICc(M0,M1,M2,M3)

#get r squared value


1-(deviance(M1)/deviance(M0))

new.data<-expand.grid(Height_RR=seq(min(AGB_ES2$Height_RR),max(AGB_ES2$Height_RR),length.out = 10),
                      CWD=seq(min(AGB_ES2$CWD,na.rm = T),max(AGB_ES2$CWD,na.rm = T),length.out = 10))
new.data$yi<-(new.data$Height_RR*0.3412)+1.1422+(0.0014*new.data$CWD)+((new.data$CWD*new.data$Height_RR)*0.0006)
new.data$SE<-0.0945+0.1409
new.data$UCI<-(new.data$Height_RR*0.7599)+0.4028
new.data$LCI<-(new.data$Height_RR*0.3896)-0.1497


#now plot this result
theme_set(theme_bw(base_size=12))
P1<-ggplot(AGB_ES2,aes(x=Height_RR,y=yi))+geom_point(shape=1,aes(size=1/vi))+geom_ribbon(data=new.data,aes(x=Height_RR,y=yi,ymax=yi+(2*SE),ymin=yi-(2*SE)),size=2,alpha=0.5)
P1<-ggplot(AGB_ES2,aes(x=Height_RR,y=yi))+geom_point(shape=1,aes(size=1/vi))+geom_ribbon(data=new.data,aes(x=Height_RR,y=yi,ymax=UCI,ymin=LCI,size=2,alpha=0.5))
P2<-P1+geom_line(data=new.data,aes(x=Height_RR,y=yi),size=1)+scale_size(range = c(2,6))+ guides(size=FALSE)
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Change in aboveground biomass \nfollowing invasion (log response ratio)")
ggsave("Figures/AGB_height.pdf",width = 6,height = 4,dpi = 400,units = "in")

P1<-ggplot(new.data,aes(x=Height_RR,y=CWD,fill=yi))+geom_raster()
ggplot(data=AGB_ES2,aes(x=Height_RR,y=CWD,colour=yi,size=yi))+geom_point()

