#this script is to analyse how the change in plant height affects aboveground biomass storage

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)

#read in data
Nit<-read.csv("Data/Nitrogen_trait_climate.csv")
head(Nit)
Nit<-subset(Nit, select = -c(Native_C3,Inv_C3) )
Nit<-Nit[complete.cases(Nit),]

Nit$Height_diff<-Nit$Inv_height-Nit$Native_height
Nit$Height_RR<-log(Nit$Inv_height)-log(Nit$Native_height)
Nit$Height_prop<-(Nit$Inv_height-Nit$Native_height)/Nit$Native_height
Nit$Woody_diff<-Nit$Inv_woodiness-Nit$Native_woodiness
Nit$Fix_diff<-Nit$Inv_Nit_fix-Nit$Nat_Nit_fix
Nit$Inv_fix<-ifelse(Nit$Inv_Nit_fix==1,"F","NF")
Nit$Diff_RR<-log(Nit$EF_I)-log(Nit$EF_UI)

#correct the variance for all measurements so that they are equal to SD
Nit$SE_UI<-ifelse(Nit$Var=="SE",Nit$SE_UI/sqrt(Nit$SS_UI),Nit$SE_UI)
Nit$SE_I<-ifelse(Nit$Var=="SE",Nit$SE_I/sqrt(Nit$SS_I),Nit$SE_I)


#do some data exploration
ggplot(Nit2,aes(x=Height_RR,y=Diff_RR,colour=Inv_fix))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Nit2,aes(x=CWD,y=Diff_RR,colour=Inv_fix))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Nit2,aes(x=Precip,y=Diff_RR,colour=Inv_fix))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Nit2,aes(x=Temp/10,y=Diff_RR,colour=Inv_fix))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Nit2,aes(x=Woody_diff,y=Diff_RR,colour=Inv_fix))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(Nit2,aes(x=Inv_fix,y=Diff_RR))+geom_point(shape=1)+ stat_summary(fun.y="mean", geom="point",size=2,colour="red")
ggplot(Nit_ES,aes(x=Height_RR,y=CWD,size=yi))+geom_point()


#now calculate effect sizes in metafor
Nit2<-subset(Nit,!is.na(Height_RR)&Terr_aqu=="Terrestrial")
Nit_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=Nit2)


M0<-rma.mv(yi~1,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M2<-rma.mv(yi~CWD,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M3<-rma.mv(yi~Precip,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M4<-rma.mv(yi~Temp,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M5<-rma.mv(yi~Height_RR*CWD,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M6<-rma.mv(yi~Height_RR*Precip,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M7<-rma.mv(yi~Height_RR*Temp,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M8<-rma.mv(yi~Height_RR*Inv_fix,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M9<-rma.mv(yi~Inv_fix,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M10<-rma.mv(yi~Height_RR+Inv_fix,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M11<-rma.mv(yi~Inv_fix*CWD,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M12<-rma.mv(yi~Inv_fix*Precip,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M13<-rma.mv(yi~Inv_fix*Temp,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M14<-rma.mv(yi~Inv_fix*Temp+Height_RR,vi,random=list(~1|Study),data=Nit_ES,method="ML")
M15<-rma.mv(yi~Height_RR+Inv_fix,vi,random=list(~1|Study),data=Nit_ES,method="ML")
AICc(M0,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15)

#get r squared value

1-(deviance(M8)/deviance(M0))

#create new dataset for predictions
new.data<-expand.grid(Inv_fix=c("F","NF"),
                      Height_RR=seq(min(Nit_ES$Height_RR,na.rm = T),max(Nit_ES$Height_RR,na.rm = T),length.out = 1000))


# Calculate yi as you did:

new.data$yi<-ifelse(new.data$Inv_fix=="NF",
                    0.6086-0.3303+(-2.1923*new.data$Height_RR)+(2.1944*new.data$Height_RR),0.6086+(-2.1923*new.data$Height_RR))
ifelse(new.data$Inv_fix="F"&)
# Plot
theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=yi,colour=Inv_fix))+geom_line()+geom_point(shape=1,data=Nit_ES,fill=NA,aes(size=1/vi))
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Difference in native and invasive nitrogen fixing ability")
ggsave("Figures/Nitrogen_traits.pdf",width = 8,height = 6,units = "in",dpi = 400)
