#script to look at changes in ecosystem services following non-native plant invasions
#looking specifically at differences in the impacts of different functional types

rm(list=ls())

#first load packages
library(ggplot2)
library(MuMIn)
library(metafor)
library(reshape)
library(plyr)


#get data
setwd("C:/Users/Phil/Dropbox/Work/Active projects/PhD/Publications, Reports and Responsibilities/Chapters/2. Invasive species meta-analysis/Invas_ES/Data")
Inv<-read.csv("Invasive_ES.csv")

#correct life form values
levels(Inv$NS_LF)[which(levels(Inv$NS_LF)=="Sedge")] <- "Forb"
levels(Inv$NS_LF)[which(levels(Inv$NS_LF)=="Forbs")] <- "Forb"

#create woody/not woody column
#for invasives
Inv$WI<-Inv$IS_LF
levels(Inv$WI)[which(levels(Inv$WI)=="Forb")] <- "No"
levels(Inv$WI)[which(levels(Inv$WI)=="Grass")] <- "No"
levels(Inv$WI)[which(levels(Inv$WI)=="Tree")] <- "Yes"
levels(Inv$WI)[which(levels(Inv$WI)=="Shrub")] <- "Yes"

#and for natives
Inv$WN<-Inv$NS_LF
levels(Inv$WN)[which(levels(Inv$WN)=="Forb")] <- "No"
levels(Inv$WN)[which(levels(Inv$WN)=="Grass")] <- "No"
levels(Inv$WN)[which(levels(Inv$WN)=="Tree")] <- "Yes"
levels(Inv$WN)[which(levels(Inv$WN)=="Shrub")] <- "Yes"


#create a loop to compare invader and native species
Inv$Trans<-NULL
for (i in 1:nrow(Inv)){
  Inv$Trans[i]<-paste(as.character(Inv$NS_LF[i]),"-",as.character(Inv$IS_LF[i]),sep ="")
}

Inv$Trans<-ifelse(is.na(Inv$NS_LF)|is.na(Inv$IS_LF),NA,Inv$Trans)

#now create a loop for a new classification of ecosystem service type
levels(Inv$Serv_type)[which(levels(Inv$Serv_type)=="Water infiltration")] <- "Water provision"
levels(Inv$Serv_type)[which(levels(Inv$Serv_type)=="Evapotranspiration")] <- "Water provision"


#subset into different service types
#first aboveground carbon
AGB<-subset(Inv,Serv_type=="Aboveground carbon storage")
AGB<-AGB[AGB$WN %in% c("Yes","No"),]
AGB$WN<-factor(AGB$WN)

AGB_RR<-escalc(data=AGB,measure="ROM",m2i=ES_UI,sd2i=SE_UI,n2i=SS_UI,m1i=ES_I,sd1i=SE_I,n1i=SS_I,append=T)

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|SiteID),method="ML",data=AGB_RR)
Model1<-rma.mv(yi,vi,mods=~WI,random=list(~1|SiteID),method="ML",data=AGB_RR)
Model2<-rma.mv(yi,vi,mods=~WN,random=list(~1|SiteID),method="ML",data=AGB_RR)
Model3<-rma.mv(yi,vi,mods=~WN+WI,random=list(~1|SiteID),method="ML",data=AGB_RR)
Model4<-rma.mv(yi,vi,mods=~WN*WI,random=list(~1|SiteID),method="ML",data=AGB_RR)
Model5<-rma.mv(yi,vi,mods=~Trans,random=list(~1|SiteID),method="ML",data=AGB_RR)

AICc(Model0,Model1,Model2,Model3,Model4,Model5)

Model_pred<-data.frame(Pred=predict(Model3)$pred,ci.lb=predict(Model3)$ci.lb,ci.ub=predict(Model3)$ci.ub)
Model_pred$WI<-AGB_RR$WI
Model_pred$WN<-AGB_RR$WN


Model_pred2<-Model_pred[!duplicated(Model_pred), ]
head(Model_pred2)

ggplot(Model_pred,aes(x=WI,y=exp(Pred)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1,colour=WN))+geom_pointrange(position=position_dodge(width=0.5))+geom_hline(y=0,lty=2)



#next belowground carbon
BGB<-subset(Inv,Serv_type=="Belowground carbon storage")
BGB<-BGB[BGB$WN %in% c("Yes","No"),]
BGB$WN<-factor(BGB$WN)

BGB_RR<-escalc(data=BGB,measure="ROM",m2i=ES_UI,sd2i=SE_UI,n2i=SS_UI,m1i=ES_I,sd1i=SE_I,n1i=SS_I,append=T)
BGB_RR<-subset(BGB_RR,!is.na(yi))

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|SiteID),method="ML",data=BGB_RR)
Model1<-rma.mv(yi,vi,mods=~WI,random=list(~1|SiteID),method="ML",data=BGB_RR)
Model2<-rma.mv(yi,vi,mods=~WN,random=list(~1|SiteID),method="ML",data=BGB_RR)
Model3<-rma.mv(yi,vi,mods=~WN+WI,random=list(~1|SiteID),method="ML",data=BGB_RR)
Model4<-rma.mv(yi,vi,mods=~WN*WI,random=list(~1|SiteID),method="ML",data=BGB_RR)
Model5<-rma.mv(yi,vi,mods=~Trans,random=list(~1|SiteID),method="ML",data=BGB_RR)

AICc(Model0,Model1,Model2,Model3,Model4,Model5)

Model_pred<-data.frame(Pred=predict(Model1)$pred,ci.lb=predict(Model1)$ci.lb,ci.ub=predict(Model1)$ci.ub)
Model_pred$WI<-BGB_RR$WI



Model_pred2<-Model_pred[!duplicated(Model_pred), ]
head(Model_pred2)

ggplot(Model_pred,aes(x=WI,y=exp(Pred)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1))+geom_pointrange(position=position_dodge(width=0.5))+geom_hline(y=0,lty=2)

#next belowground carbon
CC<-subset(Inv,Serv_type=="Carbon sequestration")
CC<-CC[CC$WN %in% c("Yes","No"),]
CC$WN<-factor(CC$WN)
CC$NS_LF<-factor(CC$NS_LF)

CC_RR<-escalc(data=CC,measure="ROM",m2i=ES_UI,sd2i=SE_UI,n2i=SS_UI,m1i=ES_I,sd1i=SE_I,n1i=SS_I,append=T)
CC_RR<-subset(CC_RR,!is.na(yi))

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|SiteID),method="ML",data=CC_RR)
Model1<-rma.mv(yi,vi,mods=~WI,random=list(~1|SiteID),method="ML",data=CC_RR)
Model2<-rma.mv(yi,vi,mods=~WN,random=list(~1|SiteID),method="ML",data=CC_RR)
Model3<-rma.mv(yi,vi,mods=~WN+WI,random=list(~1|SiteID),method="ML",data=CC_RR)
Model4<-rma.mv(yi,vi,mods=~IS_LF,random=list(~1|SiteID),method="ML",data=CC_RR)
Model5<-rma.mv(yi,vi,mods=~NS_LF,random=list(~1|SiteID),method="ML",data=CC_RR)
Model6<-rma.mv(yi,vi,mods=~NS_LF+IS_LF,random=list(~1|SiteID),method="ML",data=CC_RR)



AICc(Model0,Model1,Model2,Model3,Model4,Model5,Model6)

Model_pred<-data.frame(Pred=predict(Model4)$pred,ci.lb=predict(Model4)$ci.lb,ci.ub=predict(Model4)$ci.ub)
Model_pred$WI<-CC$IS_LF



Model_pred2<-Model_pred[!duplicated(Model_pred), ]
head(Model_pred2)

ggplot(Model_pred,aes(x=WI,y=exp(Pred)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1))+geom_pointrange(position=position_dodge(width=0.5))+geom_hline(y=0,lty=2)


#now water provision
WP<-subset(Inv,Serv_type=="Water provision")
WP<-WP[WP$WN %in% c("Yes","No"),]
WP$WN<-factor(WP$WN)
WP$NS_LF<-factor(WP$NS_LF)

WP_RR<-escalc(data=WP,measure="ROM",m2i=ES_UI,sd2i=SE_UI,n2i=SS_UI,m1i=ES_I,sd1i=SE_I,n1i=SS_I,append=T)
WP_RR<-subset(WP_RR,!is.na(yi))

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|SiteID),method="ML",data=WP_RR)
Model1<-rma.mv(yi,vi,mods=~WI,random=list(~1|SiteID),method="ML",data=WP_RR)
Model2<-rma.mv(yi,vi,mods=~WN,random=list(~1|SiteID),method="ML",data=WP_RR)
Model3<-rma.mv(yi,vi,mods=~WN+WI,random=list(~1|SiteID),method="ML",data=WP_RR)
Model4<-rma.mv(yi,vi,mods=~IS_LF,random=list(~1|SiteID),method="ML",data=WP_RR)
Model5<-rma.mv(yi,vi,mods=~NS_LF,random=list(~1|SiteID),method="ML",data=WP_RR)
Model6<-rma.mv(yi,vi,mods=~NS_LF+IS_LF,random=list(~1|SiteID),method="ML",data=WP_RR)



AICc(Model0,Model1,Model2,Model3,Model4,Model5,Model6)

Model_pred<-data.frame(Pred=predict(Model4)$pred,ci.lb=predict(Model4)$ci.lb,ci.ub=predict(Model4)$ci.ub)
Model_pred$WI<-WP$IS_LF



Model_pred2<-Model_pred[!duplicated(Model_pred), ]
head(Model_pred2)

ggplot(Model_pred,aes(x=WI,y=exp(Pred)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1))+geom_pointrange(position=position_dodge(width=0.5))+geom_hline(y=0,lty=2)

#now water quality
WQ<-subset(Inv,Serv_type=="Maintenance of water quality")
WQ<-WQ[WQ$WN %in% c("Yes","No"),]
WQ$WN<-factor(WQ$WN)
WQ$NS_LF<-factor(WQ$NS_LF)

WQ_RR<-escalc(data=WQ,measure="ROM",m2i=ES_UI,sd2i=SE_UI,n2i=SS_UI,m1i=ES_I,sd1i=SE_I,n1i=SS_I,append=T)
WQ_RR<-subset(WQ_RR,!is.na(yi))

Model0<-rma.mv(yi,vi,mods=~1,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model1<-rma.mv(yi,vi,mods=~WI,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model2<-rma.mv(yi,vi,mods=~WN,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model3<-rma.mv(yi,vi,mods=~WN+WI,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model4<-rma.mv(yi,vi,mods=~IS_LF,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model5<-rma.mv(yi,vi,mods=~NS_LF,random=list(~1|SiteID),method="ML",data=WQ_RR)
Model6<-rma.mv(yi,vi,mods=~NS_LF+IS_LF,random=list(~1|SiteID),method="ML",data=WQ_RR)



AICc(Model0,Model1,Model2,Model3,Model4,Model5,Model6)

Model5<-rma.mv(yi,vi,mods=~NS_LF,random=list(~1|SiteID),method="REML",data=WQ_RR)


Model_pred<-data.frame(Pred=predict(Model5)$pred,ci.lb=predict(Model5)$ci.lb,ci.ub=predict(Model5)$ci.ub)
Model_pred$WI<-WQ_RR$NS_LF



Model_pred2<-Model_pred[!duplicated(Model_pred), ]
head(Model_pred2)

ggplot(Model_pred,aes(x=WI,y=exp(Pred)-1,ymax=exp(ci.ub)-1,ymin=exp(ci.lb)-1))+geom_pointrange(position=position_dodge(width=0.5))+geom_hline(y=0,lty=2)

