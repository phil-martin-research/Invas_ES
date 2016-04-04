#this script is to analyse how the change in plant height affects soil moisture

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)
library(plyr)
#and functions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#read in data
SW<-read.csv("Data/Moisture_trait_climate.csv")
head(SW)
SW<-subset(SW, select = -c(Native_C3,Inv_C3) )
SW<-subset(SW, Terr_aqu=="Terrestrial" )
SW<-SW[complete.cases(SW),]

SW$Height_diff<-SW$Inv_height-SW$Native_height
SW$Height_RR<-log(SW$Inv_height)-log(SW$Native_height)
SW$Height_prop<-(SW$Inv_height-SW$Native_height)/SW$Native_height
SW$Woody_diff<-SW$Inv_woodiness-SW$Native_woodiness
SW$Diff_RR<-log(SW$EF_I)-log(SW$EF_UI)
SW$CWD2<-(SW$CWD-mean(SW$CWD))/sd(SW$CWD)

#correct the variance for all measurements so that they are equal to SD
SW$SE_UI<-ifelse(SW$Var=="SE",SW$SE_UI/sqrt(SW$SS_UI),SW$SE_UI)
SW$SE_I<-ifelse(SW$Var=="SE",SW$SE_I/sqrt(SW$SS_I),SW$SE_I)


#do some data exploration
ggplot(SW,aes(x=Height_RR,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)
ggplot(SW,aes(x=Woody_diff,y=Diff_RR))+geom_point()+geom_smooth(method="lm",se=F)


#now calculate effect sizes in metafor
SW_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=SW)
Site_unique<-unique(SW_ES$SiteID)
Model_AIC_summary<-NULL
for (i in 1:100){
  print(i)
  SW_samp<-NULL
  for (j in 1:length(Site_unique)){#sample one site for each study so that no reference site is used more than once
    SW_sub<-subset(SW_ES,SiteID==Site_unique[j])
    SW_sub<-SW_sub[sample(nrow(SW_sub), 1), ] 
    SW_samp<-rbind(SW_sub,SW_samp)
  }
  Model0<-rma.mv(yi~1,vi,random=list(~1|Study),data=SW_samp,method="ML")
  Model1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=SW_samp,method="ML")
  Model2<-rma.mv(yi~CWD2,vi,random=list(~1|Study),data=SW_samp,method="ML")
  Model3<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=SW_samp,method="ML")
  summary(Model0)
  summary(Model1)
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5]))
  Model_AIC$Vars<-c("Null","Height_RR","CWD","Height*CWD") #details of model variables
  Model_AIC$logLik<-c(Model0$fit.stats$ML[1],Model1$fit.stats$ML[1],Model2$fit.stats$ML[1],#put logLiklihood in the table
                      Model3$fit.stats$ML[1])
  Null_dev<-deviance(Model0)
  Dev<-c(deviance(Model0),deviance(Model1),deviance(Model2),deviance(Model3))
  Model_AIC$R2<-1-(Dev/Null_dev) #calculate pseudo-r squared using model deviance
  Model_AIC$R2<-ifelse(Model_AIC$R2<0,0,Model_AIC$R2)
  Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
  Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
  Model_AIC$rel_lik<-exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2)#calculate the relative likelihood of model
  Model_AIC$weight<-Model_AIC$rel_lik/(sum(Model_AIC$rel_lik))
  Model_AIC$Run<-i
  Model_AIC$Rank<-seq(1,4,1) #rank models from 1-4 in terms of parsimony
  Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}
Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)
#summarise the boostrapping routine by giving median values for model statistics - log liklihood, AICc delta AICc, R squared
Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/100,log_liklihood=median(logLik),AICc_med=median(AICc),
                      delta_med=median(delta),R2_med=median(R2))


#now boostrap the top model to get parameter estimates
Site_unique<-unique(SW_ES$SiteID)
Param_boot<-NULL
for (i in 1:100){
  print(i)
  SW_samp<-NULL
  for (j in 1:length(Site_unique)){#use same routine as previously to subsample dataset avoiding pseudo-replication
    SW_sub<-subset(SW_ES,SiteID==Site_unique[j])
    SW_sub<-SW_sub[sample(nrow(SW_sub), 1), ]
    SW_samp<-rbind(SW_sub,SW_samp)
  }
  Model1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=SW_samp,method="ML")
  Param_vals<-data.frame(Parameter=c("Intercept","Height"),
                         estimate=round(coef(summary(Model1))[1],2),
                         se=round(coef(summary(Model1))[2],2),
                         pval=round(coef(summary(Model1))[4],3),
                         ci_lb=round(coef(summary(Model1))[5],2),
                         ci_ub=round(coef(summary(Model1))[6],2))
  Param_boot<-rbind(Param_vals,Param_boot)
}

#produce summary of parameter estimates
Param_boot_sum<-ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
                      upper=median(ci.ub),med_pval=median(pval),se=median(se))

#write this table of parameter estimates
write.table(Param_boot_sum,file="Tables/Soil_moist_parameter_estimates.csv",sep=",")

#create new dataset for predictions
new.data<-expand.grid(Height_RR=seq(min(SW_ES$Height_RR),max(SW_ES$Height_RR),length.out = 1000))
#create new dataframe to produce convex hull
new.data$yi<-(new.data$Height_RR*0.02)-0.14

theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=yi)) +geom_point()+geom_point(shape=1,size=2,data=SW_ES,aes(x=Height_RR,y=yi))
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Change in soil moisture (log ratio)")
ggsave("Figures/SW_climate.pdf",width = 6,height = 4,units = "in",dpi = 400)
