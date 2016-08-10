#this script is to analyse how the change in plant height affects soil carbon

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)
library(plyr)
library(ncf)
library(scales)
library(boot)

#and functions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

twosidep<-function(data){
  p1<-sum(data>0)/length(data)
  p2<-sum(data<0)/length(data)
  p<-min(p1,p2)*2
  return(p)
}


#read in data
Carb<-read.csv("Data/SC_trait_climate2.csv")
head(Carb)
Carb<-subset(Carb, select = -c(Native_C3,Inv_C3,Temp,Precip,Native_woodiness,Inv_woodiness) )
Carb<-subset(Carb, Terr_aqu=="Terrestrial")
Carb<-Carb[complete.cases(Carb),]

Carb$Height_diff<-Carb$Inv_height-Carb$Native_height
Carb$Height_RR<-log(Carb$Inv_height)-log(Carb$Native_height)
Carb$Height_Perc<-Carb$Inv_height/Carb$Native_height
Carb$Height_prop<-(Carb$Inv_height-Carb$Native_height)/Carb$Native_height
Carb$Diff_RR<-log(Carb$EF_I)-log(Carb$EF_UI)
Carb$CWD2<-(Carb$CWD-mean(Carb$CWD))/sd(Carb$CWD)
Carb$Height_RR2<-(Carb$Height_RR-mean(Carb$Height_RR))/sd(Carb$Height_RR)

#correct the variance for all measurements so that they are equal to SD
Carb$SE_UI<-ifelse(Carb$Var=="SE",Carb$SE_UI/sqrt(Carb$SS_UI),Carb$SE_UI)
Carb$SE_I<-ifelse(Carb$Var=="SE",Carb$SE_I/sqrt(Carb$SS_I),Carb$SE_I)

#now calculate effect sizes in metafor
Carb_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=Carb)



Ref_unique<-unique(Carb_ES$Study)
Model_AIC_summary<-NULL
for (i in 1:10000) {
  print(i)
  Carb_samp<-NULL
  for (j in 1:length(Ref_unique)){#sample one site for each study so that no reference site is used more than once
    Carb_sub<-subset(Carb_ES,Study==Ref_unique[j])
    Carb_sub<-Carb_sub[sample(nrow(Carb_sub), 1), ] 
    Carb_samp<-rbind(Carb_sub,Carb_samp)
  }
  Model0<-rma(yi~1,vi,data=Carb_samp,method="ML")
  Model1<-rma(yi~Inv_height,vi,data=Carb_samp,method="ML")
  Model2<-rma(yi~Inv_height+CWD2,vi,data=Carb_samp,method="ML")
  Model3<-rma(yi~Inv_height*CWD2,vi,data=Carb_samp,method="ML")
  Model4<-rma(yi~Height_RR,vi,data=Carb_samp,method="ML")
  Model5<-rma(yi~CWD2,vi,data=Carb_samp,method="ML")
  Model6<-rma(yi~Height_RR+CWD2,vi,data=Carb_samp,method="ML")
  Model7<-rma(yi~Height_RR*CWD2,vi,data=Carb_samp,method="ML")

  Model0$tau2
  
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5],Model4$fit.stats$ML[5]))
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5],Model4$fit.stats$ML[5],Model5$fit.stats$ML[5],Model6$fit.stats$ML[5],Model7$fit.stats$ML[5]))
  Model_AIC$Vars<-c("Null","Height","Height+CWD","Height*CWD","Relative height", "CWD", #details of model variables
                    "Relative height+CWD","Relative Height*CWD")
  Model_AIC$logLik<-c(Model0$fit.stats$ML[1],Model1$fit.stats$ML[1],Model2$fit.stats$ML[1],#put logLiklihood in the table
                      Model3$fit.stats$ML[1],Model4$fit.stats$ML[1],Model5$fit.stats$ML[1],
                      Model6$fit.stats$ML[1],Model7$fit.stats$ML[1])
  Null_dev<-Model0$tau2
  Dev<-c(Model0$tau2,Model1$tau2,Model2$tau2,Model3$tau2,Model4$tau2,Model5$tau2,
         Model6$tau2,Model7$tau2)#calculate deviance of models
  Model_AIC$R2<-round(1-(Dev/Null_dev),2) #calculate pseudo-r squared using model deviance
  Model_AIC$R2<-ifelse(Model_AIC$R2<0,0,Model_AIC$R2)
  Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
  Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
  Model_AIC$rel_lik<-round(exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2),2)#calculate the relative likelihood of model
  Model_AIC$weight<-round(Model_AIC$rel_lik/(sum(Model_AIC$rel_lik)),2)
  Model_AIC$Run<-i
  Model_AIC$Rank<-seq(1,8,1) #rank models from 1-8 in terms of parsimony
  Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}
Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)
#summarise the boostrapping routine by giving median values for model statistics - log liklihood, AICc delta AICc, R squared
Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/10000,log_liklihood=median(logLik),AICc_med=median(AICc),
                      delta_med=median(delta),R2_med=median(R2))
write.csv(Model_sel_boot,"Tables/Soil_model_sel.csv")


#now boostrap the top model to get parameter estimates
Ref_unique<-unique(Carb_ES$Study)
Param_boot<-NULL
Resid_corr<-NULL
for (i in 1:10000){
  print(i)
  Carb_samp<-NULL
  for (j in 1:length(Ref_unique)){#use same routine as previously to subsample dataset avoiding pseudo-replication
    Carb_sub<-subset(Carb_ES,Study==Ref_unique[j])
    Carb_sub<-Carb_sub[sample(nrow(Carb_sub), 1), ]
    Carb_samp<-rbind(Carb_sub,Carb_samp)
  }
  Model4<-rma(yi~Height_RR*CWD,vi,data=Carb_samp,method="REML")
  Param_vals<-data.frame(Parameter=c("Intercept","Height","CWD","Height*CWD"),
                         estimate=round(coef(summary(Model4))[1],5),
                         se=round(coef(summary(Model4))[2],5),
                         pval=round(coef(summary(Model4))[4],5),
                         ci_lb=round(coef(summary(Model4))[5],5),
                         ci_ub=round(coef(summary(Model4))[6],5))
  Param_boot<-rbind(Param_vals,Param_boot)
  #Soil.cor<-spline.correlog(Carb_samp$Latitude, Carb_samp$Longitude, resid(Model4), latlon=T,resamp=1000,quiet=T)
  #Soil.corr2<-data.frame(Dist=Soil.cor$boot$boot.summary$predicted$x[1,],Cor=Soil.cor$boot$boot.summary$predicted$y[6,],UCI=Soil.cor$boot$boot.summary$predicted$y[2,],LCI=Soil.cor$boot$boot.summary$predicted$y[10,])
  #Resid_corr<-rbind(Resid_corr,Soil.corr2)
}


#produce summary of parameter estimates
Param_boot_sum<-ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
                      upper=median(ci.ub),med_pval=median(pval),se=median(se))

#write this table of parameter estimates
write.table(Param_boot_sum,file="Tables/Soil_C_parameter_estimates.csv",sep=",",row.names=F)


#create new dataset for predictions
new.data<-data.frame(Height_RR=seq(min(Carb_ES$Height_RR),max(Carb_ES$Height_RR),length.out = 1000))
# Calculate yi as you did:
new.data$yi<-(new.data$Height_RR*0.09) + 0.16
new.data$yi_lo<-(new.data$Height_RR*0) + 0.01
new.data$yi_hi<-(new.data$Height_RR*0.18) + 0.31
# Plot

theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=yi))+geom_line(size=1.5)+
  geom_point(data=Carb_ES,shape=1,size=2)
P1+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Change in community biomass \n(log response ratio)")
ggsave("Figures/Carb_height.pdf",width = 8,height = 4,units = "in",dpi = 400)
ggsave("Figures/Carb_height.png",width = 8,height = 4,units = "in",dpi = 400)
