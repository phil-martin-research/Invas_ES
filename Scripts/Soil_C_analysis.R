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

#and functions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#read in data
Carb<-read.csv("Data/SC_trait_climate2.csv")
head(Carb)
Carb<-subset(Carb, select = -c(Native_C3,Inv_C3) )
Carb<-subset(Carb, Terr_aqu=="Terrestrial" )
Carb<-Carb[complete.cases(Carb),]

Carb$Height_diff<-Carb$Inv_height-Carb$Native_height
Carb$Height_RR<-log(Carb$Inv_height)-log(Carb$Native_height)
Carb$Height_Perc<-Carb$Inv_height/Carb$Native_height
Carb$Height_prop<-(Carb$Inv_height-Carb$Native_height)/Carb$Native_height
Carb$Woody_diff<-Carb$Inv_woodiness-Carb$Native_woodiness
Carb$Diff_RR<-log(Carb$EF_I)-log(Carb$EF_UI)
Carb$CWD2<-(Carb$CWD-mean(Carb$CWD))/sd(Carb$CWD)

#correct the variance for all measurements so that they are equal to SD
Carb$SE_UI<-ifelse(Carb$Var=="SE",Carb$SE_UI/sqrt(Carb$SS_UI),Carb$SE_UI)
Carb$SE_I<-ifelse(Carb$Var=="SE",Carb$SE_I/sqrt(Carb$SS_I),Carb$SE_I)

#now calculate effect sizes in metafor
Carb_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=Carb)

Ref_unique<-unique(Carb_ES$EF_UI)
Model_AIC_summary<-NULL
for (i in 1:1000) {
  print(i)
  Carb_samp<-NULL
  for (j in 1:length(Ref_unique)){#sample one site for each study so that no reference site is used more than once
    Carb_sub<-subset(Carb_ES,EF_UI==Ref_unique[j])
    Carb_sub<-Carb_sub[sample(nrow(Carb_sub), 1), ] 
    Carb_samp<-rbind(Carb_sub,Carb_samp)
  }
  Model0<-rma.mv(yi~1,vi,random=list(~1|Study),data=Carb_samp,method="ML")
  Model1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=Carb_samp,method="ML")
  Model2<-rma.mv(yi~CWD2,vi,random=list(~1|Study),data=Carb_samp,method="ML")
  Model3<-rma.mv(yi~Height_Perc+CWD2,vi,random=list(~1|Study),data=Carb_samp,method="ML")
  Model4<-rma.mv(yi~Height_Perc*CWD2,vi,random=list(~1|Study),data=Carb_samp,method="ML")
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5],Model4$fit.stats$ML[5]))
  Model_AIC$Vars<-c("Null","Height_RR","CWD","Height+CWD","Height*CWD") #details of model variables
  Model_AIC$logLik<-c(Model0$fit.stats$ML[1],Model1$fit.stats$ML[1],Model2$fit.stats$ML[1],#put logLiklihood in the table
                      Model3$fit.stats$ML[1],Model4$fit.stats$ML[1])
  ggplot(Carb_samp,aes(x=SS_I+SS_UI,y=((1/vi)/528467.4)*100,label=SiteID))+geom_point()+geom_text(size=2)+scale_y_continuous(labels = comma)
  sum(1/Carb_samp$vi)
  
  
  plot(fitted(Model4),resid(Model4))
  Null_dev<-deviance(Model0)
  Dev<-c(deviance(Model0),deviance(Model1),deviance(Model2),deviance(Model3),deviance(Model4))
  Model_AIC$R2<-1-(Dev/Null_dev) #calculate pseudo-r squared using model deviance
  Model_AIC$R2<-ifelse(Model_AIC$R2<0,0,Model_AIC$R2)
  Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
  Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
  Model_AIC$rel_lik<-exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2)#calculate the relative likelihood of model
  Model_AIC$weight<-Model_AIC$rel_lik/(sum(Model_AIC$rel_lik))
  Model_AIC$Run<-i
  Model_AIC$Rank<-seq(1,5,1) #rank models from 1-8 in terms of parsimony
  Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}
Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)
#summarise the boostrapping routine by giving median values for model statistics - log liklihood, AICc delta AICc, R squared
Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/1000,log_liklihood=median(logLik),AICc_med=median(AICc),
                      delta_med=median(delta),R2_med=median(R2))
write.csv(Model_sel_boot,"Soil_model_sel.csv")


#now boostrap the top model to get parameter estimates
Ref_unique<-unique(Carb_ES$EF_UI)
Param_boot<-NULL
Resid_corr<-NULL
for (i in 1:1000){
  print(i)
  Carb_samp<-NULL
  for (j in 1:length(Ref_unique)){#use same routine as previously to subsample dataset avoiding pseudo-replication
    Carb_sub<-subset(Carb_ES,EF_UI==Ref_unique[j])
    Carb_sub<-Carb_sub[sample(nrow(Carb_sub), 1), ]
    Carb_samp<-rbind(Carb_sub,Carb_samp)
  }
  for (j in seq(1,3)) {Carb_samp<-subset(Carb_samp,vi>min(Carb_samp$vi))}
  Model4<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=Carb_samp,method="REML")
  Param_vals<-data.frame(Parameter=c("Intercept","Height","CWD","Height*CWD"),
                         estimate=round(coef(summary(Model4))[1],2),
                         se=round(coef(summary(Model4))[2],2),
                         pval=round(coef(summary(Model4))[4],3),
                         ci_lb=round(coef(summary(Model4))[5],2),
                         ci_ub=round(coef(summary(Model4))[6],2))
  Param_boot<-rbind(Param_vals,Param_boot)
  Soil.cor<-spline.correlog(Carb_samp$Latitude, Carb_samp$Longitude, resid(Model4), latlon=T,resamp=1000,quiet=T)
  Soil.corr2<-data.frame(Dist=Soil.cor$boot$boot.summary$predicted$x[1,],Cor=Soil.cor$boot$boot.summary$predicted$y[6,],UCI=Soil.cor$boot$boot.summary$predicted$y[2,],LCI=Soil.cor$boot$boot.summary$predicted$y[10,])
  Resid_corr<-rbind(Resid_corr,Soil.corr2)
}

#produce summary of parameter estimates
Param_boot_sum<-ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
                      upper=median(ci.ub),med_pval=median(pval),se=median(se))

#write this table of parameter estimates
write.table(Param_boot_sum,file="Tables/Soil_C_parameter_estimates.csv",sep=",")

#produce summary of analysis of spatial auto-correlation
head(Resid_corr)
SA_boot_sum<-ddply(Resid_corr,.(Dist),summarise,Cor_mean=mean(Cor),lower=mean(LCI),
                   upper=mean(UCI))

theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(SA_boot_sum,aes(x=Dist,y=Cor_mean,ymax=upper,ymin=lower))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))
ggsave("Figures/Soil_C_resid_correl.pdf",height=6,width=8,units="in",dpi=300)
ggsave("Figures/Soil_C_resid_correl.png",height=6,width=8,units="in",dpi=300)


ggplot(data=Carb_ES,aes(x=CWD2,y=Height_RR))+geom_point()


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
new.data$yi<-(new.data$Height_RR*0.10) + 0.19 + (1.29*new.data$CWD2) + ((new.data$CWD2*new.data$Height_RR)*0.03)

# Plot

theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=(CWD2*374.6319)+-551.739,fill=yi)) +
  geom_raster() + 
  scale_fill_gradient2("Change in soil carbon storage \n(log response ratio)",low="blue",mid="light grey",high="red")+geom_point(shape=1,size=2,data=Carb,fill=NA)
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Climatic water deficit (mm)")
ggsave("Figures/Carb_climate.pdf",width = 8,height = 4,units = "in",dpi = 400)
ggsave("Figures/Carb_climate.png",width = 8,height = 4,units = "in",dpi = 400)
