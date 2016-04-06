#this script is to analyse how the change in plant height affects aboveground biomass storage

rm(list=ls())

library(ggplot2)
library(lme4)
library(metafor)
library(MuMIn)
library(sp)
library(ncf)

#and functions
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

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


#now calculate effect sizes in metafor
AGB2<-subset(AGB,!is.na(Height_RR)&Terr_aqu=="Terrestrial")
AGB_ES<-escalc("ROM",m2i=EF_UI,m1i=EF_I,sd2i=SE_UI,sd1i=SE_I,n2i=SS_UI,n1i=SS_I,data=AGB2)
Site_unique<-unique(AGB_ES$SiteID)
head(AGB_ES)

#produce a correlogram to look at spatial autocorrelation of effect sizes

AGB.cor<-spline.correlog(AGB_ES$Latitude, AGB_ES$Longitude, AGB_ES$yi, latlon=T,resamp=1000,quiet=T)
AGB.cor2<-data.frame(Dist=AGB.cor$boot$boot.summary$predicted$x[1,],Cor=AGB.cor$boot$boot.summary$predicted$y[6,],UCI=AGB.cor$boot$boot.summary$predicted$y[2,],LCI=AGB.cor$boot$boot.summary$predicted$y[10,])
theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(AGB.cor2,aes(x=Dist,y=Cor,ymax=LCI,ymin=UCI))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))
ggsave("Figures/AGB_correl.pdf",height=6,width=8,units="in",dpi=300)


?spline.correlog


Model_AIC_summary<-NULL
for (i in 1:100){
  print(i)
  AGB_samp<-NULL
  for (j in 1:length(Site_unique)){#sample one site for each study so that no reference site is used more than once
    AGB_sub<-subset(AGB_ES,SiteID==Site_unique[j])
    AGB_sub<-AGB_sub[sample(nrow(AGB_sub), 1), ] 
    AGB_samp<-rbind(AGB_sub,AGB_samp)
  }
  Model0<-rma.mv(yi~1,vi,random=list(~1|Study),data=AGB_samp,method="ML")
  Model1<-rma.mv(yi~Height_RR,vi,random=list(~1|Study),data=AGB_samp,method="ML")
  Model2<-rma.mv(yi~CWD2,vi,random=list(~1|Study),data=AGB_samp,method="ML")
  Model3<-rma.mv(yi~Height_RR+CWD2,vi,random=list(~1|Study),data=AGB_samp,method="ML")
  Model4<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=AGB_samp,method="ML")
  Model_AIC<-data.frame(AICc=c(Model0$fit.stats$ML[5],Model1$fit.stats$ML[5],Model2$fit.stats$ML[5],#produce AICc values for the models
                               Model3$fit.stats$ML[5],Model4$fit.stats$ML[5]))
  Model_AIC$Vars<-c("Null","Height","CWD","Height+CWD", #details of model variables
                    "Height*CWD")
  Model_AIC$logLik<-c(Model0$fit.stats$ML[1],Model1$fit.stats$ML[1],Model2$fit.stats$ML[1],#put logLiklihood in the table
                      Model3$fit.stats$ML[1],Model4$fit.stats$ML[1])
  Null_dev<-deviance(Model0)
  Dev<-c(deviance(Model0),deviance(Model1),deviance(Model2),deviance(Model3),deviance(Model4))#calculate deviance of models
  Model_AIC$R2<-round(1-(Dev/Null_dev),2) #calculate pseudo-r squared using model deviance
  Model_AIC$R2<-ifelse(Model_AIC$R2<0,0,Model_AIC$R2)
  Model_AIC<-Model_AIC[order(Model_AIC$AICc),] #reorder from lowest to highest
  Model_AIC$delta<-Model_AIC$AICc-Model_AIC$AICc[1]#calculate AICc delta
  Model_AIC$rel_lik<-round(exp((Model_AIC$AICc[1]-Model_AIC$AICc)/2),2)#calculate the relative likelihood of model
  Model_AIC$weight<-round(Model_AIC$rel_lik/(sum(Model_AIC$rel_lik)),2)
  Model_AIC$Run<-i
  Model_AIC$Rank<-seq(1,5,1) #rank models from 1-8 in terms of parsimony
  Model_AIC_summary<-rbind(Model_AIC,Model_AIC_summary)
}

Model_AIC_summary$Rank1<-ifelse(Model_AIC_summary$Rank==1,1,0)
#summarise the boostrapping routine by giving median values for model statistics - log liklihood, AICc delta AICc, R squared
Model_sel_boot<-ddply(Model_AIC_summary,.(Vars),summarise,Modal_rank=Mode(Rank),Prop_rank=sum(Rank1)/100,log_liklihood=median(logLik),AICc_med=median(AICc),
                      delta_med=median(delta),R2_med=median(R2))


#now boostrap the top model to get parameter estimates
Site_unique<-unique(AGB_ES$SiteID)
Param_boot<-NULL
Resid_corr<-NULL
for (i in 1:1000){
  print(i)
  AGB_samp<-NULL
  for (j in 1:length(Site_unique)){#use same routine as previously to subsample dataset avoiding pseudo-replication
    AGB_sub<-subset(AGB_ES,SiteID==Site_unique[j])
    AGB_sub<-AGB_sub[sample(nrow(AGB_sub), 1), ]
    AGB_samp<-rbind(AGB_sub,AGB_samp)
  }
  Model4<-rma.mv(yi~Height_RR*CWD2,vi,random=list(~1|Study),data=AGB_samp,method="REML")
  Param_vals<-data.frame(Parameter=c("Intercept","Height","CWD","Height*CWD"),
                         estimate=round(coef(summary(Model4))[1],2),
                         se=round(coef(summary(Model4))[2],2),
                         pval=round(coef(summary(Model4))[4],3),
                         ci_lb=round(coef(summary(Model4))[5],2),
                         ci_ub=round(coef(summary(Model4))[6],2))
  Param_boot<-rbind(Param_vals,Param_boot)
  AGB.cor<-spline.correlog(AGB_samp$Latitude, AGB_samp$Longitude, resid(Model4), latlon=T,resamp=1000,quiet=T)
  AGB.cor2<-data.frame(Dist=AGB.cor$boot$boot.summary$predicted$x[1,],Cor=AGB.cor$boot$boot.summary$predicted$y[6,],UCI=AGB.cor$boot$boot.summary$predicted$y[2,],LCI=AGB.cor$boot$boot.summary$predicted$y[10,])
  Resid_corr<-rbind(Resid_corr,AGB.cor2)
}

#produce summary of parameter estimates
Param_boot_sum<-ddply(Param_boot,.(Parameter),summarise,coef_estimate=median(estimate),lower=median(ci.lb),
                      upper=median(ci.ub),med_pval=median(pval),se=median(se))



#write this table of parameter estimates
write.table(Param_boot_sum,file="Tables/AGB_parameter_estimates.csv",sep=",")

#produce summary of analysis of spatial auto-correlation
head(Resid_corr)
SA_boot_sum<-ddply(Resid_corr,.(Dist),summarise,Cor_mean=mean(Cor),lower=mean(LCI),
                      upper=mean(UCI))

theme_set(theme_bw(base_size=12))
Moran_plot1<-ggplot(SA_boot_sum,aes(x=Dist,y=Cor_mean,ymax=upper,ymin=lower))+geom_ribbon(alpha=0.2)+geom_line(size=0.5,colour="black")
Moran_plot2<-Moran_plot1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")
Moran_plot2+geom_hline(y=0,lty=2)+xlab("Distance between sites (km)")+ylab("Moran's I correlation")+scale_size_continuous(range = c(1,3))
ggsave("Figures/AGB_resid_correl.pdf",height=6,width=8,units="in",dpi=300)


##############################################
#produce figures for AGB######################
##############################################


#create new dataset for predictions
new.data<-expand.grid(Height_RR=seq(min(AGB_ES$Height_RR),max(AGB_ES$Height_RR),length.out = 1000),
                      CWD2=seq(min(AGB_ES$CWD2,na.rm = T),max(AGB_ES$CWD2,na.rm = T),length.out = 1000))
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

new.data$yi<-(new.data$Height_RR*0.88) + 0.45 + (0.*new.data$CWD2) + ((new.data$CWD2*new.data$Height_RR)*-0.16)

# Plot
theme_set(theme_bw(base_size=12))
P1<-ggplot(new.data, aes(x=Height_RR,y=(CWD2*sd(AGB$CWD))+mean(AGB$CWD),fill=yi)) +
  geom_raster() + 
  scale_fill_gradient2("Change in community biomass \n(log response ratio)",low="blue",mid="light grey",high="red")+geom_point(shape=1,size=2,data=AGB2,fill=NA)
P2<-P1+theme(panel.border = element_rect(size=1.5,colour="black",fill=NA))
P2+xlab("Difference between invasive and native species height\n(log response ratio)")+ylab("Climatic water deficit (mm)")
ggsave("Figures/AGB_climate.pdf",width = 8,height = 4,units = "in",dpi = 400)
