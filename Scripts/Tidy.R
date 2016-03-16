#script to tidy up data on the effects of invasive species on ecosystem functions

ESI<-read.csv("Data/ES_invasives2.csv")
ESI2<-unique(ESI)

write.csv(ESI2,"Data/ES_invasives3.csv",row.names = F)
