#script to sort out trait data and check for inconsistencies

library(ggplot2)

#load in data
Species_traits<-read.csv("Data/Species_resolved_traits.csv")


Species_traits$Median_height<-ifelse(!is.na(Species_traits$Single_height),Species_traits$Single_height,(Species_traits$Upper_height+Species_traits$Lower_height)/2)

write.csv(Species_traits,"Data/Species_traits_median.csv")

head(Species_traits)


ggplot(Species_traits,aes(y= Trait_db_height,x=Upper_height))+geom_point()+geom_abline()+geom_smooth(method="lm")
ggplot(Species_traits,aes(y= Trait_db_height,x=Lower_height))+geom_point()+geom_abline()+geom_smooth(method="lm")
ggplot(Species_traits,aes(y= Trait_db_height,x=Median_height))+geom_point()+geom_abline()+geom_smooth(method="lm")



#look at how close height values derived from internet are to trait databases
#first calculated as a percentage
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Upper_height)/Species_traits$Trait_db_height)^2,na.rm = T))-1
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Lower_height)/Species_traits$Trait_db_height)^2,na.rm = T))-1
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Median_height)/Species_traits$Trait_db_height)^2,na.rm = T))-1

#then in metres
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Upper_height))^2,na.rm = T))
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Lower_height))^2,na.rm = T))
sqrt(mean(((Species_traits$Trait_db_height-Species_traits$Median_height))^2,na.rm = T))


#plot how this varies by height

plot(Species_traits$Trait_db_height,Species_traits$Median_height-Species_traits$Trait_db_height)


M1<-lm(Trait_db_height~Median_height,data=Species_traits)
summary(M1)
