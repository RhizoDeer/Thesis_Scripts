#############################################################
# 
#  Code to compare chryseobacterium processing
#  Revision 07/26
#  c.arnton@dundee.ac.uk
#
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages 
library("vegan")
library ("ggplot2")
library("viridis")
library("gridExtra")

#############################################################
#set working directory
setwd("C:/Users/catar/Documents/Dundee PhD/SynCom")
getwd()

#import data
library(readxl)
ChryData <- read_excel("Chryseobacterium_Calibration_Data.xlsx")
View(ChryData)

#check structure
names(ChryData)
str(ChryData)

#Set factors
ChryData$Method <- as.factor(ChryData$Method)

#Order factors
ChryData$Method <- factor(ChryData$Method, levels = c("Liquid", "Lawn"))

#test distribution
shapiro.test(ChryData$`CFU/mL_at_OD600_1`)
#p>0.05 so parametric
hist((ChryData$`CFU/mL_at_OD600_1`), breaks=3)

####Orig. numbers####

#t-test
t.test(`CFU/mL_at_OD600_1`~Method, data=ChryData)

#Standard deviations - orig numbers
#Liquid
Liquid <- subset(ChryData, Method=="Liquid")
View(Liquid)
sd(Liquid$`CFU/mL_at_OD600_1`)
#Lawn
Lawn <- subset(ChryData, Method=="Lawn")
View(Lawn)
sd(Lawn$`CFU/mL_at_OD600_1`)

#Plot
p<-ggplot(ChryData, aes(x=Method, y=`CFU/mL_at_OD600_1`, fill=Method))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  scale_y_log10() +
  theme(legend.position="none")+
  scale_fill_manual(values = c("#44AA99", "#88CCEE"))+
  ylab("CFU/mL at OD600 = 1")+
  annotation_logticks(sides = "l")+
  coord_cartesian(ylim = c(1000000, 10000000000))
p

####log10####

#t-test
t.test(log10~Method, data=ChryData)

#Standard deviations
#Liquid
Liquid <- subset(ChryData, Method=="Liquid")
View(Liquid)
sd(Liquid$log10)
#Lawn
Lawn <- subset(ChryData, Method=="Lawn")
View(Lawn)
sd(Lawn$log10)

#Plot
p<-ggplot(ChryData, aes(x=Method, y=log10, fill=Method))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(6,10)+
  theme(legend.position="none")+
  scale_fill_manual(values = c("#44AA99", "#88CCEE"))+
  ylab("Log10(CFU/mL at OD600 = 1)")
p
