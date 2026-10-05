#############################################################
# Clean-up the memory and start a new session
#############################################################
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
rm(list=ls())
dev.off()

library ("ggplot2")
library("dplyr")
library("MASS")
library("ggfortify")
library("ggthemes")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMRplus")
library("lsmeans")
library("readxl")
library('FSA')


#Set working directory

setwd("C:/Users/KA44139/Documents/PhD Year 1/SynCom")
setwd("C:/Users/catar/Documents/Dundee PhD Year 1/SynCom")
getwd()

#Import data
library(readxl)
JH201_Info <- read_excel("JH201_Info.xlsx")
View(JH201_Info)

#Make factors
JH201_Info$Genotype <- as.factor(JH201_Info$Genotype)
JH201_Info$Syncom_Status <- as.factor(JH201_Info$Syncom_Status)
JH201_Info$Opp_Rhizo <- as.factor(JH201_Info$Opp_Rhizo)
JH201_Info$SynCom_ID <- as.factor(JH201_Info$SynCom_ID)

#Order factors
JH201_Info$Genotype <- factor(JH201_Info$Genotype, levels = c("Morex", "Barke", "124_17", "124_52", "Bulk"))
JH201_Info$Syncom_Status <- factor(JH201_Info$Syncom_Status, levels = c("Live", "Dead"))
JH201_Info$SynCom_ID <- factor(JH201_Info$SynCom_ID, levels = c("RKMP", "PBS", "SDW", "None"))

####Run 5####
run5 <- subset(JH201_Info, Date_Inoc==230831)
run5 <- subset(run5, Colony_Count_Date==230918)
View(run5)

####Dry Weight####

#test distribution
shapiro.test(run5$`Dry_Weight_(mg)`)
hist(run5$`Dry_Weight_(mg)`, breaks=100)

#rkmp vs pbs
t.test(`Dry_Weight_(mg)`~SynCom_ID , data=run5)

#plot
dw5<-ggplot(run5, aes(x=SynCom_ID, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text = element_text(size=25))

dw5 + ggtitle("Dry Weight by Inoculant: RKMP vs PBS") +
  xlab("SynCom") + ylab("Weight (mg)") +
  scale_fill_manual(values=c("#238A8D", "#FFFFFF"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

