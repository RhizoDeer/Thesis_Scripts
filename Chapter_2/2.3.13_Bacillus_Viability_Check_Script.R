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
ViaData <- read_excel("C:/Users/KA44139/Downloads/Viability_Check_Info.xlsx", 
                                   range = "A1:I43")
View(ViaData)

#check structure
names(ViaData)
str(ViaData)

#Make factors factors and not factors not factors
ViaData$Isolate <- as.factor(ViaData$Isolate)
ViaData$Round <- as.factor(ViaData$Round)
ViaData$Date_Plated <- as.factor(ViaData$Date_Plated)
str(ViaData)

#Order factors
ViaData$Isolate <- factor(ViaData$Isolate, levels = c("Bi66", "Bi84", "BAR06"))

#Plot
ggplot(ViaData, aes(x=Round, y=`CFU/mL`, shape=Date_Plated, color=Isolate, group=interaction(Isolate,Date_Plated))) +
  geom_point()+
  geom_line()

#Convert to log scale
ViaData$logCFU <- log10(ViaData$`CFU/mL`)
View(ViaData)

#Plot on log scale
ggplot(ViaData, aes(x=Round, y=`logCFU`, shape=Date_Plated, color=Isolate, group=interaction(Isolate,Date_Plated))) +
  geom_point(size=2)+
  geom_line()+
  ylab('Log10(CFU/mL)')+
  ylim(4,11)

#Subset for just the second round
r2 <- subset(ViaData, ViaData$Date_Plated=="240410")
View(r2)

#Plot on log scale
ggplot(r2, aes(x=Round, y=`logCFU`, color=Isolate, group=Isolate)) +
  geom_point(size=2)+
  geom_line()+
  ylab('Log10(CFU/mL)')+
  ylim(4,11)+
  scale_color_manual(values=c("#785EF0", "#DC267F", "#FE6100"))+
  theme(legend.position = c(0.9, 0.85))
