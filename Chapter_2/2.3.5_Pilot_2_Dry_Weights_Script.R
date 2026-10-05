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


#Set working directory

setwd("C:/Users/catar/Documents/Dundee PhD/Plants")
getwd()

#Import data 
library(readxl)
weight_data <- read_excel("Pilot_2_Dry_Weights.xlsx")
View(weight_data)
attach(weight_data)

#Test for normal distribution
shapiro.test(`Dry_Weight_(mg)`)

weight_avo <- aov(`Dry_Weight_(mg)`~Genotype)
summary(weight_avo)

weight_tukey <- TukeyHSD(weight_avo)
weight_tukey

#Plot
weight_plot <- ggplot(weight_data, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,45)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))+
  scale_fill_manual(values=c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7"))+
  scale_x_discrete(labels = c("Barke", "124-17", "124-52", "Morex"))+
  xlab("")+
  ylab("Dry Weight (mg)")+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", axis.title.x = element_blank())
weight_plot
