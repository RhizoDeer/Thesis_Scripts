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
JH34_Info <- read_excel("JH34_Info.xlsx", 
                        col_types = c("text", "numeric", "numeric", 
                                      "text", "text", "text", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric",
                                      "numeric","text", "text"))
View(JH34_Info)

#Make factors
JH34_Info$Genotype <- as.factor(JH34_Info$Genotype)
JH34_Info$Syncom_Status <- as.factor(JH34_Info$Syncom_Status)
JH34_Info$Opp_Rhizo <- as.factor(JH34_Info$Opp_Rhizo)
JH34_Info$SynCom_ID <- as.factor(JH34_Info$SynCom_ID)

#Order factors
JH34_Info$Genotype <- factor(JH34_Info$Genotype, levels = c("Morex", "Barke", "124_17", "124_52", "Bulk"))
JH34_Info$Syncom_Status <- factor(JH34_Info$Syncom_Status, levels = c("Live", "Dead"))

#############
####Run 4####
#############
run4 <- subset(JH34_Info, Date_Inoc==230525)
View(run4)

#Remove plant that got tangled in its seed coat and failed to grow properly
remove_fail4 <- subset(run4, Sample_ID!='KA193')
View(remove_fail4)

#Remove outlier plate
remove_fail4 <- subset(remove_fail4, Sample_ID!='KA197')
View(remove_fail4)

###########################
####Dry weight#############
###########################

#Subset plants only
plants_only4 <- subset(remove_fail4, Genotype!='Bulk')
View(plants_only4)

####By SynCom, and SynCom Status####
plants_only4$`Dry_Weight_(mg)`

shapiro.test(plants_only4$`Dry_Weight_(mg)`)
#must split up for further analysis

#Plot
r<-ggplot(plants_only4, aes(x=SynCom_ID, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by SynCom, SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight (mg)") + 
  scale_fill_manual(values=c("#FDE725", "#55C667", "#808080"))

####Subset by syncom####

#RKMP
rkmp_only <- subset(plants_only4, SynCom_ID=='RKMP')
View(rkmp_only)
shapiro.test(rkmp_only$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status , data=rkmp_only)

#Plot
r<-ggplot(rkmp_only, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by SynCom Status: RKMP") +
  xlab("SynCom Status") + ylab("Dry Weight (mg)") + 
  scale_fill_manual(values=c("#FDE725", "#55C667", "#808080"))

#RKMP_NoBi27
rkmpnb_only <- subset(plants_only4, SynCom_ID=='RKMP_NoBi27')
View(rkmpnb_only)
shapiro.test(rkmpnb_only$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status , data=rkmpnb_only)

#plot
r<-ggplot(rkmpnb_only, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by SynCom Status: RKMP_No27") +
  xlab("SynCom Status") + ylab("Dry Weight (mg)") + 
  scale_fill_manual(values=c("#FDE725", "#55C667", "#808080"))

####Subset by status + uninoculated####

#live and control
live_only <- subset(plants_only4, Syncom_Status!='Dead'| is.na(plants_only4$Syncom_Status))
View(live_only)
shapiro.test(live_only$`Dry_Weight_(mg)`)

kruskal.test(`Dry_Weight_(mg)`~SynCom_ID , data=live_only)
kwAllPairsDunnTest(live_only$`Dry_Weight_(mg)`~live_only$SynCom_ID, p.adjust.method="BH")

####THESIS####
#plot
r<-ggplot(live_only, aes(x=SynCom_ID, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,50)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c("#40B0A6", "#E1BE6A", "#FFFFFF"))+
  scale_x_discrete(labels = c("RKMP", "RKMP-27","SDW"))
r

#dead and control
dead_only <- subset(plants_only4, Syncom_Status!='Live'| is.na(plants_only4$Syncom_Status))
View(dead_only)
shapiro.test(dead_only$`Dry_Weight_(mg)`)

aov(dead_only$`Dry_Weight_(mg)`~dead_only$SynCom_ID)
summary(aov(dead_only$`Dry_Weight_(mg)`~dead_only$SynCom_ID))
TukeyHSD(aov(dead_only$`Dry_Weight_(mg)`~dead_only$SynCom_ID))

####THESIS####
#plot
r<-ggplot(dead_only, aes(x=SynCom_ID, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,50)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c("#40B0A6", "#E1BE6A", "#FFFFFF"))+
  scale_x_discrete(labels = c("RKMP", "RKMP-27","SDW"))
r

##################################
####RKMP_NoBi27 ~ Runs 2 and 4####
##################################
rkmpnbruns <- subset(JH34_Info, Date_Inoc!=230216)
rkmpnbruns <- subset(rkmpnbruns, Date_Inoc!=230406)
View(rkmpnbruns)

#Remove failed samples
rkmpnbruns <- subset(rkmpnbruns, Sample_ID!='KA193')
rkmpnbruns <- subset(rkmpnbruns, Sample_ID!='KA83')
rkmpnbruns <- subset(rkmpnbruns, Sample_ID!='KA197')
View(rkmpnbruns)

#barke and bulk only
rkmpnbruns <- subset(rkmpnbruns, Genotype!='124_17')
rkmpnbruns <- subset(rkmpnbruns, Genotype!='124_52')
rkmpnbruns <- subset(rkmpnbruns, Genotype!='Morex')
View(rkmpnbruns)

#rkmpnb only
rkmpnbruns <- subset(rkmpnbruns, SynCom_ID =='RKMP_NoBi27')
View(rkmpnbruns)

rkmpnbruns$Date_Inoc <- as.factor(rkmpnbruns$Date_Inoc)

#test distribution
shapiro.test(rkmpnbruns$CFU)
hist(rkmpnbruns$CFU, breaks=100)

#adjust to log scale and test again
rkmpnbruns$logCFU <- log10(rkmpnbruns$CFU+1)
View(rkmpnbruns)
hist(rkmpnbruns$logCFU, breaks=100)
shapiro.test(rkmpnbruns$logCFU)

##################
####Dry Weight####
##################

#planted only
plant_rkmpnbruns <- subset(rkmpnbruns, Genotype!='Bulk')
View(plant_rkmpnbruns)

####By run and status####
shapiro.test(plant_rkmpnbruns$`Dry_Weight_(mg)`)
hist(plant_rkmpnbruns$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Date_Inoc*Syncom_Status , data=plant_rkmpnbruns)
summary(aov(`Dry_Weight_(mg)`~Date_Inoc*Syncom_Status , data=plant_rkmpnbruns))
TukeyHSD(aov(`Dry_Weight_(mg)`~Date_Inoc*Syncom_Status , data=plant_rkmpnbruns))

#plot
r<-ggplot(plant_rkmpnbruns, aes(x=Date_Inoc, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Run and Status") +
  xlab("Inoculation Date") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))

##################################
####RKMP ~ Runs 3 and 4####
##################################
rkmpruns <- subset(JH34_Info, Date_Inoc!=230216)
rkmpruns <- subset(rkmpruns, Date_Inoc!=230223)
View(rkmpruns)

#Remove failed samples
rkmpruns <- subset(rkmpruns, Sample_ID!='KA193')
rkmpruns <- subset(rkmpruns, Sample_ID!='KA83')
rkmpruns <- subset(rkmpruns, Sample_ID!='KA197')
View(rkmpruns)

#barke and bulk only
rkmpruns <- subset(rkmpruns, Genotype!='124_17')
rkmpruns <- subset(rkmpruns, Genotype!='124_52')
rkmpruns <- subset(rkmpruns, Genotype!='Morex')
View(rkmpruns)

#rkmpnb only
rkmpruns <- subset(rkmpruns, SynCom_ID =='RKMP')
View(rkmpruns)

rkmpruns$Date_Inoc <- as.factor(rkmpruns$Date_Inoc)

#test distribution
shapiro.test(rkmpruns$CFU)
hist(rkmpruns$CFU, breaks=100)

#adjust to log scale and test again
rkmpruns$logCFU <- log10(rkmpruns$CFU+1)
View(rkmpruns)
hist(rkmpruns$logCFU, breaks=100)
shapiro.test(rkmpruns$logCFU)

##################
####Dry Weight####
##################

#planted only
plant_rkmpruns <- subset(rkmpruns, Genotype!='Bulk')
View(plant_rkmpruns)

####By run and status####
shapiro.test(plant_rkmpruns$`Dry_Weight_(mg)`)
hist(plant_rkmpruns$`Dry_Weight_(mg)`)

#plot
r<-ggplot(plant_rkmpruns, aes(x=Date_Inoc, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Run and Status") +
  xlab("Inoculation Date") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))

####live only####
live_plant_rkmpruns <- subset(plant_rkmpruns, Syncom_Status!='Dead')
View(live_plant_rkmpruns)

shapiro.test(live_plant_rkmpruns$`Dry_Weight_(mg)`)
t.test(live_plant_rkmpruns$`Dry_Weight_(mg)`~live_plant_rkmpruns$Date_Inoc)

#plot
r<-ggplot(live_plant_rkmpruns, aes(x=Date_Inoc, y=`Dry_Weight_(mg)`, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Run: Live") +
  xlab("Inoculation Date") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6"))

####dead only####
dead_plant_rkmpruns <- subset(plant_rkmpruns, Syncom_Status!='Live')
View(dead_plant_rkmpruns)

shapiro.test(dead_plant_rkmpruns$`Dry_Weight_(mg)`)
t.test(dead_plant_rkmpruns$`Dry_Weight_(mg)`~dead_plant_rkmpruns$Date_Inoc)

#plot
r<-ggplot(dead_plant_rkmpruns, aes(x=Date_Inoc, y=`Dry_Weight_(mg)`, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(0,50)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Run: Dead") +
  xlab("Inoculation Date") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#E1BE6A", "#40B0A6"))
