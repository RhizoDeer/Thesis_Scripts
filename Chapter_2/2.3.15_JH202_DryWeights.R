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
library("gridExtra")


#Set working directory

#setwd("~/PhD Files/SynCom")
setwd("C:/Users/catar/Documents/Dundee PhD/SynCom/JH202")
#setwd("~/PhD Files/SynCom")
getwd()


#Import data
library(readxl)
#JH202_Info <- read_excel("C:/Users/KA44139/Documents/PhD Files/SynCom/JH202_Info.xlsx")
JH202_Info <- read_excel("C:/Users/catar/Documents/Dundee PhD/SynCom/JH202/JH202_Info.xlsx")
View(JH202_Info)

#check structure
names(JH202_Info)
str(JH202_Info)

#Set factors
JH202_Info$Genotype <- as.factor(JH202_Info$Genotype)
JH202_Info$Syncom_Status <- as.factor(JH202_Info$Syncom_Status)
JH202_Info$Date_Inoc <- as.factor(JH202_Info$Date_Inoc)
JH202_Info$Date_Harvest <- as.factor(JH202_Info$Date_Harvest)
JH202_Info$SynCom_ID <- as.factor(JH202_Info$SynCom_ID)


#Order factors
JH202_Info$Genotype <- factor(JH202_Info$Genotype, levels = c("Bulk", "Barke", "124-17", "124-52", "Morex"))
JH202_Info$Syncom_Status <- factor(JH202_Info$Syncom_Status, levels = c("Live", "HK", "Neither"))
JH202_Info$SynCom_ID <- factor(JH202_Info$SynCom_ID, levels = c("RKMP_Old", "RKMP_New", "KAS1", "PBS"))


#Remove unplanted
planted <- subset(JH202_Info, Genotype!="Bulk")
View(planted)

##################
#####All Runs#####
##################

#Root shoot ratio

#test distribution
shapiro.test(planted$Root_Shoot_Ratio)
#p<0.05 so non-parametric
hist((planted$Root_Shoot_Ratio), breaks=10)

#By syncom status
kruskal.test(Root_Shoot_Ratio ~ Syncom_Status, data = planted)
#p=0.02015
kwAllPairsDunnTest(planted$Root_Shoot_Ratio~planted$Syncom_Status, p.adjust.method="BH")

#Plot
p<-ggplot(planted, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#By syncom ID
kruskal.test(Root_Shoot_Ratio ~ SynCom_ID, data = planted)
#p=0.4693
kwAllPairsDunnTest(planted$Root_Shoot_Ratio~planted$SynCom_ID, p.adjust.method="BH")

#Plot
p<-ggplot(planted, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p


#By Genotype
kruskal.test(Root_Shoot_Ratio ~ Genotype, data = planted)
#p=6.801e-8
kwAllPairsDunnTest(planted$Root_Shoot_Ratio~planted$Genotype, p.adjust.method="BH")

#Plot
p<-ggplot(planted, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p


#By Run
kruskal.test(Root_Shoot_Ratio ~ Date_Inoc, data = planted)
#p=2.378e-05
kwAllPairsDunnTest(planted$Root_Shoot_Ratio~planted$Date_Inoc, p.adjust.method="BH")

#Plot
p<-ggplot(planted, aes(x=Date_Inoc, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

###############
#####Run 1#####
###############

#Subset for run 1
run1 <- subset(planted, Date_Inoc=="241003")
View(run1)

####Root shoot ratio####

#test distribution
shapiro.test(run1$Root_Shoot_Ratio)
#p>0.05 so parametric
hist((run1$Root_Shoot_Ratio), breaks=10)

#By Status and ID
aov(Root_Shoot_Ratio~Syncom_Status*SynCom_ID, data=run1)
summary(aov(Root_Shoot_Ratio~Syncom_Status*SynCom_ID, data=run1))

TukeyHSD(aov(Root_Shoot_Ratio~Syncom_Status*SynCom_ID, data=run1))

#Plot
p<-ggplot(run1, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by SynCom ID

#Old
Old1 <- subset(run1, SynCom_ID=="RKMP_Old")
View(Old1)
#test distribution
shapiro.test(Old1$Root_Shoot_Ratio)
hist((Old1$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~Syncom_Status, data = Old1)
a<-ggplot(Old1, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Old RKMP")
a

#New
New1 <- subset(run1, SynCom_ID=="RKMP_New")
View(New1)
#test distribution
shapiro.test(New1$Root_Shoot_Ratio)
hist((New1$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~Syncom_Status, data = New1)
b<-ggplot(New1, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("New RKMP")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by SynCom Status

#Live
Live1 <- subset(run1, Syncom_Status=="Live")
View(Live1)
#test distribution
shapiro.test(Live1$Root_Shoot_Ratio)
hist((Live1$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~SynCom_ID, data = Live1)
a<-ggplot(Live1, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Live")
a

#HK
HK1 <- subset(run1, Syncom_Status=="HK")
View(HK1)
#test distribution
shapiro.test(HK1$Root_Shoot_Ratio)
hist((HK1$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~SynCom_ID, data = HK1)
b<-ggplot(HK1, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("HK")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

####Above ground####

#test distribution
shapiro.test(run1$Dry_Weight_Above_mg)
hist((run1$Dry_Weight_Above_mg), breaks=10)

#By Status and ID
aov(Dry_Weight_Above_mg~Syncom_Status*SynCom_ID, data=run1)
summary(aov(Dry_Weight_Above_mg~Syncom_Status*SynCom_ID, data=run1))

TukeyHSD(aov(Dry_Weight_Above_mg~Syncom_Status*SynCom_ID, data=run1))

#Plot
p<-ggplot(run1, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

####PAPER####
#Subset by SynCom ID

#Old
Old1 <- subset(run1, SynCom_ID=="RKMP_Old")
View(Old1)
#test distribution
shapiro.test(Old1$Dry_Weight_Above_mg)
hist((Old1$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~Syncom_Status, data = Old1)
fillcol <- c("#0072B2", "#FFFFFF")
colcol <- c("#000000", "#0072B2")
a<-ggplot(Old1, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,45)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
a

#New
New1 <- subset(run1, SynCom_ID=="RKMP_New")
View(New1)
#test distribution
shapiro.test(New1$Dry_Weight_Above_mg)
hist((New1$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~Syncom_Status, data = New1)
b<-ggplot(New1, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,45)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by SynCom Status

#Live
Live1 <- subset(run1, Syncom_Status=="Live")
View(Live1)
#test distribution
shapiro.test(Live1$Dry_Weight_Above_mg)
hist((Live1$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~SynCom_ID, data = Live1)
c<-ggplot(Live1, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,45)+
  theme(legend.position="none")+
  ggtitle("Live")
c

#HK
HK1 <- subset(run1, Syncom_Status=="HK")
View(HK1)
#test distribution
shapiro.test(HK1$Dry_Weight_Above_mg)
hist((HK1$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~SynCom_ID, data = HK1)
d<-ggplot(HK1, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,45)+
  theme(legend.position="none")+
  ggtitle("HK")
d

#faceted graph of both
grid.arrange(c, d, ncol=2)

####Below ground####

#test distribution
shapiro.test(run1$Dry_Weight_Below_mg)
#p<0.05 so non-parametric
hist((run1$Dry_Weight_Below_mg), breaks=10)

#Plot
p<-ggplot(run1, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

####PAPER####
#Subset by SynCom ID

#Old
Old1 <- subset(run1, SynCom_ID=="RKMP_Old")
View(Old1)
#test distribution
shapiro.test(Old1$Dry_Weight_Below_mg)
hist((Old1$Dry_Weight_Below_mg), breaks=5)
wilcox.test(Dry_Weight_Below_mg~Syncom_Status, data = Old1)
a<-ggplot(Old1, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,110)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
a

#New
New1 <- subset(run1, SynCom_ID=="RKMP_New")
View(New1)
#test distribution
shapiro.test(New1$Dry_Weight_Below_mg)
hist((New1$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~Syncom_Status, data = New1)
b<-ggplot(New1, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,110)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by SynCom Status

#Live
Live1 <- subset(run1, Syncom_Status=="Live")
View(Live1)
#test distribution
shapiro.test(Live1$Dry_Weight_Below_mg)
hist((Live1$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = Live1)
a<-ggplot(Live1, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,110)+
  theme(legend.position="none")+
  ggtitle("Live")
a

#HK
HK1 <- subset(run1, Syncom_Status=="HK")
View(HK1)
#test distribution
shapiro.test(HK1$Dry_Weight_Below_mg)
hist((HK1$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = HK1)
b<-ggplot(HK1, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(25,110)+
  theme(legend.position="none")+
  ggtitle("HK")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

###############
#####Run 2#####
###############

#Subset for run 2
run2 <- subset(planted, Date_Inoc=="241010")
View(run2)

####Root shoot ratio####

#test distribution
shapiro.test(run2$Root_Shoot_Ratio)
hist((run2$Root_Shoot_Ratio), breaks=10)

#By Genotype and ID
#Plot
p<-ggplot(run2, aes(x=Genotype, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by SynCom ID

#Old
Old2 <- subset(run2, SynCom_ID=="RKMP_Old")
View(Old2)
#test distribution
shapiro.test(Old2$Root_Shoot_Ratio)
hist((Old2$Root_Shoot_Ratio), breaks=5)
kruskal.test(Root_Shoot_Ratio ~ Genotype, data = Old2)
kwAllPairsDunnTest(Old2$Root_Shoot_Ratio~Old2$Genotype, p.adjust.method="BH")
a<-ggplot(Old2, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("Old RKMP")
a

#New
New2 <- subset(run2, SynCom_ID=="RKMP_New")
View(New2)
#test distribution
shapiro.test(New2$Root_Shoot_Ratio)
hist((New2$Root_Shoot_Ratio), breaks=5)
kruskal.test(Root_Shoot_Ratio ~ Genotype, data = New2)
kwAllPairsDunnTest(New2$Root_Shoot_Ratio~New2$Genotype, p.adjust.method="BH")
b<-ggplot(New2, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("New RKMP")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by Genotype

#Morex
Morex2 <- subset(run2, Genotype=="Morex")
View(Morex2)
#test distribution
shapiro.test(Morex2$Root_Shoot_Ratio)
hist((Morex2$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~SynCom_ID, data = Morex2)
#Plot
a<-ggplot(Morex2, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("Morex")
a

#Barke
barke2 <- subset(run2, Genotype=="Barke")
View(barke2)
#test distribution
shapiro.test(barke2$Root_Shoot_Ratio)
hist((barke2$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~SynCom_ID, data = barke2)
#Plot
b<-ggplot(barke2, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("Barke")
b

#124-17
Int172 <- subset(run2, Genotype=="124-17")
View(Int172)
#test distribution
shapiro.test(Int172$Root_Shoot_Ratio)
hist((Int172$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~SynCom_ID, data = Int172)
#Plot
c<-ggplot(Int172, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("124-17")
c

#124-52
Int522 <- subset(run2, Genotype=="124-52")
View(Int522)
#test distribution
shapiro.test(Int522$Root_Shoot_Ratio)
hist((Int522$Root_Shoot_Ratio), breaks=5)
wilcox.test(Root_Shoot_Ratio~SynCom_ID, data = Int522)
#Plot
d<-ggplot(Int522, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,2)+
  theme(legend.position="none")+
  ggtitle("124-52")
d

#faceted graph of all four
grid.arrange(a, b, c, d, ncol=2)

####Above####

#test distribution
shapiro.test(run2$Dry_Weight_Above_mg)
hist((run2$Dry_Weight_Above_mg), breaks=10)

#By Genotype and ID
aov(Dry_Weight_Above_mg~Genotype*SynCom_ID, data=run2)
summary(aov(Dry_Weight_Above_mg~Genotype*SynCom_ID, data=run2))

TukeyHSD(aov(Dry_Weight_Above_mg~Genotype*SynCom_ID, data=run2))
#Plot
p<-ggplot(run2, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by SynCom ID

#Old
Old2 <- subset(run2, SynCom_ID=="RKMP_Old")
View(Old2)
#test distribution
shapiro.test(Old2$Dry_Weight_Above_mg)
hist((Old2$Dry_Weight_Above_mg), breaks=5)
aov(Dry_Weight_Above_mg~Genotype, data=Old2)
summary(aov(Dry_Weight_Above_mg~Genotype, data=Old2))
TukeyHSD(aov(Dry_Weight_Above_mg~Genotype, data=Old2))
a<-ggplot(Old2, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("Old RKMP")
a

#New
New2 <- subset(run2, SynCom_ID=="RKMP_New")
View(New2)
#test distribution
shapiro.test(New2$Dry_Weight_Above_mg)
hist((New2$Dry_Weight_Above_mg), breaks=5)
aov(Dry_Weight_Above_mg~Genotype, data=New2)
summary(aov(Dry_Weight_Above_mg~Genotype, data=New2))
TukeyHSD(aov(Dry_Weight_Above_mg~Genotype, data=New2))
b<-ggplot(New2, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("New RKMP")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by Genotype

#Morex
Morex2 <- subset(run2, Genotype=="Morex")
View(Morex2)
#test distribution
shapiro.test(Morex2$Dry_Weight_Above_mg)
hist((Morex2$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~SynCom_ID, data = Morex2)
#Plot
a<-ggplot(Morex2, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("Morex")
a

#Barke
barke2 <- subset(run2, Genotype=="Barke")
View(barke2)
#test distribution
shapiro.test(barke2$Dry_Weight_Above_mg)
hist((barke2$Dry_Weight_Above_mg), breaks=5)
wilcox.test(Dry_Weight_Above_mg~SynCom_ID, data = barke2)
#Plot
b<-ggplot(barke2, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("Barke")
b

#124-17
Int172 <- subset(run2, Genotype=="124-17")
View(Int172)
#test distribution
shapiro.test(Int172$Dry_Weight_Above_mg)
hist((Int172$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~SynCom_ID, data = Int172)
#Plot
c<-ggplot(Int172, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("124-17")
c

#124-52
Int522 <- subset(run2, Genotype=="124-52")
View(Int522)
#test distribution
shapiro.test(Int522$Dry_Weight_Above_mg)
hist((Int522$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~SynCom_ID, data = Int522)
#Plot
d<-ggplot(Int522, aes(x=SynCom_ID, y=Dry_Weight_Above_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,55)+
  theme(legend.position="none")+
  ggtitle("124-52")
d

#faceted graph of all four
grid.arrange(a, b, c, d, ncol=2)

####Below####

#test distribution
shapiro.test(run2$Dry_Weight_Below_mg)
hist((run2$Dry_Weight_Below_mg), breaks=10)

#By Genotype and ID
#Plot
p<-ggplot(run2, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by SynCom ID

#Old
Old2 <- subset(run2, SynCom_ID=="RKMP_Old")
View(Old2)
#test distribution
shapiro.test(Old2$Dry_Weight_Below_mg)
hist((Old2$Dry_Weight_Below_mg), breaks=5)
kruskal.test(Dry_Weight_Below_mg ~ Genotype, data = Old2)
kwAllPairsDunnTest(Old2$Dry_Weight_Below_mg~Old2$Genotype, p.adjust.method="BH")
a<-ggplot(Old2, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("Old RKMP")
a

#New
New2 <- subset(run2, SynCom_ID=="RKMP_New")
View(New2)
#test distribution
shapiro.test(New2$Dry_Weight_Below_mg)
hist((New2$Dry_Weight_Below_mg), breaks=5)
kruskal.test(Dry_Weight_Below_mg ~ Genotype, data = New2)
kwAllPairsDunnTest(New2$Dry_Weight_Below_mg~New2$Genotype, p.adjust.method="BH")
b<-ggplot(New2, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("New RKMP")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

#Subset by Genotype

#Morex
Morex2 <- subset(run2, Genotype=="Morex")
View(Morex2)
#test distribution
shapiro.test(Morex2$Dry_Weight_Below_mg)
hist((Morex2$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = Morex2)
#Plot
a<-ggplot(Morex2, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("Morex")
a

#Barke
barke2 <- subset(run2, Genotype=="Barke")
View(barke2)
#test distribution
shapiro.test(barke2$Dry_Weight_Below_mg)
hist((barke2$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = barke2)
#Plot
b<-ggplot(barke2, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("Barke")
b

#124-17
Int172 <- subset(run2, Genotype=="124-17")
View(Int172)
#test distribution
shapiro.test(Int172$Dry_Weight_Below_mg)
hist((Int172$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = Int172)
#Plot
c<-ggplot(Int172, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("124-17")
c

#124-52
Int522 <- subset(run2, Genotype=="124-52")
View(Int522)
#test distribution
shapiro.test(Int522$Dry_Weight_Below_mg)
hist((Int522$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~SynCom_ID, data = Int522)
#Plot
d<-ggplot(Int522, aes(x=SynCom_ID, y=Dry_Weight_Below_mg, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,75)+
  theme(legend.position="none")+
  ggtitle("124-52")
d

#faceted graph of all four
grid.arrange(a, b, c, d, ncol=2)

###############
#####Run 3#####
###############

#Subset for run 3
run3 <- subset(planted, Date_Inoc=="241031")
View(run3)

#Root shoot ratio

#test distribution
shapiro.test(run3$Root_Shoot_Ratio)
hist((run3$Root_Shoot_Ratio), breaks=10)

#Plot
p<-ggplot(run3, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

###############
#####Run 4#####
###############

#Subset for run 4
run4 <- subset(planted, Date_Inoc=="241107")
View(run4)

#Root shoot ratio

#test distribution
shapiro.test(run4$Root_Shoot_Ratio)
hist((run4$Root_Shoot_Ratio), breaks=10)

#Plot
p<-ggplot(run4, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

####################
#####Runs 3 & 4#####
####################

#Subset for runs 3 & 4
runs34 <- subset(planted, Date_Inoc!="241003")
runs34 <- subset(runs34, Date_Inoc!="241010")
View(runs34)

####Root shoot ratio####

#test distribution
shapiro.test(runs34$Root_Shoot_Ratio)
hist((runs34$Root_Shoot_Ratio), breaks=10)

####Plot for runs####
p<-ggplot(runs34, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by syncom status

#Live
Live34 <- subset(runs34, Syncom_Status=="Live")
View(Live34)
#test distribution
shapiro.test(Live34$Root_Shoot_Ratio)
hist((Live34$Root_Shoot_Ratio), breaks=5)
aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = Live34)
summary(aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = Live34))
TukeyHSD(aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = Live34))
a<-ggplot(Live34, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Live")
a

#HK
HK34 <- subset(runs34, Syncom_Status=="HK")
View(HK34)
#test distribution
shapiro.test(HK34$Root_Shoot_Ratio)
hist((HK34$Root_Shoot_Ratio), breaks=5)
aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = HK34)
summary(aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = HK34))
TukeyHSD(aov(Root_Shoot_Ratio~Genotype*Date_Inoc, data = HK34))
b<-ggplot(HK34, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Heat-Killed")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

####Plot for syncom status####
p<-ggplot(runs34, aes(x=Genotype, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Subset by genotype

#Morex
Morex34 <- subset(runs34, Genotype=="Morex")
View(Morex34)
#test distribution
shapiro.test(Morex34$Root_Shoot_Ratio)
hist((Morex34$Root_Shoot_Ratio), breaks=5)
wilcox.test(Root_Shoot_Ratio~Syncom_Status, data = Morex34)
#Plot
a<-ggplot(Morex34, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Morex")
a

#Barke
barke34 <- subset(runs34, Genotype=="Barke")
View(barke34)
#test distribution
shapiro.test(barke34$Root_Shoot_Ratio)
hist((barke34$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~Syncom_Status, data = barke34)
#Plot
b<-ggplot(barke34, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("Barke")
b

#124-17
Int1734 <- subset(runs34, Genotype=="124-17")
View(Int1734)
#test distribution
shapiro.test(Int1734$Root_Shoot_Ratio)
hist((Int1734$Root_Shoot_Ratio), breaks=5)
t.test(Root_Shoot_Ratio~Syncom_Status, data = Int1734)
#Plot
c<-ggplot(Int1734, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("124-17")
c

#124-52
Int5234 <- subset(runs34, Genotype=="124-52")
View(Int5234)
#test distribution
shapiro.test(Int5234$Root_Shoot_Ratio)
hist((Int5234$Root_Shoot_Ratio), breaks=5)
wilcox.test(Root_Shoot_Ratio~Syncom_Status, data = Int5234)
#Plot
d<-ggplot(Int5234, aes(x=Syncom_Status, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3)+
  theme(legend.position="none")+
  ggtitle("124-52")
d

#faceted graph of all four
grid.arrange(a, b, c, d, ncol=2)

####Above ground####

#test distribution
shapiro.test(runs34$Dry_Weight_Above_mg)
hist((runs34$Dry_Weight_Above_mg), breaks=10)

aov(Dry_Weight_Above_mg~Genotype*Syncom_Status, data=runs34)
summary(aov(Dry_Weight_Above_mg~Genotype*Syncom_Status, data=runs34))

#Plot
p<-ggplot(runs34, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

####PAPER####
#Subset by genotype

#Morex
mfillcol <- c("#CC79A7", "#FFFFFF")
mcolcol <- c("#000000", "#CC79A7")
Morex34 <- subset(runs34, Genotype=="Morex")
View(Morex34)
#test distribution
shapiro.test(Morex34$Dry_Weight_Above_mg)
hist((Morex34$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~Syncom_Status, data = Morex34)
#Plot
a<-ggplot(Morex34, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = mcolcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,60)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = mfillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
a

#Barke
bfillcol <- c("#0072B2", "#FFFFFF")
bcolcol <- c("#000000", "#0072B2")
barke34 <- subset(runs34, Genotype=="Barke")
View(barke34)
#test distribution
shapiro.test(barke34$Dry_Weight_Above_mg)
hist((barke34$Dry_Weight_Above_mg), breaks=5)
wilcox.test(Dry_Weight_Above_mg~Syncom_Status, data = barke34)
#Plot
b<-ggplot(barke34, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = bcolcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,60)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = bfillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
b

#124-17
fillcol17 <- c("#56B4E9", "#FFFFFF")
colcol17 <- c("#000000", "#56B4E9")
Int1734 <- subset(runs34, Genotype=="124-17")
View(Int1734)
#test distribution
shapiro.test(Int1734$Dry_Weight_Above_mg)
hist((Int1734$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~Syncom_Status, data = Int1734)
#Plot
c<-ggplot(Int1734, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol17)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,60)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol17)+
  scale_x_discrete(labels = c("Live", "HI"))
c

#124-52
fillcol52 <- c("#E69F00", "#FFFFFF")
colcol52 <- c("#000000", "#E69F00")
Int5234 <- subset(runs34, Genotype=="124-52")
View(Int5234)
#test distribution
shapiro.test(Int5234$Dry_Weight_Above_mg)
hist((Int5234$Dry_Weight_Above_mg), breaks=5)
t.test(Dry_Weight_Above_mg~Syncom_Status, data = Int5234)
#Plot
d<-ggplot(Int5234, aes(x=Syncom_Status, y=Dry_Weight_Above_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol52)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,60)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol52)+
  scale_x_discrete(labels = c("Live", "HI"))
d

#faceted graph of all four
grid.arrange(b, c, d, a, ncol=2)

#subset by SynCom Status
#set the colours
Plant_Col <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Live
Live34 <- subset(runs34, Syncom_Status=="Live")
View(Live34)
#test distribution
shapiro.test(Live34$Dry_Weight_Above_mg)
hist((Live34$Dry_Weight_Above_mg), breaks=5)
aov(Dry_Weight_Above_mg~Genotype, data=Live34)
summary(aov(Dry_Weight_Above_mg~Genotype, data=Live34))
TukeyHSD(aov(Dry_Weight_Above_mg~Genotype, data=Live34))
a<-ggplot(Live34, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,50)+
  theme(legend.position="none")+
  ggtitle("Live")+
  scale_fill_manual(values=Plant_Col)
a

#HK
HK34 <- subset(runs34, Syncom_Status=="HK")
View(HK34)
#test distribution
shapiro.test(HK34$Dry_Weight_Above_mg)
hist((HK34$Dry_Weight_Above_mg), breaks=5)
aov(Dry_Weight_Above_mg~Genotype, data=HK34)
summary(aov(Dry_Weight_Above_mg~Genotype, data=HK34))
TukeyHSD(aov(Dry_Weight_Above_mg~Genotype, data=HK34))
b<-ggplot(HK34, aes(x=Genotype, y=Dry_Weight_Above_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(15,50)+
  theme(legend.position="none")+
  ggtitle("HK")+
  scale_fill_manual(values=Plant_Col)
b

#faceted graph of both
grid.arrange(b, a, ncol=2)

####Below ground####

#test distribution
shapiro.test(runs34$Dry_Weight_Below_mg)
hist((runs34$Dry_Weight_Below_mg), breaks=10)

#Plot
p<-ggplot(runs34, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

####PAPER####
#Subset by genotype

#Morex
mfillcol <- c("#CC79A7", "#FFFFFF")
mcolcol <- c("#000000", "#CC79A7")
Morex34 <- subset(runs34, Genotype=="Morex")
View(Morex34)
#test distribution
shapiro.test(Morex34$Dry_Weight_Below_mg)
hist((Morex34$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~Syncom_Status, data = Morex34)
#Plot
a<-ggplot(Morex34, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = mcolcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,120)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = mfillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
a

#Barke
bfillcol <- c("#0072B2", "#FFFFFF")
bcolcol <- c("#000000", "#0072B2")
barke34 <- subset(runs34, Genotype=="Barke")
View(barke34)
#test distribution
shapiro.test(barke34$Dry_Weight_Below_mg)
hist((barke34$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~Syncom_Status, data = barke34)
#Plot
b<-ggplot(barke34, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = bcolcol)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,120)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = bfillcol)+
  scale_x_discrete(labels = c("Live", "HI"))
b

#124-17
fillcol17 <- c("#56B4E9", "#FFFFFF")
colcol17 <- c("#000000", "#56B4E9")
Int1734 <- subset(runs34, Genotype=="124-17")
View(Int1734)
#test distribution
shapiro.test(Int1734$Dry_Weight_Below_mg)
hist((Int1734$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~Syncom_Status, data = Int1734)
#Plot
c<-ggplot(Int1734, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol17)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,120)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol17)+
  scale_x_discrete(labels = c("Live", "HI"))
c

#124-52
fillcol52 <- c("#E69F00", "#FFFFFF")
colcol52 <- c("#000000", "#E69F00")
Int5234 <- subset(runs34, Genotype=="124-52")
View(Int5234)
#test distribution
shapiro.test(Int5234$Dry_Weight_Below_mg)
hist((Int5234$Dry_Weight_Below_mg), breaks=5)
wilcox.test(Dry_Weight_Below_mg~Syncom_Status, data = Int5234)
#Plot
d<-ggplot(Int5234, aes(x=Syncom_Status, y=Dry_Weight_Below_mg, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8), colour = colcol52)+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,120)+
  theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = fillcol52)+
  scale_x_discrete(labels = c("Live", "HI"))
d

#faceted graph of all four
grid.arrange(b, c, d, a, ncol=2)

#Subset live only
Live34 <- subset(runs34, Syncom_Status=="Live")
View(Live34)
#Subset Barke and Morex
BMLive34 <- subset(Live34, Genotype %in% c("Barke", "Morex"))
View(BMLive34)
#test distribution
shapiro.test(BMLive34$Dry_Weight_Below_mg)
hist((BMLive34$Dry_Weight_Below_mg), breaks=5)
t.test(Dry_Weight_Below_mg~Genotype, data = BMLive34)
#Plot
b<-ggplot(BMLive34, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  #ylim(10,110)+
  theme(legend.position="none")+
  ggtitle("Barke vs Morex, Live Only")
b

#subset by SynCom Status
#set the colours
Plant_Col <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Live
Live34 <- subset(runs34, Syncom_Status=="Live")
View(Live34)
#test distribution
shapiro.test(Live34$Dry_Weight_Below_mg)
hist((Live34$Dry_Weight_Below_mg), breaks=5)
aov(Dry_Weight_Below_mg~Genotype, data=Live34)
summary(aov(Dry_Weight_Below_mg~Genotype, data=Live34))
TukeyHSD(aov(Dry_Weight_Below_mg~Genotype, data=Live34))
a<-ggplot(Live34, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,110)+
  theme(legend.position="none")+
  ggtitle("Live")+
  scale_fill_manual(values=Plant_Col)
a

#HK
HK34 <- subset(runs34, Syncom_Status=="HK")
View(HK34)
#test distribution
shapiro.test(HK34$Dry_Weight_Below_mg)
hist((HK34$Dry_Weight_Below_mg), breaks=5)
kruskal.test(Dry_Weight_Below_mg~Genotype, data=HK34)
kwAllPairsDunnTest(Dry_Weight_Below_mg~Genotype, data=HK34, p.adjust.method="BH")
b<-ggplot(HK34, aes(x=Genotype, y=Dry_Weight_Below_mg, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(10,110)+
  theme(legend.position="none")+
  ggtitle("HK")+
  scale_fill_manual(values=Plant_Col)
b

#faceted graph of both
grid.arrange(b, a, ncol=2)

###################
#####All Barke#####
###################

#Subset for Barke
AllBarke <- subset(planted, Genotype=="Barke")
View(AllBarke)

####Root shoot ratio####

#test distribution
shapiro.test(AllBarke$Root_Shoot_Ratio)
hist((AllBarke$Root_Shoot_Ratio), breaks=10)

p<-ggplot(AllBarke, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))
p

#Split by SynCom status

#Live
BarkeLive <- subset(AllBarke, Syncom_Status=="Live")
View(BarkeLive)
#test distribution
shapiro.test(BarkeLive$Root_Shoot_Ratio)
hist((BarkeLive$Root_Shoot_Ratio), breaks=10)
aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeLive)
summary(aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeLive))
TukeyHSD(aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeLive))
#Plot
a<-ggplot(BarkeLive, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3.5)+
  #theme(legend.position="none")+
  ggtitle("Live")
a

#HK
BarkeHK <- subset(AllBarke, Syncom_Status=="HK")
View(BarkeHK)
#test distribution
shapiro.test(BarkeHK$Root_Shoot_Ratio)
hist((BarkeHK$Root_Shoot_Ratio), breaks=10)
aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeHK)
summary(aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeHK))
TukeyHSD(aov(Root_Shoot_Ratio~SynCom_ID*Date_Inoc, data = BarkeHK))
#Plot
b<-ggplot(BarkeHK, aes(x=SynCom_ID, y=Root_Shoot_Ratio, fill=Date_Inoc))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+
  ylim(0,3.5)+
  #theme(legend.position="none")+
  ggtitle("Heat-Killed")
b

#faceted graph of both
grid.arrange(a, b, ncol=2)

