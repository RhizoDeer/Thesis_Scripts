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
#setwd("C:/Users/catar/Documents/Dundee PhD Year 1/SynCom")
getwd()

#Import data
library(readxl)
JH34_Info <- read_excel("JH34_Info.xlsx", 
                        col_types = c("text", "numeric", "numeric", 
                                      "text", "text", "text", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "numeric","text", "text"))
View(JH34_Info)

#Make factors
JH34_Info$Genotype <- as.factor(JH34_Info$Genotype)
JH34_Info$Syncom_Status <- as.factor(JH34_Info$Syncom_Status)
JH34_Info$Opp_Rhizo <- as.factor(JH34_Info$Opp_Rhizo)

#Order factors
JH34_Info$Genotype <- factor(JH34_Info$Genotype, levels = c("Morex", "Barke", "124_17", "124_52", "Bulk"))
JH34_Info$Syncom_Status <- factor(JH34_Info$Syncom_Status, levels = c("Live", "Dead"))

####Separate out the runs####
run1 <- subset(JH34_Info, Date_Inoc==230216)
#View(run1)

run2 <- subset(JH34_Info, Date_Inoc==230223)
#View(run2)

run3 <- subset(JH34_Info, Date_Inoc==230406)
#View(run3)

#####################
########Run 1########
#####################

####Dry Weight####

#Remove bulk and failed samples
plants_only1 <- subset(run1, Genotype!='Bulk'& Sample_ID!='KA30')
View(plants_only1)

#Test for normal distribution
shapiro.test(plants_only1$`Dry_Weight_(mg)`)
hist(plants_only1$`Dry_Weight_(mg)`)

####By genotype####

kruskal.test(`Dry_Weight_(mg)`~Genotype , data=plants_only1)
#kruskal.test(plants_only1$`Dry_Weight_(mg)`~plants_only1$Genotype)

kwAllPairsDunnTest(plants_only1$`Dry_Weight_(mg)`~plants_only1$Genotype, p.adjust.method="BH")

#Plot
r<-ggplot(plants_only1, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####By SynCom status####

wilcox.test(`Dry_Weight_(mg)`~Syncom_Status , data=plants_only1)

#Plot
r<-ggplot(plants_only1, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Syncom_Status") +
  xlab("Syncom_Status") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####By both####

weight_model1 <- lm(`Dry_Weight_(mg)`~Genotype*Syncom_Status, plants_only1)
summary(weight_model1)

#simplify model
simple_weight_model1 <- step(weight_model1)
summary(simple_weight_model1)

AIC(weight_model1)

#individual genotype

barke1 <- subset(plants_only1, Genotype=='Barke')
View(barke1)
shapiro.test(barke1$`Dry_Weight_(mg)`)
hist(barke1$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=barke1)

morex1 <- subset(plants_only1, Genotype=='Morex')
View(morex1)
shapiro.test(morex1$`Dry_Weight_(mg)`)
hist(morex1$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=morex1)

int17_1 <- subset(plants_only1, Genotype=='124_17')
View(int17_1)
shapiro.test(int17_1$`Dry_Weight_(mg)`)
hist(int17_1$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int17_1)

int52_1 <- subset(plants_only1, Genotype=='124_52')
View(int52_1)
shapiro.test(int52_1$`Dry_Weight_(mg)`)
hist(int52_1$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int52_1)

f<-ggplot(plants_only1, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

f  +geom_jitter( size=5,shape=21, position=position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha=2)

f + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Live only

live_plants_only1 <- subset(plants_only1, Syncom_Status=='Live')
View(live_plants_only1)

shapiro.test(live_plants_only1$`Dry_Weight_(mg)`)
hist(live_plants_only1$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only1)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only1))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only1))

#kruskal.test(`Dry_Weight_(mg)`~Genotype , data=live_plants_only1)

#kwAllPairsDunnTest(live_plants_only1$`Dry_Weight_(mg)`~live_plants_only1$Genotype, p.adjust.method="BH")

#Plot
r<-ggplot(live_plants_only1, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Live Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Dead only

dead_plants_only1 <- subset(plants_only1, Syncom_Status=='Dead')
View(dead_plants_only1)

shapiro.test(dead_plants_only1$`Dry_Weight_(mg)`)
hist(dead_plants_only1$`Dry_Weight_(mg)`)

kruskal.test(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only1)

kwAllPairsDunnTest(dead_plants_only1$`Dry_Weight_(mg)`~dead_plants_only1$Genotype, p.adjust.method="BH")

#Plot
r<-ggplot(dead_plants_only1, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Dead Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()


#####################
########Run 2########
#####################

####Dry Weight####

#Remove bulk
plants_only2 <- subset(run2, Genotype!='Bulk')
View(plants_only2)

#Test for normal distribution
shapiro.test(plants_only2$`Dry_Weight_(mg)`)
hist(plants_only2$`Dry_Weight_(mg)`)

####By genotype####

aov(`Dry_Weight_(mg)`~Genotype , data=plants_only2)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=plants_only2))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=plants_only2))

#Plot
r<-ggplot(plants_only2, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

####By SynCom status####

t.test(`Dry_Weight_(mg)`~Syncom_Status , data=plants_only2)

#Plot
r<-ggplot(plants_only2, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Syncom_Status") +
  xlab("Syncom_Status") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667", "#238A8D", "#453781FF"))


####By both####

weight_model2 <- lm(`Dry_Weight_(mg)`~Genotype+Syncom_Status, plants_only2)
summary(weight_model2)

#simplify model
simple_weight_model2 <- step(weight_model2)
summary(simple_weight_model2)

AIC(weight_model2)

aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only2)
summary(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only2))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only2))

#individual genotype

barke2 <- subset(plants_only2, Genotype=='Barke')
View(barke2)
shapiro.test(barke2$`Dry_Weight_(mg)`)
hist(barke2$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=barke2)

b2<-ggplot(barke2, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

b2 + ggtitle("Barke: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

morex2 <- subset(plants_only2, Genotype=='Morex')
View(morex2)
shapiro.test(morex2$`Dry_Weight_(mg)`)
hist(morex2$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=morex2)

m2<-ggplot(morex2, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

m2 + ggtitle("Morex: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

int17_2 <- subset(plants_only2, Genotype=='124_17')
View(int17_2)
shapiro.test(int17_2$`Dry_Weight_(mg)`)
hist(int17_2$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=int17_2)

int172<-ggplot(int17_2, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

int172 + ggtitle("124_17: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

int52_2 <- subset(plants_only2, Genotype=='124_52')
View(int52_2)
shapiro.test(int52_2$`Dry_Weight_(mg)`)
hist(int52_2$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int52_2)

int522<-ggplot(int52_2, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

int522 + ggtitle("124_52: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

f<-ggplot(plants_only2, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

f + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667", "#238A8D", "#453781FF"))

#Live only

live_plants_only2 <- subset(plants_only2, Syncom_Status=='Live')
View(live_plants_only2)
shapiro.test(live_plants_only2$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only2)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only2))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only2))

#Plot
r<-ggplot(live_plants_only2, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Live Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

#Dead only

dead_plants_only2 <- subset(plants_only2, Syncom_Status=='Dead')
View(dead_plants_only2)
shapiro.test(dead_plants_only2$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only2)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only2))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only2))

#Plot
r<-ggplot(dead_plants_only2, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Dead Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

################
####Both 1&2####
################

####dry weight####

plants_only_both <- subset(JH34_Info, Genotype!='Bulk'& Sample_ID!='KA30')
shapiro.test(plants_only_both$`Dry_Weight_(mg)`)

####by genotype####

kruskal.test(`Dry_Weight_(mg)`~Genotype , data=plants_only_both)
kwAllPairsDunnTest(plants_only_both$`Dry_Weight_(mg)`~plants_only_both$Genotype, p.adjust.method="BH")

#Plot
DWxGB<-ggplot(plants_only_both, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
DWxGB + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####by syncom status####

wilcox.test(`Dry_Weight_(mg)`~Syncom_Status , data=plants_only_both)

#Plot
DWxSSB<-ggplot(plants_only_both, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
DWxSSB + ggtitle("Dry Weight by Syncom_Status") +
  xlab("Syncom_Status") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####by both####

weight_modelB <- lm(`Dry_Weight_(mg)`~Genotype*Syncom_Status, plants_only_both)
summary(weight_modelB)

#simplify model
simple_weight_modelB <- step(weight_modelB)
summary(simple_weight_modelB)

AIC(weight_modelB)

#individual genotype
barkeB <- subset(plants_only_both, Genotype=='Barke')
#View(barkeB)
shapiro.test(barkeB$`Dry_Weight_(mg)`)
hist(barkeB$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=barkeB)

morexB <- subset(plants_only_both, Genotype=='Morex')
#View(morexB)
shapiro.test(morexB$`Dry_Weight_(mg)`)
hist(morexB$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=morexB)

int17_B <- subset(plants_only_both, Genotype=='124_17')
#View(int17_B)
shapiro.test(int17_B$`Dry_Weight_(mg)`)
hist(int17_B$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int17_B)

int52_B <- subset(plants_only_both, Genotype=='124_52')
#View(int52_1)
shapiro.test(int52_B$`Dry_Weight_(mg)`)
hist(int52_B$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int52_B)

#Plot
DWxGSSB<-ggplot(plants_only_both, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

DWxGSSB + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Live only

live_plants_only_both <- subset(plants_only_both, Syncom_Status=='Live')
View(live_plants_only_both)

shapiro.test(live_plants_only_both$`Dry_Weight_(mg)`)
hist(live_plants_only_both$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_both)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_both))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_both))

#Plot
r<-ggplot(live_plants_only_both, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Live Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Dead only

dead_plants_only_both <- subset(plants_only_both, Syncom_Status=='Dead')
View(dead_plants_only_both)

shapiro.test(dead_plants_only_both$`Dry_Weight_(mg)`)
hist(dead_plants_only_both$`Dry_Weight_(mg)`)

kruskal.test(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only_both)

kwAllPairsDunnTest(dead_plants_only_both$`Dry_Weight_(mg)`~dead_plants_only_both$Genotype, p.adjust.method="BH")

#Plot
r<-ggplot(dead_plants_only_both, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Dead Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#################
####Run 3########
#################

####Dry Weight####

#Remove bulk
plants_only3 <- subset(run3, Genotype!='Bulk')
View(plants_only3)

#Test for normal distribution
shapiro.test(plants_only3$`Dry_Weight_(mg)`)
hist(plants_only3$`Dry_Weight_(mg)`)

####By genotype####

aov(`Dry_Weight_(mg)`~Genotype , data=plants_only3)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=plants_only3))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=plants_only3))

#Plot
r<-ggplot(plants_only3, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

####By SynCom status####

t.test(`Dry_Weight_(mg)`~Syncom_Status , data=plants_only3)

#Plot
r<-ggplot(plants_only3, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Syncom_Status") +
  xlab("Syncom_Status") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#238A8D", "#453781FF"))

####By both####

aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only3)
summary(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only3))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status , data=plants_only3))

#individual genotype

barke3 <- subset(plants_only3, Genotype=='Barke')
View(barke3)
shapiro.test(barke3$`Dry_Weight_(mg)`)
hist(barke2$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=barke3)

b3<-ggplot(barke3, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

b3 + ggtitle("Barke: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

morex3 <- subset(plants_only3, Genotype=='Morex')
View(morex3)
shapiro.test(morex3$`Dry_Weight_(mg)`)
hist(morex3$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=morex3)

m3<-ggplot(morex3, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

m3 + ggtitle("Morex: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

int17_3 <- subset(plants_only3, Genotype=='124_17')
View(int17_3)
shapiro.test(int17_3$`Dry_Weight_(mg)`)
hist(int17_3$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=int17_3)

int173<-ggplot(int17_3, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

int173 + ggtitle("124_17: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

int52_3 <- subset(plants_only3, Genotype=='124_52')
View(int52_3)
shapiro.test(int52_3$`Dry_Weight_(mg)`)
hist(int52_3$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int52_3)

int523<-ggplot(int52_3, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  ylim(15,50)+
  theme(axis.text = element_text(size=25))

int523 + ggtitle("124_52: Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight(mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "none")

f<-ggplot(plants_only3, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

f + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#238A8D", "#453781FF"))

#Live only

live_plants_only3 <- subset(plants_only3, Syncom_Status=='Live')
View(live_plants_only3)
shapiro.test(live_plants_only3$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only3)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only3))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only3))

#Plot
r<-ggplot(live_plants_only3, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Live Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

#Dead only

dead_plants_only3 <- subset(plants_only3, Syncom_Status=='Dead')
View(dead_plants_only3)
shapiro.test(dead_plants_only3$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only3)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only3))
TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only3))

#Plot
r<-ggplot(dead_plants_only3, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Dead Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

####################
#####Runs 1,2,&3####
####################

####dry weight####

plants_only_123 <- subset(JH34_Info, Genotype!='Bulk'& Sample_ID!='KA30')
shapiro.test(plants_only_123$`Dry_Weight_(mg)`)

plants_only_123$SynComWithStatus <- paste(plants_only_123$SynCom_ID, plants_only_123$Syncom_Status, sep="_")
View(plants_only_123)

####all together####
####by genotype####

kruskal.test(`Dry_Weight_(mg)`~Genotype , data=plants_only_123)
kwAllPairsDunnTest(plants_only_123$`Dry_Weight_(mg)`~plants_only_123$Genotype, p.adjust.method="BH")

#Plot
DWxGB<-ggplot(plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
DWxGB + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####by syncom status####

wilcox.test(`Dry_Weight_(mg)`~Syncom_Status , data=plants_only_123)

#Plot
DWxSSB<-ggplot(plants_only_123, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
DWxSSB + ggtitle("Dry Weight by Syncom_Status") +
  xlab("Syncom_Status") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####by both####

#individual genotype
barke123 <- subset(plants_only_123, Genotype=='Barke')
#View(barke123)
shapiro.test(barke123$`Dry_Weight_(mg)`)
hist(barke123$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=barke123)

morex123 <- subset(plants_only_123, Genotype=='Morex')
#View(morex123)
shapiro.test(morex123$`Dry_Weight_(mg)`)
hist(morex123$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=morex123)

int17_123 <- subset(plants_only_123, Genotype=='124_17')
#View(int17_123)
shapiro.test(int17_123$`Dry_Weight_(mg)`)
hist(int17_123$`Dry_Weight_(mg)`)
t.test(`Dry_Weight_(mg)`~Syncom_Status,data=int17_123)

int52_123 <- subset(plants_only_123, Genotype=='124_52')
#View(int52_123)
shapiro.test(int52_123$`Dry_Weight_(mg)`)
hist(int52_123$`Dry_Weight_(mg)`)
wilcox.test(`Dry_Weight_(mg)`~Syncom_Status,data=int52_123)

#Plot
DWxGSSB<-ggplot(plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

DWxGSSB + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Live only

live_plants_only_123 <- subset(plants_only_123, Syncom_Status=='Live')
View(live_plants_only_123)

shapiro.test(live_plants_only_123$`Dry_Weight_(mg)`)
hist(live_plants_only_123$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_123)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_123))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=live_plants_only_123))

#Plot
r<-ggplot(live_plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Live Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Dead only

dead_plants_only_123 <- subset(plants_only_123, Syncom_Status=='Dead')
View(dead_plants_only_123)

shapiro.test(dead_plants_only_123$`Dry_Weight_(mg)`)
hist(dead_plants_only_123$`Dry_Weight_(mg)`)

aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only_123)
summary(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only_123))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype , data=dead_plants_only_123))

#Plot
r<-ggplot(dead_plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))
r + ggtitle("Dry Weight by Genotype, Dead Only") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####separate by syncom####

#Plot by SynCom and SynCom Status
r<-ggplot(plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=SynComWithStatus))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by Genotype, SynCom, & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

#Plot by SynCom only
r<-ggplot(plants_only_123, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by Genotype & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_viridis_d()

####################
#####Runs 2&3#######
####################

####dry weight####

plants_only_23 <- subset(JH34_Info, Genotype!='Bulk'& Sample_ID!='KA30' & Date_Inoc!=230216)
View(plants_only_23)
shapiro.test(plants_only_23$`Dry_Weight_(mg)`)
hist(plants_only_23$`Dry_Weight_(mg)`, breaks=50)

plants_only_23$SynComWithStatus <- paste(plants_only_23$SynCom_ID, plants_only_23$Syncom_Status, sep="_")
plants_only_23$SynComWithStatus <- factor(plants_only_23$SynComWithStatus, levels = c("RKMP_NoBi27_Live", "RKMP_NoBi27_Dead", "RKMP_Live", "RKMP_Dead"))
plants_only_23$SynCom_ID <- factor(plants_only_23$SynCom_ID, levels = c("RKMP_NoBi27", "RKMP"))
View(plants_only_23)

#by SynCom and genotype

aov(`Dry_Weight_(mg)`~Genotype*SynCom_ID , data=plants_only_23)
summary(aov(`Dry_Weight_(mg)`~Genotype*SynCom_ID , data=plants_only_23))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype*SynCom_ID , data=plants_only_23))

r<-ggplot(plants_only_23, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by Genotype & SynCom") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#238A8D"))

#subset

b12 <- subset(plants_only_23, Genotype=='Barke')
View(b12)
shapiro.test(b12$`Dry_Weight_(mg)`)
wilcox.test(b12$`Dry_Weight_(mg)`~b12$SynCom_ID)
plot(b12$`Dry_Weight_(mg)`~b12$SynCom_ID)

m12 <- subset(plants_only_23, Genotype=='Morex')
View(m12)
shapiro.test(m12$`Dry_Weight_(mg)`)
t.test(m12$`Dry_Weight_(mg)`~m12$SynCom_ID)
plot(m12$`Dry_Weight_(mg)`~m12$SynCom_ID)

int1712 <- subset(plants_only_23, Genotype=='124_17')
View(int1712)
shapiro.test(int1712$`Dry_Weight_(mg)`)
t.test(int1712$`Dry_Weight_(mg)`~int1712$SynCom_ID)
plot(int1712$`Dry_Weight_(mg)`~int1712$SynCom_ID)

int5212 <- subset(plants_only_23, Genotype=='124_52')
View(int5212)
shapiro.test(int5212$`Dry_Weight_(mg)`)
t.test(int5212$`Dry_Weight_(mg)`~int5212$SynCom_ID)
plot(int5212$`Dry_Weight_(mg)`~int5212$SynCom_ID)

#by syncom only
t.test(plants_only_23$`Dry_Weight_(mg)`~plants_only_23$SynCom_ID)

r<-ggplot(plants_only_23, aes(x=SynCom_ID, y=`Dry_Weight_(mg)`, fill=SynCom_ID))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by SynCom") +
  xlab("SynCom") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#238A8D"))

#by genotype only
aov(plants_only_23$`Dry_Weight_(mg)`~plants_only_23$Genotype)
summary(aov(plants_only_23$`Dry_Weight_(mg)`~plants_only_23$Genotype))
TukeyHSD(aov(plants_only_23$`Dry_Weight_(mg)`~plants_only_23$Genotype))

r<-ggplot(plants_only_23, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=Genotype))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by Genotype") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#CC79A7", "#0072B2", "#56B4E9", "#E69F00", "#FFFFFF"))

#by syncom status only
t.test(plants_only_23$`Dry_Weight_(mg)`~plants_only_23$Syncom_Status)

r<-ggplot(plants_only_23, aes(x=Syncom_Status, y=`Dry_Weight_(mg)`, fill=Syncom_Status))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by SynCom Status") +
  xlab("SynCom Status") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667"))

#all together
aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status*SynCom_ID , data=plants_only_23)
summary(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status*SynCom_ID , data=plants_only_23))

TukeyHSD(aov(`Dry_Weight_(mg)`~Genotype*Syncom_Status*SynCom_ID , data=plants_only_23))

r<-ggplot(plants_only_23, aes(x=Genotype, y=`Dry_Weight_(mg)`, fill=SynComWithStatus))+ 
  geom_boxplot(position=position_dodge(0.8))+ 
  geom_jitter(position=position_dodge(0.8))+ 
  #ylim(0,700)+
  theme(axis.text.x = element_text(size=14, angle=45, hjust=1))

r + ggtitle("Dry Weight by Genotype, SynCom, & SynCom Status") +
  xlab("Genotype") + ylab("Dry Weight (mg)") +
  scale_fill_manual(values=c("#FDE725", "#55C667", "#238A8D", "#453781FF"))