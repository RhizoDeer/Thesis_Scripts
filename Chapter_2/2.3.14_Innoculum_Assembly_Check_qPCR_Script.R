#####################################################################################
#Univariate stats
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in Katie's TC report 08/24
#  Revision 08/24
#  c.arnton@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk
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
library("ggplot2")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Katie/")
#set working directory-Katie CPU
#TO BE COMPLETED
#############################################################

#import the dataset
KA_info <-read.delim("KA_qPCR_0824_2.txt", row.names = 1, header = T)

#inspect the files
KA_info

#check the type of columns: independent variables, e.g., treatment, should be factor while dependent variables, e.g., biomass, should be numeric 
#Independent variables
class(KA_info$SynCom)
class(KA_info$Preparation)
#Dependent variables
class(KA_info$Ct)

#subset for SynCom
KA_info_SynCom <- subset(KA_info, SynCom != "NTC")
KA_info_SynCom

#convert to factor
KA_info_SynCom$SynCom <- as.factor(KA_info_SynCom$SynCom)
class(KA_info_SynCom$SynCom)

KA_info_SynCom$Preparation <- as.factor(KA_info_SynCom$Preparation)
class(KA_info_SynCom$Preparation)

###################################################################################################
#Let's re-order factors and colors
###################################################################################################

#SynCom
levels(KA_info_SynCom$SynCom)
KA_info_SynCom$SynCom <- factor(KA_info_SynCom$SynCom, levels=c("Bacillus+", "Bacillus-"))

#Preparations
levels(KA_info_SynCom$Preparation)
KA_info_SynCom$Preparation <- factor(KA_info_SynCom$Preparation, levels=c("P1", "P2"))

#color coding
DB_cols <- c("#56B4E9", "#D55E00")

###################################################################################################
#Ct data visualization using boxplot 
###################################################################################################
#useful website http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization 

#Type effect
p <- ggplot(KA_info_SynCom, aes(x=SynCom, y=Ct)) + 
  geom_boxplot() +
  facet_wrap( ~Preparation)
p

dev.off()

#our colours
p <- ggplot(KA_info_SynCom, aes(x=SynCom, y=Ct, color = SynCom)) + 
  geom_boxplot(position=position_dodge(0.8)) + scale_color_manual(values = DB_cols) +
  facet_wrap( ~Preparation, labeller = labeller(Preparation = c("P1"="Agitated","P2"="Settled")))+
  theme(legend.position = "none")+ 
  geom_jitter(position=position_dodge(0.8))
p

#formal assessment of data distribution to identify the appropriate test 
hist(KA_info_SynCom$Ct)
#Shapiro test: https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test
shapiro.test(KA_info_SynCom$Ct)
# p value above .05 denotes data distribution comparable to a  normal one, use parametric tests

#statistical test: using a two-way anova
#general formula before the sign goes the name fo the column of dependent variable, after the name of the column of independent you want to test
res.aov <-aov(Ct ~ SynCom * Preparation, data = KA_info_SynCom)
summary(res.aov)

#identify significant differences
TukeyHSD(res.aov)

#end