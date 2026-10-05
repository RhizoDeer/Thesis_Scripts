#####################################################################################
#BETA DIVERSITY CALCULATION
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in Katie's work
#  Revision 01/24
#  carnton@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
library("phyloseq")
library("vegan")
library("ggplot2")
#############################################################

#set working directory: DB cpu
#setwd("/cluster/db/R_shared/JH201/")

#set working directory: Katie
setwd("/cluster/db/carnton/JH201")

#check working directory
getwd()

#import the datasets
JH201_tube <-readRDS("Original_data/JH201_tube_rare.rds")

#inspect the files
JH201_tube
sample_data(JH201_tube)

#remove control and dead samples
JH201_tube_nocontrol <- subset_samples(JH201_tube, SynCom.ID != 'SDW')
JH201 <- subset_samples(JH201_tube_nocontrol, SynCom.ID != 'PBS')
JH201 <- subset_samples(JH201, Syncom.Status != 'Dead')
JH201
sample_data(JH201)

#Convert date variables to factors
sample_data(JH201)$Date.Harvest <- as.factor(sample_data(JH201)$Date.Harvest)
sample_data(JH201)$Date.Inoc <- as.factor(sample_data(JH201)$Date.Inoc)

###################################################################################################
#Beta diversity calculation
# Rationale: we are deling with a multivariate dataset, i.e, the number of dependent variables (i.e., observations) is greater
# than the independent ones (i.e., the individual factors in our experimental design)
# This is common theme in ecology, and we use ecological theories to summarise this diversity into a bidimensional space
# Effectively we will convert microbiota features, such as ASVs presence/absence and abundance, in distances.
# This has two advantages: the first is provide us (and your readers) with an immediated overview of community relatedness
# The second one is that those "distances" can be treated statistically and this will allow us to infer whether one of the independent variables has more impact than others
# This manuscript is a must read prior delving with the code: https://academic.oup.com/femsec/article/62/2/142/434668?login=true 
# Another useful manuscript is this one dedicated to constrained ordinations, which we refer to "CAP", one of the approaches we will be implementing:https://doi.org/10.1890/0012-9658(2003)084[0511:CAOPCA]2.0.CO;2
# More info for options in phyloseq here: https://joey711.github.io/phyloseq/distance.html
###################################################################################################

###################################################################################################
#arrange data for plotting
#note that here we have two factors
#1) the SynCom
#2) the genotype
##################################################################################################

#SynCom
sample_data(JH201)$SynCom.ID <- factor(sample_data(JH201)$SynCom.ID, levels=c("RKMP", "RKMP.NoBi27"))

#Genotype
sample_data(JH201)$Genotype <- factor(sample_data(JH201)$Genotype, levels=c("Bulk", "Barke", "SL17", "SL52", "Morex"))

#color coding
DB_cols <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

###################################################################################################
#First type of ordination, unconstrained one: here we ask the data to inform us on sequence relatedness
#The code has two elements: the generation of the distance and its visualisation, respectively
###################################################################################################

#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
JH201.ord <- ordinate(JH201, "NMDS", "bray")
#note the elements of the function: 1) a phyloseq object, 2) the type of ordination we want to implement and 3) the type of distance

#now the visualisation
plot_ordination(JH201, JH201.ord, type="samples", color="Genotype", shape="SynCom.ID")

#ggplots function to increase effectivness of the visualisation
p = plot_ordination(JH201, JH201.ord, type="samples", color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("All live samples, unconstrained")

###################################################################################################
#Statistical analysis: here we ask two questions. 
#First, on the basis of microbiota composition, can we statistically discriminate level of one or more factors?
#Second, if so, what is the proportion of variance explained by that factor(s)?
#We apply a permutational analysis of variance (see attached powerpoint presentation)
###################################################################################################

#this calculation should use the same distance used to build the graphical output
BC <- phyloseq::distance(JH201, "bray")
BC

#here we can use a formula ANOVA-like to identify factors of interest
Stat <- adonis2(BC ~ Genotype * SynCom.ID, data= as.data.frame(as.matrix(sample_data(JH201))), permutations = 5000)
#Stat <- adonis2(BC ~ SynCom.ID * Genotype, data= as.data.frame(as.matrix(sample_data(JH201))), permutations = 5000)
#Exact results are dependent on variable order in the model, though the differences are minimal -> pick one way round and stick to it
Stat
#note the script
#first a formula
#second where to source who is who
#third the number of permutations


###################################################################################################
#Second type of ordination, constrained one: here we ask to maximise the differnces as function of a given variable
#The code will have the same structure, with one important difference in the distance calculation
###################################################################################################

JH201.cap <- ordinate(JH201, "CAP", "bray", ~ Genotype * SynCom.ID)
#note the formula to specify what factor(s) look for

plot_ordination(JH201, JH201.cap, color="Genotype", shape="SynCom.ID")

#ggplots function to increase effectivness of the visualisation
p = plot_ordination(JH201, JH201.cap,  color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("All live samples, constrained")

#########################################################
####Repeat analysis after subsetting for planted only####
#########################################################

#subset
JH201_planted <- subset_samples(JH201, Genotype != "Bulk")
sample_data(JH201_planted)

#modify colour scheme to remove bulk
DB_cols_planted <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#unconstrained ordination
JH201_planted.ord <- ordinate(JH201_planted, "NMDS", "bray")

#visualisation with SynCom
p = plot_ordination(JH201_planted, JH201_planted.ord, type="samples", color="Genotype", shape="SynCom.ID")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("Planted live samples, unconstrained")

#visualisation with date inoculated
p = plot_ordination(JH201_planted, JH201_planted.ord, type="samples", color="Genotype", shape="Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("Planted live samples, unconstrained")

#Statistical analysis
BC_planted <- phyloseq::distance(JH201_planted, "bray")
BC_planted
Stat <- adonis2(BC_planted ~ Genotype * SynCom.ID * Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_planted))), permutations = 5000)
Stat

#constrained ordination
JH201_planted.cap <- ordinate(JH201_planted, "CAP", "bray", ~ Genotype * SynCom.ID)

#visualisation - SynCom
p = plot_ordination(JH201_planted, JH201_planted.cap,  color="Genotype", shape="SynCom.ID") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("Planted live samples, constrained")

#visualisation - Run
p = plot_ordination(JH201_planted, JH201_planted.cap,  color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("Planted live samples, constrained")

###################################################
####Repeat analysis after subsetting for SynCom####
###################################################

####RKMP####

#Subset
JH201_RKMP <- subset_samples(JH201, SynCom.ID=="RKMP")
sample_data(JH201_RKMP)

#unconstrained ordination - using run as a variable instead of SynCom as there's only one SynCom here
#Using date inoculated for now as I've spotted some errors in harvest date for some samples **FIX THIS**
JH201_RKMP.ord <- ordinate(JH201_RKMP, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_RKMP, JH201_RKMP.ord, type="samples", color="Genotype", shape="Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("NMDS 16S data Tubes, Bray distance")

#Statistical analysis
BC_RKMP <- phyloseq::distance(JH201_RKMP, "bray")
BC_RKMP
Stat <- adonis2(BC_RKMP ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_RKMP))), permutations = 5000)
Stat

#constrained ordination
JH201_RKMP.cap <- ordinate(JH201_RKMP, "CAP", "bray", ~ Genotype * Date.Inoc)
#visualisation
p = plot_ordination(JH201_RKMP, JH201_RKMP.cap,  color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 16, 17))
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("CAP 16S data, Bray distance")

####RKMP.NoBi27####

#Subset
JH201_RKMP.NoBi27 <- subset_samples(JH201, SynCom.ID=="RKMP.NoBi27")
sample_data(JH201_RKMP.NoBi27)

#unconstrained ordination - using run as a variable instead of SynCom as there's only one SynCom here
#Using date inoculated for now as I've spotted some errors in harvest date for some samples **FIX THIS**
JH201_RKMP.NoBi27.ord <- ordinate(JH201_RKMP.NoBi27, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_RKMP.NoBi27, JH201_RKMP.NoBi27.ord, type="samples", color="Genotype", shape="Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("NMDS 16S data Tubes, Bray distance")

#Statistical analysis
BC_RKMP.NoBi27 <- phyloseq::distance(JH201_RKMP.NoBi27, "bray")
BC_RKMP.NoBi27
Stat <- adonis2(BC_RKMP.NoBi27 ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_RKMP.NoBi27))), permutations = 5000)
Stat

#constrained ordination
JH201_RKMP.NoBi27.cap <- ordinate(JH201_RKMP.NoBi27, "CAP", "bray", ~ Genotype * Date.Inoc)
#visualisation
p = plot_ordination(JH201_RKMP.NoBi27, JH201_RKMP.NoBi27.cap,  color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 16, 17))
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("CAP 16S data, Bray distance")

####################################################################
####Repeat analysis after subsetting for planted only and SynCom####
####################################################################

####RKMP####

#subset
JH201_planted_RKMP <- subset_samples(JH201_planted, SynCom.ID == "RKMP")
sample_data(JH201_planted_RKMP)

#unconstrained ordination
JH201_planted_RKMP.ord <- ordinate(JH201_planted_RKMP, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_planted_RKMP, JH201_planted_RKMP.ord, type="samples", color="Genotype", shape="Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("NMDS 16S data Tubes, Bray distance")

#Statistical analysis
BC_planted_RKMP <- phyloseq::distance(JH201_planted_RKMP, "bray")
BC_planted_RKMP
Stat <- adonis2(BC_planted_RKMP ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_planted_RKMP))), permutations = 5000)
Stat

#constrained ordination
JH201_planted_RKMP.cap <- ordinate(JH201_planted_RKMP, "CAP", "bray", ~ Genotype * Date.Inoc)
#visualisation
p = plot_ordination(JH201_planted_RKMP, JH201_planted_RKMP.cap,  color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 16, 17))
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("CAP 16S data, Bray distance")

####RKMP.NoBi27####

#subset
JH201_planted_RKMP.NoBi27 <- subset_samples(JH201_planted, SynCom.ID == "RKMP.NoBi27")
sample_data(JH201_planted_RKMP.NoBi27)

#unconstrained ordination
JH201_planted_RKMP.NoBi27.ord <- ordinate(JH201_planted_RKMP.NoBi27, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_planted_RKMP.NoBi27, JH201_planted_RKMP.NoBi27.ord, type="samples", color="Genotype", shape="Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("NMDS 16S data Tubes, Bray distance")

#Statistical analysis
BC_planted_RKMP.NoBi27 <- phyloseq::distance(JH201_planted_RKMP.NoBi27, "bray")
BC_planted_RKMP.NoBi27
Stat <- adonis2(BC_planted_RKMP.NoBi27 ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_planted_RKMP.NoBi27))), permutations = 5000)
Stat

#constrained ordination
JH201_planted_RKMP.NoBi27.cap <- ordinate(JH201_planted_RKMP.NoBi27, "CAP", "bray", ~ Genotype * Date.Inoc)
#visualisation
p = plot_ordination(JH201_planted_RKMP.NoBi27, JH201_planted_RKMP.NoBi27.cap,  color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 16, 17))
p = p + scale_colour_manual(values = DB_cols_planted)
p + ggtitle("CAP 16S data, Bray distance")

#####################
####Subset by Run####
#####################

####RKMP and all genotypes####

#subset
JH201_r3 <- subset_samples(JH201, Date.Inoc == "230406")
sample_data(JH201_r3)

#unconstrained ordination
JH201_r3.ord <- ordinate(JH201_r3, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r3, JH201_r3.ord, type="samples", color="Genotype")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("RKMP and all genotypes, unconstrained")

#Statistical analysis
BC_r3 <- phyloseq::distance(JH201_r3, "bray")
BC_r3
#Using genotype
Stat <- adonis2(BC_r3 ~ Genotype, data= as.data.frame(as.matrix(sample_data(JH201_r3))), permutations = 5000)
Stat
#Using plant
Stat <- adonis2(BC_r3 ~ Plant, data= as.data.frame(as.matrix(sample_data(JH201_r3))), permutations = 5000)
Stat

####THESIS####
#constrained ordination
JH201_r3.cap <- ordinate(JH201_r3, "CAP", "bray", ~ Genotype)
#visualisation
p = plot_ordination(JH201_r3, JH201_r3.cap,  color="Genotype") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols, labels = c("Bulk", "Barke", "124-17", "124-52", "Morex"))
p + theme(legend.position="bottom")

####RKMP_NoBi27 and all genotypes####

#subset
JH201_r2 <- subset_samples(JH201, Date.Inoc == "230223")
sample_data(JH201_r2)
#only one bulk sample - two were lost to the pilot and two were lost before rarifying due to low read count

#unconstrained ordination
JH201_r2.ord <- ordinate(JH201_r2, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r2, JH201_r2.ord, type="samples", color="Genotype")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("RKMP_NoBi27 and all genotypes unconstrained")

#Statistical analysis
BC_r2 <- phyloseq::distance(JH201_r2, "bray")
BC_r2
Stat <- adonis2(BC_r2 ~ Genotype, data= as.data.frame(as.matrix(sample_data(JH201_r2))), permutations = 5000)
Stat

#constrained ordination
JH201_r2.cap <- ordinate(JH201_r2, "CAP", "bray", ~ Genotype)
#visualisation
p = plot_ordination(JH201_r2, JH201_r2.cap,  color="Genotype") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 16, 17))
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("RKMP_NoBi27 and all genotypes constrained")

####Barke and both SynComs####

#subset
JH201_r4 <- subset_samples(JH201, Date.Inoc == "230525")
sample_data(JH201_r4)

#modify colour scheme to remove all bar barke and bulk
DB_cols_BB <- c("#0072B2", "#000000")

#unconstrained ordination
JH201_r4.ord <- ordinate(JH201_r4, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r4, JH201_r4.ord, type="samples", color="Genotype", shape = "SynCom.ID")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_BB)
p + ggtitle("Barke and both SynComs, Unconstrained")

#Statistical analysis
BC_r4 <- phyloseq::distance(JH201_r4, "bray")
BC_r4
Stat <- adonis2(BC_r4 ~ Genotype*SynCom.ID, data= as.data.frame(as.matrix(sample_data(JH201_r4))), permutations = 5000)
Stat

#constrained ordination
JH201_r4.cap <- ordinate(JH201_r4, "CAP", "bray", ~ Genotype*SynCom.ID)
#visualisation
p = plot_ordination(JH201_r4, JH201_r4.cap,  color="Genotype", shape = "SynCom.ID") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols_BB)
p + ggtitle("Barke and both SynComs, constrained")

#stats for RKMP_NoBi27 only
#subset
JH201_r4_No27 <- subset_samples(JH201_r4, SynCom.ID == "RKMP.NoBi27")
sample_data(JH201_r4_No27)
#Statistical analysis
BC_r4_No27 <- phyloseq::distance(JH201_r4_No27, "bray")
BC_r4_No27
Stat <- adonis2(BC_r4_No27 ~ Plant, data= as.data.frame(as.matrix(sample_data(JH201_r4_No27))), permutations = 5000)
Stat

####################################
####Subset by Run - planted only####
####################################

####RKMP and all genotypes####

#subset
JH201_r3_planted <- subset_samples(JH201_planted, Date.Inoc == "230406")
sample_data(JH201_r3_planted)

#unconstrained ordination
JH201_r3_planted.ord <- ordinate(JH201_r3_planted, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r3_planted, JH201_r3_planted.ord, type="samples", color="Genotype")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("RKMP and all genotypes, planted, unconstrained")

#Statistical analysis
BC_r3_planted <- phyloseq::distance(JH201_r3_planted, "bray")
BC_r3_planted
Stat <- adonis2(BC_r3_planted ~ Genotype, data= as.data.frame(as.matrix(sample_data(JH201_r3_planted))), permutations = 5000)
Stat

####THESIS####
#constrained ordination
JH201_r3_planted.cap <- ordinate(JH201_r3_planted, "CAP", "bray", ~ Genotype)
#visualisation
p = plot_ordination(JH201_r3_planted, JH201_r3_planted.cap,  color="Genotype") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols_planted, labels = c("Barke", "124-17", "124-52", "Morex"))
p + theme(legend.position="bottom")


####RKMP_NoBi27 and all genotypes####

#subset
JH201_r2_planted <- subset_samples(JH201_planted, Date.Inoc == "230223")
sample_data(JH201_r2_planted)

#unconstrained ordination
JH201_r2_planted.ord <- ordinate(JH201_r2_planted, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r2_planted, JH201_r2_planted.ord, type="samples", color="Genotype")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols)
p + ggtitle("RKMP_NoBi27 and all genotypes, planted, unconstrained")

#Statistical analysis
BC_r2_planted <- phyloseq::distance(JH201_r2_planted, "bray")
BC_r2_planted
Stat <- adonis2(BC_r2_planted ~ Genotype, data= as.data.frame(as.matrix(sample_data(JH201_r2_planted))), permutations = 5000)
Stat

####THESIS####
#constrained ordination
JH201_r2_planted.cap <- ordinate(JH201_r2_planted, "CAP", "bray", ~ Genotype)
#visualisation
p = plot_ordination(JH201_r2_planted, JH201_r2_planted.cap,  color="Genotype") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols_planted, labels = c("Barke", "124-17", "124-52", "Morex"))
p + theme(legend.position="bottom")

####Barke and both SynComs####

#subset
JH201_r4_planted <- subset_samples(JH201_planted, Date.Inoc == "230525")
sample_data(JH201_r4_planted)

#modify colour scheme to remove all bar barke and bulk
DB_cols_BB <- c("#0072B2", "#000000")

#unconstrained ordination
JH201_r4_planted.ord <- ordinate(JH201_r4_planted, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_r4_planted, JH201_r4_planted.ord, type="samples", color="Genotype", shape = "SynCom.ID")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_BB)
p + ggtitle("Barke and both SynComs, planted, unconstrained")

#Statistical analysis
BC_r4_planted <- phyloseq::distance(JH201_r4_planted, "bray")
BC_r4_planted
Stat <- adonis2(BC_r4_planted ~ SynCom.ID, data= as.data.frame(as.matrix(sample_data(JH201_r4_planted))), permutations = 5000)
Stat

#constrained ordination
JH201_r4_planted.cap <- ordinate(JH201_r4_planted, "CAP", "bray", ~ SynCom.ID)
#visualisation
p = plot_ordination(JH201_r4_planted, JH201_r4_planted.cap,  color="Genotype", shape = "SynCom.ID") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = DB_cols_BB)
p + ggtitle("Barke and both SynComs, planted, constrained")

##################
####Run Effect####
##################

####Barke and RKMP####

#subset
JH201_Barke_RKMP <- subset_samples(JH201_planted, Genotype == "Barke")
JH201_Barke_RKMP <- subset_samples(JH201_Barke_RKMP, SynCom.ID == "RKMP")
sample_data(JH201_Barke_RKMP)

#modify colour scheme to remove all bar barke
DB_cols_Ba <- c("#0072B2")

#unconstrained ordination
JH201_Barke_RKMP.ord <- ordinate(JH201_Barke_RKMP, "NMDS", "bray")
#visualisation
p = plot_ordination(JH201_Barke_RKMP, JH201_Barke_RKMP.ord, type="samples", color="Genotype", shape = "Date.Inoc")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_Ba)
p + ggtitle("Barke and RKMP, planted, unconstrained")

#Statistical analysis
BC_JH201_Barke_RKMP <- phyloseq::distance(JH201_Barke_RKMP, "bray")
BC_JH201_Barke_RKMP
Stat <- adonis2(BC_JH201_Barke_RKMP ~ Date.Inoc, data= as.data.frame(as.matrix(sample_data(JH201_Barke_RKMP))), permutations = 5000)
Stat

#constrained ordination
JH201_Barke_RKMP.cap <- ordinate(JH201_Barke_RKMP, "CAP", "bray", ~ Date.Inoc)
#visualisation
p = plot_ordination(JH201_Barke_RKMP, JH201_Barke_RKMP.cap,  color="Genotype", shape = "Date.Inoc") 
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = DB_cols_Ba)
p + ggtitle("Barke and RKMP, planted, constrained")

###################################################################################################
#Key points to remember
###################################################################################################

#1) Microbiota data are multivariate datasets
#2) Main differences between constrained and unconstrained ordination
#3) Ordinations are an effective to visualise data, not a statistical analysis per se 

###################################################################################################
