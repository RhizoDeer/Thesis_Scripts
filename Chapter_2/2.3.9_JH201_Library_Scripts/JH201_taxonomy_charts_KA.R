#####################################################################################
#TAXONOMY CHARTS
####################################################################################
#############################################################
#
# Ref to the ARTICLE 
# 
#  Code to compute calculations presented in Katie's work
#  Revision 02/24
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
library("viridis")

#############################################################
#set working directory: DB cpu
#setwd("/cluster/db/R_shared/JH201/")
#set working directory: Katie
setwd("/cluster/db/carnton/JH201")
#check working directory
getwd()

#import the datasets
#JH201 <-readRDS("JH201_genus_50_reads.rds")
JH201 <-readRDS("Original_data/JH201_genus_50_reads.rds")

#inspect the files
JH201
sample_data(JH201)

#summary
min(sample_sums(JH201))
max(sample_sums(JH201))
mean(sample_sums(JH201))

#split for sample type
#inoculants
JH201_inoculants <- subset_samples(JH201, Sample.Type ==  "Inoculant")
JH201_inoculants
#reads distribution
sort(sample_sums(JH201_inoculants))

#non inoculants
JH201_tube <- subset_samples(JH201, Sample.Type ==  "Tube")
JH201_tube
#reads distribution
sort(sample_sums(JH201_tube))

##############################################################################################################################################
# Generate a new phyloseq object for calculation: the problem is that we do have differences in sequencing depth (i.e., some samples have more some have less)
# A solution could be "down-sampling" everything to the lowest count, a process called rarefying. 
# However, rarefying is a "contentious" matter in microbe science
# https://doi.org/10.1371/journal.pcbi.1003531 
# as we have a nearly 10x difference between the sample with the largest number of reads and the one with the lowest we will rarefy the phyloseq object
# https://doi.org/10.1186/s40168-017-0237-y
#############################################################################################################################################

#rarefy at 10K
JH201_inoculants_rare <- rarefy_even_depth(JH201_inoculants, 1500)
#ignore the warning, the object will be saved
JH201_inoculants_rare

###################################################################################################
#save the rarefied object for reproducibility of the code
#saveRDS(JH201_inoculants_rare, file = "JH201_inoculants_rares.rds")
###################################################################################################

#load saved rare
#JH201_inoculants_rare <-readRDS("Original_data/JH201_inoculants_rares.rds")

#Plotting inoculant composition
p <- plot_bar(JH201_inoculants_rare, fill="Genus")
p + scale_fill_viridis(discrete=TRUE)

#Task: double-check for consistency with predicted inoculum

#re-arrange order
sample_data(JH201_inoculants_rare)$SynCom.ID <- factor(sample_data(JH201_inoculants_rare)$SynCom.ID, levels=c("RKMP", "RKMP.NoBi27"))
sample_data(JH201_inoculants_rare)$Syncom.Status <- factor(sample_data(JH201_inoculants_rare)$Syncom.Status, levels=c("Live", "Dead"))

#faceted version
sample_data(JH201_inoculants_rare)
p <- plot_bar(JH201_inoculants_rare, "Genus", fill="Genus", facet_grid=SynCom.ID~Syncom.Status)
p + scale_fill_viridis(discrete=TRUE)

#Remove KA-IR5-1, its the only live inoculant which doesn't have a direct 'dead' pair
JH201_inoculants_rare_noR5 <- subset_samples(JH201_inoculants_rare, Date.Inoc != "230831")
sample_data(JH201_inoculants_rare_noR5)
#Plot
p <- plot_bar(JH201_inoculants_rare_noR5, "Genus", fill="Genus", facet_grid=SynCom.ID~Syncom.Status)
p + scale_fill_viridis(discrete=TRUE)

#subset for bacillus
JH201_inoculants_rare_bacillus <- subset_taxa(JH201_inoculants_rare, Genus == "Bacillus")
p <- plot_bar(JH201_inoculants_rare_bacillus)
p + scale_fill_viridis(discrete=TRUE)

#subset for RKMP Live only
JH201_inoculants_rare_RKMP_live <- subset_samples(JH201_inoculants_rare, SynCom.ID == "RKMP")
JH201_inoculants_rare_RKMP_live <- subset_samples(JH201_inoculants_rare_RKMP_live, Syncom.Status == "Live")
sample_data(JH201_inoculants_rare_RKMP_live)
#plot
p <- plot_bar(JH201_inoculants_rare_RKMP_live, fill="Genus")
p + scale_fill_viridis(discrete=TRUE)

##############################################################
####Re-try with the separately pre-processed inoculum data####
##############################################################

JH201_inoculants2 <-readRDS("Original_data/JH201_inoculant_genus_10r20p.rds")
JH201_inoculants2

#rarefy at 1500
JH201_inoculants2_rare <- rarefy_even_depth(JH201_inoculants2, 1500)
#ignore the warning, the object will be saved
JH201_inoculants2_rare

#save the rarefied object for reproducibility of the code
#saveRDS(JH201_inoculants2_rare, file = "Original_data/JH201_inoculants2_rares.rds")

#re-arrange order
sample_data(JH201_inoculants2_rare)$SynCom.ID <- factor(sample_data(JH201_inoculants2_rare)$SynCom.ID, levels=c("RKMP", "RKMP.NoBi27"))
sample_data(JH201_inoculants2_rare)$Syncom.Status <- factor(sample_data(JH201_inoculants2_rare)$Syncom.Status, levels=c("Live", "Dead"))

#Plotting inoculant composition
p <- plot_bar(JH201_inoculants2_rare, fill="Genus")
p + scale_fill_viridis(discrete=TRUE)

sample_data(JH201_inoculants2_rare)
p <- plot_bar(JH201_inoculants2_rare, "Genus", fill="Genus", facet_grid=SynCom.ID~Syncom.Status)
p + scale_fill_viridis(discrete=TRUE)

####USED IN THESIS####
#Remove KA-IR5-1, its the only live inoculant which doesn't have a direct 'dead' pair
JH201_inoculants2_rare_noR5 <- subset_samples(JH201_inoculants2_rare, Date.Inoc != "230831")
sample_data(JH201_inoculants2_rare_noR5)
#Plot
ID_labels <- c("RKMP"="RKMP", "RKMP.NoBi27"="RKMP-27")
status_labels <- c("Live"="Live","Dead"="HI")
p <- plot_bar(JH201_inoculants2_rare_noR5, "Genus", fill="Genus")
p + scale_fill_viridis(discrete=TRUE)+facet_grid(SynCom.ID~Syncom.Status, labeller = labeller(SynCom.ID = ID_labels, Syncom.Status = status_labels))+
  theme(legend.position = "right", axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

#subset for bacillus
JH201_inoculants2_rare_bacillus <- subset_taxa(JH201_inoculants2_rare, Genus == "Bacillus")
p <- plot_bar(JH201_inoculants2_rare_bacillus)
p + scale_fill_viridis(discrete=TRUE)

#subset for RKMP Live only
JH201_inoculants2_rare_RKMP_live <- subset_samples(JH201_inoculants2_rare, SynCom.ID == "RKMP")
JH201_inoculants2_rare_RKMP_live <- subset_samples(JH201_inoculants2_rare_RKMP_live, Syncom.Status == "Live")
sample_data(JH201_inoculants2_rare_RKMP_live)
#plot
p <- plot_bar(JH201_inoculants2_rare_RKMP_live, fill="Genus")
p + scale_fill_viridis(discrete=TRUE)

###################################################################################################
#arrange data for plotting tubes
#a) check the low count samples
sort(sample_sums(JH201_tube))
sample_data(JH201_tube)
#retain samples with 15,000 reads or more
JH201_10K <- prune_samples(sample_sums(JH201_tube) > 10000, JH201_tube)
JH201_10K
sort(sample_sums(JH201_10K))

#rarefy at 10K
JH201_tube_rare <- rarefy_even_depth(JH201_10K, 10000)
#ignore the warning, the object will be saved
JH201_tube_rare

###################################################################################################
#save the rarefied object for reproducibility of the code
#saveRDS(JH201_tube_rare, file = "JH201_tube_rare.rds")
###################################################################################################

#load saved rare
#JH201_tube_rare <-readRDS("Original_data/JH201_tube_rare.rds")

sample_data(JH201_tube_rare)
#subest for live samples
JH201_tube_rare_live <- subset_samples(JH201_tube_rare, Syncom.Status == "Live")

#re-arrange samples
#SynCom
sample_data(JH201_tube_rare_live)$SynCom.ID <- factor(sample_data(JH201_tube_rare_live)$SynCom.ID, levels=c("RKMP", "RKMP.NoBi27"))
#Genotype
sample_data(JH201_tube_rare_live)$Genotype <- factor(sample_data(JH201_tube_rare_live)$Genotype, levels=c("Morex", "Barke", "SL17", "SL52", "Bulk"))

#plotting
p <- plot_bar(JH201_tube_rare_live, "Genus", fill="Genus", facet_grid=SynCom.ID~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subest for "HI" samples
JH201_tube_rare_dead <- subset_samples(JH201_tube_rare, Syncom.Status == "Dead")

#re-arrange samples
#SynCom
sample_data(JH201_tube_rare_dead)$SynCom.ID <- factor(sample_data(JH201_tube_rare_dead)$SynCom.ID, levels=c("RKMP", "RKMP.NoBi27"))
#Genotype
sample_data(JH201_tube_rare_dead)$Genotype <- factor(sample_data(JH201_tube_rare_dead)$Genotype, levels=c("Bulk", "Morex", "Barke", "SL17", "SL52"))

#plotting
p <- plot_bar(JH201_tube_rare_dead, "Genus", fill="Genus", facet_grid=SynCom.ID~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subest for controls
JH201_tube_rare_control <- subset_samples(JH201_tube_rare, Syncom.Status == "N.A")

#re-arrange samples
#Genotype
sample_data(JH201_tube_rare_control)$Genotype <- factor(sample_data(JH201_tube_rare_control)$Genotype, levels=c("Bulk", "Morex", "Barke", "SL17", "SL52"))

#plotting
p <- plot_bar(JH201_tube_rare_control, "Genus", fill="Genus", facet_grid=SynCom.ID~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subset for run with RKMP and all genotypes
JH201_tube_rare_r3 <- subset_samples(JH201_tube_rare, Date.Inoc == "230406")
sample_data(JH201_tube_rare_r3)
sample_data(JH201_tube_rare_r3)$Genotype <- factor(sample_data(JH201_tube_rare_r3)$Genotype, levels=c("Morex", "Barke", "SL17", "SL52", "Bulk"))
p <- plot_bar(JH201_tube_rare_r3, "Genus", fill="Genus", facet_grid=~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subset for run with RKMP_NoBi27 and all genotypes
JH201_tube_rare_r2 <- subset_samples(JH201_tube_rare, Date.Inoc == "230223")
sample_data(JH201_tube_rare_r2)
sample_data(JH201_tube_rare_r2)$Genotype <- factor(sample_data(JH201_tube_rare_r2)$Genotype, levels=c("Morex", "Barke", "SL17", "SL52", "Bulk"))
p <- plot_bar(JH201_tube_rare_r2, "Genus", fill="Genus", facet_grid=~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subset for run with Barke and both SynComs - Live only
JH201_tube_rare_r4 <- subset_samples(JH201_tube_rare, Date.Inoc == "230525")
JH201_tube_rare_r4 <- subset_samples(JH201_tube_rare_r4, Syncom.Status == "Live")
sample_data(JH201_tube_rare_r4)
sample_data(JH201_tube_rare_r4)$Genotype <- factor(sample_data(JH201_tube_rare_r4)$Genotype, levels=c("Barke", "Bulk"))
p <- plot_bar(JH201_tube_rare_r4, "Genus", fill="Genus", facet_grid=SynCom.ID~Genotype)
p + scale_fill_viridis(discrete=TRUE)

#subset for Barke live across all runs
JH201_tube_rare_barke <- subset_samples(JH201_tube_rare, Genotype == "Barke")
JH201_tube_rare_barke <- subset_samples(JH201_tube_rare_barke, Syncom.Status == "Live")
sample_data(JH201_tube_rare_barke)
p <- plot_bar(JH201_tube_rare_barke, "Genus", fill="Genus", facet_grid=SynCom.ID~Date.Inoc)
p + scale_fill_viridis(discrete=TRUE)

####THESIS####
#subset for Barke RKMP live across all runs
JH201_tube_rare_barke_RKMP <- subset_samples(JH201_tube_rare_barke, SynCom.ID == "RKMP")
sample_data(JH201_tube_rare_barke_RKMP)
RunLabels <- c("230406"="Experiment 1", "230525"="Experiment 2", "230831"="Experiment 3")
p <- plot_bar(JH201_tube_rare_barke_RKMP, "Genus", fill="Genus", facet_grid=)
p + scale_fill_viridis(discrete=TRUE)+facet_grid(SynCom.ID~Date.Inoc, labeller = labeller(Date.Inoc = RunLabels))+
  theme(legend.position = "right", axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

ID_labels <- c("RKMP"="RKMP", "RKMP.NoBi27"="RKMP-27")
status_labels <- c("Live"="Live","Dead"="HI")
p <- plot_bar(JH201_inoculants2_rare_noR5, "Genus", fill="Genus")
p + scale_fill_viridis(discrete=TRUE)

#subset for Barke and Bulk RKMP live across all runs
JH201_tube_rare_bb_RKMP <- subset_samples(JH201_tube_rare, Genotype != "Morex")
JH201_tube_rare_bb_RKMP <- subset_samples(JH201_tube_rare_bb_RKMP, Genotype != "SL17")
JH201_tube_rare_bb_RKMP <- subset_samples(JH201_tube_rare_bb_RKMP, Genotype != "SL52")
JH201_tube_rare_bb_RKMP <- subset_samples(JH201_tube_rare_bb_RKMP, SynCom.ID == "RKMP")
JH201_tube_rare_bb_RKMP <- subset_samples(JH201_tube_rare_bb_RKMP, Syncom.Status != "Dead")
sample_data(JH201_tube_rare_bb_RKMP)
p <- plot_bar(JH201_tube_rare_bb_RKMP, "Genus", fill="Genus", facet_grid=Genotype~Date.Inoc)
p + scale_fill_viridis(discrete=TRUE)
