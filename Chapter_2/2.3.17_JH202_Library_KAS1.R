#####################################################################################
#BUILT-IN CODE
####################################################################################
#############################################################
#
# Ref to the Thesis
# 
#  Code to compute calculations SynComs KAS1
#  Revision 05/25
#  c.arnton@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk 

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################
#required packages 
library("phyloseq")
library("vegan")
library ("ggplot2")
library ("DESeq2")
library("viridis")

#############################################################
#set working directory: DB cpu
#setwd("/cluster/db/R_shared/NSQ04_v0525/")
#set working directory: Katie home
setwd("C:/Users/catar/Documents/Dundee PhD/SynCom/JH202")
#set working dirctory: Katie work
#setwd("~/PhD Files/SynCom/JH202")

#skip to line 164 for pre-processed data 

#import the datasets
NSQ04 <-readRDS("NSQ04_dada2_silva_138.1.rds")

#inspect the files
NSQ04

#total number of reads
sum(sample_sums(NSQ04))

#and its three "constituents"

#ASV counts
otu_table(NSQ04)

#Taxonomy information
tax_table(NSQ04)

#mapping files
sample_data(NSQ04)

#stats overall library
sum(sample_sums(NSQ04))

#proceed with independent filtering in both objects prior merging
##################################################################
#Pre-processing: remove Chloroplast and Mitochondria but retain NA
#rationale: 16S rRNA primers may amplify plant-derived sequences 
#those may interfere with data analysis (we are after bacteria, not plants)
#################################################################

#NSQ04
NSQ04_no_chlor <-subset_taxa(NSQ04, (Order!="Chloroplast") | is.na(Order))
NSQ04_no_chlor

NSQ04_no_plants <-subset_taxa(NSQ04_no_chlor, (Family!="Mitochondria") | is.na(Family))
NSQ04_no_plants

##################################################################
#Prune putative contaminant ASVs
#Rationale: bacteria are everywhere and 16S rRNA amplification is very efficient
#We use a list of bacteria likely representing the contamination we may face in our lab as baseline
#More info here: https://www.frontiersin.org/articles/10.3389/fmicb.2018.01650/full 
##################################################################

#Import the list of contaminat ASVs from JH06 library
Contaminant_ASVs <- read.delim("JH06_contaminant_ASVs_ids.txt", header = FALSE)

#identify the proportion of putative contaminants in the merged object
NSQ04_contaminants <- intersect(taxa_names(NSQ04_no_plants),Contaminant_ASVs)
NSQ04_contaminants

#########################################################################
# Remove ASVs assigned to NA at phylum level
# Rationale: if a sequence cannot be assigned at a very high level such as Phylum is not of use us
#########################################################################

NSQ04_no_plants_1 <- subset_taxa(NSQ04_no_plants, Phylum!= "NA")
NSQ04_no_plants_1

#########################################################################
#Subsest for KAS1
#remove editing experiment
#remove ASV with 0 counts
#########################################################################

NSQ04_KAS1 <- subset_samples(NSQ04_no_plants_1, SynCom.ID == "KAS1")
NSQ04_KAS1 <- subset_samples(NSQ04_KAS1, Date.Harvest != "241129")
NSQ04_KAS1_integer <- prune_taxa(taxa_sums(NSQ04_KAS1) > 0, NSQ04_KAS1)
NSQ04_KAS1_integer

#########################################################################
# inspect reads distribution across objects
# Rationale: if we want to compare like-with-like we should know sample distribution
#########################################################################

#summary
sum(sample_sums(NSQ04_KAS1))
min(sample_sums(NSQ04_KAS1))
max(sample_sums(NSQ04_KAS1))
mean(sample_sums(NSQ04_KAS1))

#graphical outputs
sort(sample_sums(NSQ04_KAS1))
hist(sample_sums(NSQ04_KAS1), main = paste("NSQ04 reads distribution; KAS1"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: Heat Killed
NSQ04_KAS1_HK <- subset_samples(NSQ04_KAS1, Syncom.Status == "HK")
sort(sample_sums(NSQ04_KAS1_HK))
hist(sample_sums(NSQ04_KAS1_HK), main = paste("NSQ04 reads distribution; KAS1_HK"), xlab = paste("reads"), ylab = paste ("number of samples"))

#inspect data distribution of various experiments: Live
NSQ04_KAS1_Live <- subset_samples(NSQ04_KAS1, Syncom.Status == "Live")
sort(sample_sums(NSQ04_KAS1_Live))
hist(sample_sums(NSQ04_KAS1_Live), main = paste("NSQ04 reads distribution; KAS1_Live"), xlab = paste("reads"), ylab = paste ("number of samples"))

#########################################################################
# work with live samples and set criteria to work with unevenly sequenced samples
# a) remove samples with less than 1,000 reads
# b) keep ASV with more than 15 reads in at least 10% of the samples
# c) rarefy at even sequencing depth
# d) agglomerate at genus level
#########################################################################

#remove poorly sequenced samples
NSQ04_KAS1_Live_30K <- prune_samples(sample_sums(NSQ04_KAS1_Live) > 30000, NSQ04_KAS1_Live)

#set an arbitrary threshold of 10 reads in 10% of samples
NSQ04_KAS1_Live_30K_threshold = filter_taxa(NSQ04_KAS1_Live_30K, function(x) sum(x > 15) > (0.1 *length(x)), TRUE)
NSQ04_KAS1_Live_30K_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_KAS1_Live_30K_threshold))

#agglomerate at genus level
NSQ04_KAS1_Live_genus <- tax_glom(NSQ04_KAS1_Live_30K_threshold, taxrank="Genus")

#inspect the files
NSQ04_KAS1_Live_genus
sort(sample_sums(NSQ04_KAS1_Live_genus))

#rarefy at 29K
NSQ04_KAS1_Live_threshold_29K <- rarefy_even_depth(NSQ04_KAS1_Live_genus, 29000)

#save the file for the reproducibility of the code
#saveRDS(NSQ04_KAS1_Live_threshold_29K, file = "NSQ04_KAS1_Live_calculation.rds")
NSQ04_KAS1_calculation <-readRDS("NSQ04_KAS1_Live_calculation.rds")
NSQ04_KAS1_calculation

#########################################################################
#Rhizo6 poster ordination with genotype & replicate effect
#########################################################################


#data visualisation
#re-order the factors
sample_data(NSQ04_KAS1_calculation)$Genotype <- factor(sample_data(NSQ04_KAS1_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_KAS1_calculation)$Date.Inoc <- factor(sample_data(NSQ04_KAS1_calculation)$Date.Inoc, levels=c("241031", "241107"))

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_KAS1_calculation, x="Genotype", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove Observed, this is just a confirmation
p = plot_richness(NSQ04_KAS1_calculation, x="Genotype", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove inoculum for plotting and statistical analysis
NSQ04_KAS1_calculation_ni <- subset_samples(NSQ04_KAS1_calculation, Genotype != "Inoculum")

#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#without inoculum
p = plot_richness(NSQ04_KAS1_calculation_ni, x="Genotype", color="Genotype", measures=c("Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#ordination
#https://strata.uga.edu/software/pdf/mdsTutorial.pdf
NSQ04_KAS1.ord <- ordinate(NSQ04_KAS1_calculation, "NMDS", "bray")
p = plot_ordination(NSQ04_KAS1_calculation, NSQ04_KAS1.ord, type="samples", color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

####PAPER####
#Constrain for genotype effect
NSQ04_KAS1.cap <- ordinate(NSQ04_KAS1_calculation, "CAP", "bray", ~ Genotype)
#note the formula to specify what factor(s) look for
a = plot_ordination(NSQ04_KAS1_calculation, NSQ04_KAS1.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "Replicate", labels = c("1", "2")) + theme(legend.position="bottom")+ guides(colour = "none")
a

#Stacked bar
p = plot_bar(NSQ04_KAS1_calculation, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#no inoculum
p = plot_bar(NSQ04_KAS1_calculation_ni, "Genus", facet_grid=Date.Inoc ~ Genotype, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Only inoculum
NSQ04_KAS1_calculation_i <- subset_samples(NSQ04_KAS1_calculation, Genotype == "Inoculum")
p = plot_bar(NSQ04_KAS1_calculation_i, "Genus", facet_grid=Genotype ~ Date.Inoc, fill='Genus')
p = p + scale_fill_viridis(discrete=TRUE)
p

#Constrain genotype effect - no inoculum
NSQ04_KAS1_ni.cap <- ordinate(NSQ04_KAS1_calculation_ni, "CAP", "bray", ~ Genotype)
#note the formula to specify what factor(s) look for
b = plot_ordination(NSQ04_KAS1_calculation_ni, NSQ04_KAS1_ni.cap, type="samples", color="Genotype", shape="Date.Inoc") 
b = b + geom_point(size = 4, alpha = 0.75)
b = b + scale_shape_manual(values = c(16, 17))
b = b + scale_colour_manual(values = SynCom_col_ni)
b

#Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_calculation_ni, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Plant * Date.Inoc, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_ni))), permutations = 5000)
Stat

#genotype effect
#remove unplanted for plotting and statistical analysis
NSQ04_KAS1_calculation_plant <- subset_samples(NSQ04_KAS1_calculation_ni, Genotype != "Bulk")

#set the colours: no inoculum; no bulk
SynCom_col_plant <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Constrain for genotype effect
NSQ04_KAS1_plant.cap <- ordinate(NSQ04_KAS1_calculation_plant, "CAP", "bray", ~ Genotype)
#Plot
c = plot_ordination(NSQ04_KAS1_calculation_plant, NSQ04_KAS1_plant.cap, type="samples", color="Genotype", shape="Date.Inoc") 
c = c + geom_point(size = 4, alpha = 0.75)
c = c + scale_shape_manual(values = c(16, 17))
c = c + scale_colour_manual(values = SynCom_col_plant)+
  scale_shape_discrete(name = "Assembly Method", labels = c("By OD", "By CFU")) + theme(legend.position="bottom")
c

#split the runs
c = plot_ordination(NSQ04_KAS1_calculation_plant, NSQ04_KAS1_plant.cap, type="samples", color="Genotype", shape="Date.Inoc") 
c = c + geom_point(size = 4, alpha = 0.75)
c = c + scale_shape_manual(values = c(16, 17))
c = c + scale_colour_manual(values = SynCom_col_plant)
c = c + facet_wrap(~Date.Inoc, 1)
c

#Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_calculation_plant, "bray")
BC_bacteria 

#here we can use a formula ANOVA-like to identify factors of interest
Stat <- adonis2(BC_bacteria ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_plant))), permutations = 5000)
Stat

#faceted graph of all three
grid.arrange(a, b, c, ncol=3)

#Run 3 planted only
NSQ04_KAS1_calculation_plant_r3 <- subset_samples(NSQ04_KAS1_calculation_plant, Date.Inoc == "241031")
NSQ04_KAS1_plant_r3.cap <- ordinate(NSQ04_KAS1_calculation_plant_r3, "CAP", "bray", ~ Genotype)
c = plot_ordination(NSQ04_KAS1_calculation_plant_r3, NSQ04_KAS1_plant_r3.cap, type="samples", color="Genotype") 
c = c + geom_point(size = 4, alpha = 0.75)
c = c + scale_shape_manual(values = c(16))
c = c + scale_colour_manual(values = SynCom_col_plant)
c
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_calculation_plant_r3, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_plant_r3))), permutations = 5000)
Stat

#Run 4 planted only
NSQ04_KAS1_calculation_plant_r4 <- subset_samples(NSQ04_KAS1_calculation_plant, Date.Inoc == "241107")
NSQ04_KAS1_plant_r4.cap <- ordinate(NSQ04_KAS1_calculation_plant_r4, "CAP", "bray", ~ Genotype)
c = plot_ordination(NSQ04_KAS1_calculation_plant_r4, NSQ04_KAS1_plant_r4.cap, type="samples", color="Genotype") 
c = c + geom_point(size = 4, alpha = 0.75)
c = c + scale_shape_manual(values = c(17))
c = c + scale_colour_manual(values = SynCom_col_plant)
c
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_calculation_plant_r4, "bray")
BC_bacteria 
Stat <- adonis2(BC_bacteria ~ Genotype, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_plant_r4))), permutations = 5000)
Stat

#run DESeq
#extract count data 
NSQ04_genotype_plant <- otu_table(NSQ04_KAS1_calculation_plant)
countData = as.data.frame(NSQ04_genotype_plant)
colnames(NSQ04_genotype_plant)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_plant)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_plant_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_plant_test <- DESeq(NSQ04_plant_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in plant compartments samples
rhizo_genotype <- results(NSQ04_plant_test, contrast = c("Genotype", "Barke", "Morex")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [inf any]
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_FDR005

#17
rhizo_genotype <- results(NSQ04_plant_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

#52
rhizo_genotype <- results(NSQ04_plant_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

#extract taxonomic information

NSQ04_Morex <- prune_taxa(rownames(rhizo_genotype_FDR005), NSQ04_KAS1_calculation)
NSQ04_Morex

#who is there?
tax_table(NSQ04_Morex)

#data visualisation
#% of enriched Chryseo vs. total
Chryseo_enriched <- as.data.frame(sample_sums(NSQ04_Morex)/sample_sums(NSQ04_KAS1_calculation))
colnames(Chryseo_enriched) <- c("Chryseobacterium_proportion")

#extract the mapping file
Chryseo_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))

#merge the dataset
Chryseo_plotting <- cbind(Chryseo_map, Chryseo_enriched)

#Order the factors
Chryseo_plotting$Genotype <- ordered(Chryseo_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))

#set the colours: with inoculum, white bulk
SynCom_col_wb <- c("#999999", "#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Plotting Chryseo
#dev.off()
p <- ggplot(Chryseo_plotting, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75))

#Without inoculum
NSQ04_Morex_ni <- prune_taxa(rownames(rhizo_genotype_FDR005), NSQ04_KAS1_calculation_ni)
Chryseo_enriched_ni <- as.data.frame(sample_sums(NSQ04_Morex_ni)/sample_sums(NSQ04_KAS1_calculation_ni))
colnames(Chryseo_enriched_ni) <- c("Chryseobacterium_proportion")

#extract the mapping file
Chryseo_map_ni <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation_ni)))

#merge the dataset
Chryseo_plotting_ni <- cbind(Chryseo_map_ni, Chryseo_enriched_ni)

#Order the factors
Chryseo_plotting_ni$Genotype <- ordered(Chryseo_plotting_ni$Genotype, levels=c("Bulk", "Barke", "124-17", "124-52", "Morex"))

#set the colours: no inoculum, white bulk
SynCom_col_ni_wb <- c("#FFFFFF", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#Plotting Chryseo
#dev.off()
r3 <- ggplot(Chryseo_plotting_ni, aes(x=Genotype, y=Chryseobacterium_proportion, fill = Genotype)) 
r3 <- r3 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Chryseobacterium_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r3

#Now for Bulk vs Genotypes

#extract count data 
NSQ04_genotype <- otu_table(NSQ04_KAS1_calculation)
countData = as.data.frame(NSQ04_genotype)
colnames(NSQ04_genotype)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_test <- DESeq(NSQ04_cds, fitType="local", betaPrior = FALSE) 

#Barke
rhizo_genotype <- results(NSQ04_test, contrast = c("Genotype", "Bulk", "Barke")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

#Morex
rhizo_genotype <- results(NSQ04_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

#17
rhizo_genotype <- results(NSQ04_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

#52
rhizo_genotype <- results(NSQ04_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_FDR005 <- rhizo_genotype[(rownames(rhizo_genotype)[which(rhizo_genotype$padj <0.05)]), ]
rhizo_genotype_FDR005

####PAPER####
#Try to make the proportion charts for everything

#define the Genera differentially enriched in plant compartments samples
rhizo_genotype <- results(NSQ04_plant_test, contrast = c("Genotype", "Barke", "Morex"))

#extract  first row
rhizo_genotype_r1 <- rhizo_genotype[1, ]
#inspect the generated file
rhizo_genotype_r1
#Link back to full dataset
KAS1_r1 <- prune_taxa(rownames(rhizo_genotype_r1), NSQ04_KAS1_calculation)
KAS1_r1
#who is there?
tax_table(KAS1_r1)
#Pseudomonas
#% of enriched vs. total
r1_enriched <- as.data.frame(sample_sums(KAS1_r1)/sample_sums(NSQ04_KAS1_calculation))
colnames(r1_enriched) <- c("Pseudomonas_proportion")
#extract the mapping file
r1_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r1_plotting <- cbind(r1_map, r1_enriched)
#Order the factors
r1_plotting$Genotype <- ordered(r1_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r1_plotting, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r1_plotting_ni <- subset(r1_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r1 <- ggplot(r1_plotting_ni, aes(x=Genotype, y=Pseudomonas_proportion, fill = Genotype)) 
r1 <- r1 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pseudomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r1

#extract row 2
rhizo_genotype_r2 <- rhizo_genotype[2, ]
#inspect the generated file
rhizo_genotype_r2
#Link back to full dataset
KAS1_r2 <- prune_taxa(rownames(rhizo_genotype_r2), NSQ04_KAS1_calculation)
KAS1_r2
#who is there?
tax_table(KAS1_r2)
#Pedobacter
#% of enriched vs. total
r2_enriched <- as.data.frame(sample_sums(KAS1_r2)/sample_sums(NSQ04_KAS1_calculation))
colnames(r2_enriched) <- c("Pedobacter_proportion")
#extract the mapping file
r2_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r2_plotting <- cbind(r2_map, r2_enriched)
#Order the factors
r2_plotting$Genotype <- ordered(r2_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r2_plotting, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r2_plotting_ni <- subset(r2_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r2 <- ggplot(r2_plotting_ni, aes(x=Genotype, y=Pedobacter_proportion, fill = Genotype)) 
r2 <- r2 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pedobacter_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r2

#extract row 3
rhizo_genotype_r3 <- rhizo_genotype[3, ]
#inspect the generated file
rhizo_genotype_r3
#Link back to full dataset
KAS1_r3 <- prune_taxa(rownames(rhizo_genotype_r3), NSQ04_KAS1_calculation)
KAS1_r3
#who is there?
tax_table(KAS1_r3)
#Chryseobacterium - can skip as already have this one

#extract row 4
rhizo_genotype_r4 <- rhizo_genotype[4, ]
#inspect the generated file
rhizo_genotype_r4
#Link back to full dataset
KAS1_r4 <- prune_taxa(rownames(rhizo_genotype_r4), NSQ04_KAS1_calculation)
KAS1_r4
#who is there?
tax_table(KAS1_r4)
#Stenotrophomonas
#% of enriched vs. total
r4_enriched <- as.data.frame(sample_sums(KAS1_r4)/sample_sums(NSQ04_KAS1_calculation))
colnames(r4_enriched) <- c("Stenotrophomonas_proportion")
#extract the mapping file
r4_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r4_plotting <- cbind(r4_map, r4_enriched)
#Order the factors
r4_plotting$Genotype <- ordered(r4_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r4_plotting, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r4_plotting_ni <- subset(r4_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r4 <- ggplot(r4_plotting_ni, aes(x=Genotype, y=Stenotrophomonas_proportion, fill = Genotype)) 
r4 <- r4 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Stenotrophomonas_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r4

#extract row 5
rhizo_genotype_r5 <- rhizo_genotype[5, ]
#inspect the generated file
rhizo_genotype_r5
#Link back to full dataset
KAS1_r5 <- prune_taxa(rownames(rhizo_genotype_r5), NSQ04_KAS1_calculation)
KAS1_r5
#who is there?
tax_table(KAS1_r5)
#Bacillus
#% of enriched vs. total
r5_enriched <- as.data.frame(sample_sums(KAS1_r5)/sample_sums(NSQ04_KAS1_calculation))
colnames(r5_enriched) <- c("Bacillus_proportion")
#extract the mapping file
r5_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r5_plotting <- cbind(r5_map, r5_enriched)
#Order the factors
r5_plotting$Genotype <- ordered(r5_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r5_plotting, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r5_plotting_ni <- subset(r5_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r5 <- ggplot(r5_plotting_ni, aes(x=Genotype, y=Bacillus_proportion, fill = Genotype)) 
r5 <- r5 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Bacillus_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r5

#extract row 6
rhizo_genotype_r6 <- rhizo_genotype[6, ]
#inspect the generated file
rhizo_genotype_r6
#Link back to full dataset
KAS1_r6 <- prune_taxa(rownames(rhizo_genotype_r6), NSQ04_KAS1_calculation)
KAS1_r6
#who is there?
tax_table(KAS1_r6)
#Pseudarthrobacter
#% of enriched vs. total
r6_enriched <- as.data.frame(sample_sums(KAS1_r6)/sample_sums(NSQ04_KAS1_calculation))
colnames(r6_enriched) <- c("Pseudarthrobacter_proportion")
#extract the mapping file
r6_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r6_plotting <- cbind(r6_map, r6_enriched)
#Order the factors
r6_plotting$Genotype <- ordered(r6_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r6_plotting, aes(x=Genotype, y=Pseudarthrobacter_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Pseudarthrobacter_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r6_plotting_ni <- subset(r6_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r6 <- ggplot(r6_plotting_ni, aes(x=Genotype, y=Pseudarthrobacter_proportion, fill = Genotype)) 
r6 <- r6 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Pseudarthrobacter_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r6

#extract row 7
rhizo_genotype_r7 <- rhizo_genotype[7, ]
#inspect the generated file
rhizo_genotype_r7
#Link back to full dataset
KAS1_r7 <- prune_taxa(rownames(rhizo_genotype_r7), NSQ04_KAS1_calculation)
KAS1_r7
#who is there?
tax_table(KAS1_r7)
#Rhodococcus
#% of enriched vs. total
r7_enriched <- as.data.frame(sample_sums(KAS1_r7)/sample_sums(NSQ04_KAS1_calculation))
colnames(r7_enriched) <- c("Rhodococcus_proportion")
#extract the mapping file
r7_map <- as.data.frame(as.matrix(sample_data(NSQ04_KAS1_calculation)))
#merge the dataset
r7_plotting <- cbind(r7_map, r7_enriched)
#Order the factors
r7_plotting$Genotype <- ordered(r7_plotting$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
#Plotting
#dev.off()
p <- ggplot(r7_plotting, aes(x=Genotype, y=Rhodococcus_proportion, fill = Genotype)) 
p <- p + geom_boxplot() + scale_fill_manual(values = SynCom_col_wb)
p + geom_point(aes(y=Rhodococcus_proportion, group=Genotype), position = position_dodge(width=0.75))
#Without inoculum
r7_plotting_ni <- subset(r7_plotting, Genotype!="Inoculum")
#Plot
#dev.off()
r7 <- ggplot(r7_plotting_ni, aes(x=Genotype, y=Rhodococcus_proportion, fill = Genotype)) 
r7 <- r7 + geom_boxplot() + scale_fill_manual(values = SynCom_col_ni_wb) + 
  geom_point(aes(y=Rhodococcus_proportion, group=Genotype), position = position_dodge(width=0.75))+
  ylim(0,0.5) + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
r7

#faceted graph of all
grid.arrange(r1, r2, r3, r4, r5, r6, r7, ncol=4)

#########################################################################
# switch to HK and set criteria to work with unevenly sequenced samples
# a) remove samples with low reads
# b) keep ASV with more than 10 reads in at least 10% of the samples
# c) rarefy at even sequencing depth
# d) agglomerate at genus level
#Skip and load file in line 501
#########################################################################

#remove poorly sequenced samples
NSQ04_KAS1_HK_2.5K <- prune_samples(sample_sums(NSQ04_KAS1_HK) > 2500, NSQ04_KAS1_HK)

#set an arbitrary threshold of 10 reads in 10% of samples
NSQ04_KAS1_HK_2.5K_threshold = filter_taxa(NSQ04_KAS1_HK_2.5K, function(x) sum(x > 10) > (0.1 *length(x)), TRUE)
NSQ04_KAS1_HK_2.5K_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_KAS1_HK_2.5K_threshold))

#agglomerate at genus level
NSQ04_KAS1_HK_genus <- tax_glom(NSQ04_KAS1_HK_2.5K_threshold, taxrank="Genus")

#inspect the files
NSQ04_KAS1_HK_genus
sort(sample_sums(NSQ04_KAS1_HK_genus))
#Some of these have lost a massive amount of reads in the agglomeration, spread is massive

#rarefy at 500
NSQ04_KAS1_HK_threshold_500 <- rarefy_even_depth(NSQ04_KAS1_HK_genus, 500)

#data visualisation
p = plot_bar(NSQ04_KAS1_HK_threshold_500, "Genus", facet_grid=Date.Inoc ~ Genotype, fill="Genus") 
p = p + scale_fill_viridis(discrete=TRUE)
p

#proportion of retained reads
ratio <- sum(sample_sums(NSQ04_KAS1_HK_threshold_500))/sum(sample_sums(NSQ04_KAS1_HK_genus))*100
ratio

#Try again getting rid of everything under 10k - a good chunk of samples
NSQ04_KAS1_HK_10K <- prune_samples(sample_sums(NSQ04_KAS1_HK) > 10000, NSQ04_KAS1_HK)

#set an arbitrary threshold of 10 reads in 10% of samples
NSQ04_KAS1_HK_10K_threshold = filter_taxa(NSQ04_KAS1_HK_10K, function(x) sum(x > 10) > (0.1 *length(x)), TRUE)
NSQ04_KAS1_HK_10K_threshold

#inspect read distribution upon thresholding
sort(sample_sums(NSQ04_KAS1_HK_10K_threshold))

#agglomerate at genus level
NSQ04_KAS1_HK_genus_10k <- tax_glom(NSQ04_KAS1_HK_10K_threshold, taxrank="Genus")

#inspect the files
NSQ04_KAS1_HK_genus_10k
sort(sample_sums(NSQ04_KAS1_HK_genus_10k))

#rarefy at 500
NSQ04_KAS1_HK_10k_threshold_500 <- rarefy_even_depth(NSQ04_KAS1_HK_genus_10k, 500)

#data visualisation
p = plot_bar(NSQ04_KAS1_HK_10k_threshold_500, "Genus", facet_grid=Date.Inoc ~ Genotype, fill="Genus")
p = p + scale_fill_viridis(discrete=TRUE)
p

#proportion of retained reads
ratio <- sum(sample_sums(NSQ04_KAS1_HK_10k_threshold_500))/sum(sample_sums(NSQ04_KAS1_HK_genus_10k))*100
ratio

#Take the original 2.5k and prune again prior to rarefication
#inspect the files
NSQ04_KAS1_HK_genus
sort(sample_sums(NSQ04_KAS1_HK_genus))

#Remove samples below 2k
NSQ04_KAS1_HK_genus_2K <- prune_samples(sample_sums(NSQ04_KAS1_HK_genus) > 2000, NSQ04_KAS1_HK_genus)

#rarefy at 2100
NSQ04_KAS1_HK_threshold_2100 <- rarefy_even_depth(NSQ04_KAS1_HK_genus_2K, 2100)

#data visualisation
p = plot_bar(NSQ04_KAS1_HK_threshold_2100, "Genus", facet_grid=Date.Inoc ~ Genotype, fill="Genus")
p = p + scale_fill_viridis(discrete=TRUE)
p

#proportion of retained reads
ratio <- sum(sample_sums(NSQ04_KAS1_HK_threshold_2100))/sum(sample_sums(NSQ04_KAS1_HK_genus_2K))*100
ratio

#continuing forward with the first version
#save the file for the reproducibility of the code
#saveRDS(NSQ04_KAS1_HK_threshold_500, file = "NSQ04_KAS1_HK_calculation.rds")
NSQ04_KAS1_HK_calculation <-readRDS("NSQ04_KAS1_HK_calculation.rds")
NSQ04_KAS1_HK_calculation

####HK Data Graphs and Analysis####

#re-order the factors
sample_data(NSQ04_KAS1_HK_calculation)$Genotype <- factor(sample_data(NSQ04_KAS1_HK_calculation)$Genotype, levels=c("Inoculum", "Bulk", "Barke", "124-17", "124-52", "Morex"))
sample_data(NSQ04_KAS1_HK_calculation)$Date.Inoc <- factor(sample_data(NSQ04_KAS1_HK_calculation)$Date.Inoc, levels=c("241031", "241107"))

#Stacked bar plot - with factors in order
p = plot_bar(NSQ04_KAS1_HK_calculation, "Genus", facet_grid=Date.Inoc ~ Genotype, fill="Genus")
p = p + scale_fill_viridis(discrete=TRUE)
p

#set the colours: with inoculum
SynCom_col <- c("#999999", "#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#alpha diversity calculation
#https://www.nature.com/articles/s41598-024-77864-y
p = plot_richness(NSQ04_KAS1_HK_calculation, x="Genotype", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)
p

#remove inoculum for plotting and statistical analysis
NSQ04_KAS1_HK_calculation_ni <- subset_samples(NSQ04_KAS1_HK_calculation, Genotype != "Inoculum")

#set the colours: no inoculum
SynCom_col_ni <- c("#000000", "#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#without inoculum
p = plot_richness(NSQ04_KAS1_HK_calculation_ni, x="Genotype", color="Genotype", measures=c("Observed", "Shannon"))
p + geom_point(size=5, alpha=0.7)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#ordination

####PAPER####
#Unconstrained
NSQ04_KAS1_HK.ord <- ordinate(NSQ04_KAS1_HK_calculation, "NMDS", "bray")
#plot
p = plot_ordination(NSQ04_KAS1_HK_calculation, NSQ04_KAS1_HK.ord, type="samples", color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col)+
  scale_shape_discrete(name = "Replicate", labels = c("1", "2")) + theme(legend.position="bottom")+ guides(colour = "none")
p

#Constrain for genotype effect
NSQ04_KAS1_HK_G.cap <- ordinate(NSQ04_KAS1_HK_calculation, "CAP", "bray", ~ Genotype)
#note the formula to specify what factor(s) look for
#Plot
a = plot_ordination(NSQ04_KAS1_HK_calculation, NSQ04_KAS1_HK_G.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Constrain for run effect
NSQ04_KAS1_HK_D.cap <- ordinate(NSQ04_KAS1_HK_calculation, "CAP", "bray", ~ Date.Inoc)
#Plot
a = plot_ordination(NSQ04_KAS1_HK_calculation, NSQ04_KAS1_HK_D.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col)
a

#Repeat ordination without inoculum
NSQ04_KAS1_HK_ni.ord <- ordinate(NSQ04_KAS1_HK_calculation_ni, "NMDS", "bray")

#and plot
p = plot_ordination(NSQ04_KAS1_HK_calculation_ni, NSQ04_KAS1_HK_ni.ord, type="samples", color="Genotype", shape="Date.Inoc") 
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_shape_manual(values = c(16, 17))
p = p + scale_colour_manual(values = SynCom_col_ni)
p

#alternative plotting where we "constrain for" genotype effect
NSQ04_KAS1_HK_ni_G.cap <- ordinate(NSQ04_KAS1_HK_calculation_ni, "CAP", "bray", ~ Genotype)

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_KAS1_HK_calculation_ni, NSQ04_KAS1_HK_ni_G.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_ni)
a

#alternative plotting where we "constrain for" run effect
NSQ04_KAS1_HK_ni_D.cap <- ordinate(NSQ04_KAS1_HK_calculation_ni, "CAP", "bray", ~ Date.Inoc)

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_KAS1_HK_calculation_ni, NSQ04_KAS1_HK_ni_D.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_ni)
a

#Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_HK_calculation_ni, "bray")
BC_bacteria 

#here we can use a formula ANOVA-like to identify factors of interest
Stat <- adonis2(BC_bacteria ~ Plant * Date.Inoc, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_HK_calculation_ni))), permutations = 5000)
Stat

#remove unplanted for plotting and statistical analysis
NSQ04_KAS1_HK_calculation_plant <- subset_samples(NSQ04_KAS1_HK_calculation_ni, Genotype != "Bulk")

#repeat ordination with planted only
NSQ04_KAS1_HK_plant.ord <- ordinate(NSQ04_KAS1_HK_calculation_plant, "NMDS", "bray")

#set the colours: no inoculum; no bulk
SynCom_col_plant <- c("#0072B2", "#56B4E9", "#E69F00", "#CC79A7")

#ggplots function to increase effectivness of the visualisation
c = plot_ordination(NSQ04_KAS1_HK_calculation_plant, NSQ04_KAS1_HK_plant.ord, type="samples", color="Genotype", shape="Date.Inoc") 
c = c + geom_point(size = 4, alpha = 0.75)
c = c + scale_shape_manual(values = c(16, 17))
c = c + scale_colour_manual(values = SynCom_col_plant)
c

#Constrained

#alternative plotting where we "constrain for" genotype effect
NSQ04_KAS1_HK_plant_G.cap <- ordinate(NSQ04_KAS1_HK_calculation_plant, "CAP", "bray", ~ Genotype)
#note the formula to specify what factor(s) look for

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_KAS1_HK_calculation_plant, NSQ04_KAS1_HK_plant_G.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_plant)
a

#alternative plotting where we "constrain for" run effect
NSQ04_KAS1_HK_plant_D.cap <- ordinate(NSQ04_KAS1_HK_calculation_plant, "CAP", "bray", ~ Date.Inoc)

#ggplots function to increase effectivness of the visualisation
a = plot_ordination(NSQ04_KAS1_HK_calculation_plant, NSQ04_KAS1_HK_plant_D.cap, type="samples", color="Genotype", shape="Date.Inoc") 
a = a + geom_point(size = 4, alpha = 0.75)
a = a + scale_shape_manual(values = c(16, 17))
a = a + scale_colour_manual(values = SynCom_col_plant)
a

#Permanova
BC_bacteria  <- phyloseq::distance(NSQ04_KAS1_HK_calculation_plant, "bray")
BC_bacteria 

#here we can use a formula ANOVA-like to identify factors of interest
Stat <- adonis2(BC_bacteria ~ Genotype * Date.Inoc, data= as.data.frame(as.matrix(sample_data(NSQ04_KAS1_HK_calculation_plant))), permutations = 5000)
Stat

#Deseq
#extract count data 
NSQ04_HK_genotype_ni <- otu_table(NSQ04_KAS1_HK_calculation_ni)
countData = as.data.frame(NSQ04_HK_genotype_ni)
colnames(NSQ04_HK_genotype_ni)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_KAS1_HK_calculation_ni)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_HK_ni_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_HK_ni_test <- DESeq(NSQ04_HK_ni_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_HK_NI_Ba <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [inf any]
rhizo_genotype_HK_NI_Ba_FDR005 <- rhizo_genotype_HK_NI_Ba[(rownames(rhizo_genotype_HK_NI_Ba)[which(rhizo_genotype_HK_NI_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_HK_NI_Ba_FDR005

#Repeat for all genotypes

#Morex
rhizo_genotype_HK_NI_M <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_HK_NI_M_FDR005 <- rhizo_genotype_HK_NI_M[(rownames(rhizo_genotype_HK_NI_M)[which(rhizo_genotype_HK_NI_M$padj <0.05)]), ]
rhizo_genotype_HK_NI_M_FDR005

#17
rhizo_genotype_HK_NI_17 <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_HK_NI_17_FDR005 <- rhizo_genotype_HK_NI_17[(rownames(rhizo_genotype_HK_NI_17)[which(rhizo_genotype_HK_NI_17$padj <0.05)]), ]
rhizo_genotype_HK_NI_17_FDR005

#52
rhizo_genotype_HK_NI_52 <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_HK_NI_52_FDR005 <- rhizo_genotype_HK_NI_52[(rownames(rhizo_genotype_HK_NI_52)[which(rhizo_genotype_HK_NI_52$padj <0.05)]), ]
rhizo_genotype_HK_NI_52_FDR005

#Barke vs other genotypes?

#Morex
rhizo_genotype_HK_NI_M <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_HK_NI_M_FDR005 <- rhizo_genotype_HK_NI_M[(rownames(rhizo_genotype_HK_NI_M)[which(rhizo_genotype_HK_NI_M$padj <0.05)]), ]
rhizo_genotype_HK_NI_M_FDR005

#17
rhizo_genotype_HK_NI_17 <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_HK_NI_17_FDR005 <- rhizo_genotype_HK_NI_17[(rownames(rhizo_genotype_HK_NI_17)[which(rhizo_genotype_HK_NI_17$padj <0.05)]), ]
rhizo_genotype_HK_NI_17_FDR005

#52
rhizo_genotype_HK_NI_52 <- results(NSQ04_HK_ni_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_HK_NI_52_FDR005 <- rhizo_genotype_HK_NI_52[(rownames(rhizo_genotype_HK_NI_52)[which(rhizo_genotype_HK_NI_52$padj <0.05)]), ]
rhizo_genotype_HK_NI_52_FDR005

#Repeat for just run 3

#subset for run 3
NSQ04_KAS1_HK_calculation_r3 <- subset_samples(NSQ04_KAS1_HK_calculation, Date.Inoc == "241031")

#extract count data 
NSQ04_HK_genotype_r3 <- otu_table(NSQ04_KAS1_HK_calculation_r3)
countData = as.data.frame(NSQ04_HK_genotype_r3)
colnames(NSQ04_HK_genotype_r3)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_KAS1_HK_calculation_r3)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_HK_r3_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_HK_r3_test <- DESeq(NSQ04_HK_r3_cds, fitType="local", betaPrior = FALSE) 

#No bulk

#Barke vs other genotypes?

#Morex
rhizo_genotype_HK_r3_M <- results(NSQ04_HK_r3_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_HK_r3_M_FDR005 <- rhizo_genotype_HK_r3_M[(rownames(rhizo_genotype_HK_r3_M)[which(rhizo_genotype_HK_r3_M$padj <0.05)]), ]
rhizo_genotype_HK_r3_M_FDR005

#17
rhizo_genotype_HK_r3_17 <- results(NSQ04_HK_r3_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_HK_r3_17_FDR005 <- rhizo_genotype_HK_r3_17[(rownames(rhizo_genotype_HK_r3_17)[which(rhizo_genotype_HK_r3_17$padj <0.05)]), ]
rhizo_genotype_HK_r3_17_FDR005

#52
rhizo_genotype_HK_r3_52 <- results(NSQ04_HK_r3_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_HK_r3_52_FDR005 <- rhizo_genotype_HK_r3_52[(rownames(rhizo_genotype_HK_r3_52)[which(rhizo_genotype_HK_r3_52$padj <0.05)]), ]
rhizo_genotype_HK_r3_52_FDR005

#Repeat for just run 4

#subset for run 4
NSQ04_KAS1_HK_calculation_r4 <- subset_samples(NSQ04_KAS1_HK_calculation, Date.Inoc == "241107")

#extract count data 
NSQ04_HK_genotype_r4 <- otu_table(NSQ04_KAS1_HK_calculation_r4)
countData = as.data.frame(NSQ04_HK_genotype_r4)
colnames(NSQ04_HK_genotype_r4)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(NSQ04_KAS1_HK_calculation_r4)))
rownames(colData)

#construct a DESeq dataset combining count data and sample information and specify the factor we will investigate
NSQ04_HK_r4_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Genotype)

#execute the differential count analysis with the function DESeq 
NSQ04_HK_r4_test <- DESeq(NSQ04_HK_r4_cds, fitType="local", betaPrior = FALSE) 

#define the Genera differentially enriched in Barke vs Bulk
rhizo_genotype_HK_r4_Ba <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Bulk", "Barke")) 

#extract  ASVs whose adjusted p.value in a given comparison is below 0.05 [inf any]
rhizo_genotype_HK_r4_Ba_FDR005 <- rhizo_genotype_HK_r4_Ba[(rownames(rhizo_genotype_HK_r4_Ba)[which(rhizo_genotype_HK_r4_Ba$padj <0.05)]), ]

#inspect the generated file
rhizo_genotype_HK_r4_Ba_FDR005

#Repeat for all genotypes

#Morex
rhizo_genotype_HK_r4_M <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Bulk", "Morex")) 
rhizo_genotype_HK_r4_M_FDR005 <- rhizo_genotype_HK_r4_M[(rownames(rhizo_genotype_HK_r4_M)[which(rhizo_genotype_HK_r4_M$padj <0.05)]), ]
rhizo_genotype_HK_r4_M_FDR005

#17
rhizo_genotype_HK_r4_17 <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Bulk", "124-17")) 
rhizo_genotype_HK_r4_17_FDR005 <- rhizo_genotype_HK_r4_17[(rownames(rhizo_genotype_HK_r4_17)[which(rhizo_genotype_HK_r4_17$padj <0.05)]), ]
rhizo_genotype_HK_r4_17_FDR005

#52
rhizo_genotype_HK_r4_52 <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Bulk", "124-52")) 
rhizo_genotype_HK_r4_52_FDR005 <- rhizo_genotype_HK_r4_52[(rownames(rhizo_genotype_HK_r4_52)[which(rhizo_genotype_HK_r4_52$padj <0.05)]), ]
rhizo_genotype_HK_r4_52_FDR005

#Barke vs other genotypes?

#Morex
rhizo_genotype_HK_r4_M <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Barke", "Morex")) 
rhizo_genotype_HK_r4_M_FDR005 <- rhizo_genotype_HK_r4_M[(rownames(rhizo_genotype_HK_r4_M)[which(rhizo_genotype_HK_r4_M$padj <0.05)]), ]
rhizo_genotype_HK_r4_M_FDR005

#17
rhizo_genotype_HK_r4_17 <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Barke", "124-17")) 
rhizo_genotype_HK_r4_17_FDR005 <- rhizo_genotype_HK_r4_17[(rownames(rhizo_genotype_HK_r4_17)[which(rhizo_genotype_HK_r4_17$padj <0.05)]), ]
rhizo_genotype_HK_r4_17_FDR005

#52
rhizo_genotype_HK_r4_52 <- results(NSQ04_HK_r4_test, contrast = c("Genotype", "Barke", "124-52")) 
rhizo_genotype_HK_r4_52_FDR005 <- rhizo_genotype_HK_r4_52[(rownames(rhizo_genotype_HK_r4_52)[which(rhizo_genotype_HK_r4_52$padj <0.05)]), ]
rhizo_genotype_HK_r4_52_FDR005
